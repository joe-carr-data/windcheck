"""W3 corpus driver: run bench/repair_multi.py over the whole corpus,
safely in parallel, with checkpointing and a resumable summary.

Pure orchestration -- repair_multi.py's behavior is never modified. Each
segment runs as its own subprocess (transactional semantics mean a killed
child leaves only an inert partial workdir, never corrupted state).

    uv run python bench/repair_corpus.py --jobs 3
    uv run python bench/repair_corpus.py --corpus "Scroll 1" --limit 5
    uv run python bench/repair_corpus.py --segments picks.txt --jobs 1
    uv run python bench/repair_corpus.py --shard 0/2   # machine A
    uv run python bench/repair_corpus.py --shard 1/2   # machine B

Concurrency cap: each repair_multi census launches the selfcross engine
with threads=0 (= all cores). Use --census-threads N (exported to the
child as WINDCHECK_CENSUS_THREADS, read by crossing_census.census_one)
and keep jobs * census-threads <= physical cores on shared boxes.

Checkpointing: a segment is skipped when its certificate exists AND its
recorded hashes.input x/y/z match the current mesh files AND its
code_commit matches `git rev-parse HEAD` -- unless --force.

Summary: out/repaired/multi/corpus_summary.jsonl, one record per segment,
rewritten atomically after every completion (idempotent across resumes).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, "bench")
from repair_segment import CORPORA                              # noqa: E402

OUT = Path("out/repaired/multi")
DRIVER_CMD = "uv run python bench/repair_multi.py"
AXES = ("x", "y", "z")
LOG_TAIL_LINES = 15


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def head_commit() -> str:
    return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                          text=True).stdout.strip()


# ---------------------------------------------------------------- selection

def enumerate_segments(corpora=None) -> list[dict]:
    """(corpus, segment) rows for every segment directory that has a mesh
    matching the audited volume -- the same mesh-glob discipline as
    repair_multi.resolve_segment (mesh/*<volume>*.tifxyz; sorted()[0] is
    what the executor will pick). Sorted by segment name (the sharding
    key), deterministically."""
    rows = []
    for corpus, root, volume, _work in (CORPORA if corpora is None
                                        else corpora):
        rootp = Path(root)
        if not rootp.exists():
            continue
        for d in sorted(rootp.iterdir()):
            if not d.is_dir():
                continue
            meshes = sorted(d.glob(f"mesh/*{volume}*.tifxyz"))
            if not meshes:
                continue
            rows.append({"corpus": corpus, "segment": d.name,
                         "mesh": meshes[0]})
    rows.sort(key=lambda r: r["segment"])
    return rows


def parse_shard(spec: str) -> tuple[int, int]:
    try:
        i, n = (int(x) for x in spec.split("/"))
    except ValueError:
        raise SystemExit(f"--shard must be i/N, got {spec!r}")
    if not (n >= 1 and 0 <= i < n):
        raise SystemExit(f"--shard out of range: {spec!r}")
    return i, n


def select_segments(rows: list[dict], corpus: str | None,
                    segments_file: Path | None, shard: str | None,
                    limit: int | None) -> list[dict]:
    """Filters applied in a fixed order: corpus -> name list -> shard ->
    limit. Rows arrive (and stay) sorted by segment name, so shard i/N
    is index-mod-N over the sorted survivors: deterministic, disjoint,
    and the union over i covers everything."""
    if corpus:
        rows = [r for r in rows if corpus.lower() in r["corpus"].lower()]
    if segments_file:
        want = {ln.strip() for ln in Path(segments_file).read_text()
                .splitlines() if ln.strip() and not ln.startswith("#")}
        rows = [r for r in rows if r["segment"] in want]
        missing = want - {r["segment"] for r in rows}
        if missing:
            print(f"WARNING: {len(missing)} requested segment(s) not "
                  f"enumerable (no dir or no mesh): {sorted(missing)}",
                  flush=True)
    if shard:
        i, n = parse_shard(shard)
        rows = [r for k, r in enumerate(rows) if k % n == i]
    if limit is not None:
        rows = rows[:limit]
    return rows


# -------------------------------------------------------------- checkpoint

def checkpoint(row: dict, commit: str, force: bool = False,
               out: Path = OUT) -> tuple[str, dict | None]:
    """('skip', cert) when the certificate exists, its recorded input
    hashes match the current mesh files, and its code_commit matches the
    current HEAD; otherwise ('run'/'rerun_*', None)."""
    cpath = out / f"{row['segment']}_multi_certificate.json"
    if not cpath.exists():
        return "run", None
    if force:
        return "rerun_forced", None
    try:
        cert = json.loads(cpath.read_text())
    except (json.JSONDecodeError, OSError):
        return "rerun_unreadable_cert", None
    recorded = (cert.get("hashes") or {}).get("input") or {}
    for ax in AXES:
        f = row["mesh"] / f"{ax}.tif"
        if not f.exists() or recorded.get(f"{ax}.tif") != sha(f):
            return "rerun_input_hash_mismatch", None
    if cert.get("code_commit") != commit:
        return "rerun_code_commit_mismatch", None
    return "skip", cert


def cert_fields(cert: dict) -> dict:
    fin = cert.get("final_events") or {}
    instr = cert.get("instrumentation") or {}
    return {"class": cert.get("segment_class"),
            "repaired": len(cert.get("transactions") or []),
            "residual": {"d0": fin.get(0, fin.get("0")),
                         "d1": fin.get(1, fin.get("1"))},
            "wall_s": instr.get("wall_seconds"),
            "code_commit": cert.get("code_commit")}


# ----------------------------------------------------------------- summary

def load_summary(path: Path) -> dict[str, dict]:
    recs: dict[str, dict] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                recs[r["segment"]] = r
    return recs


def write_summary(path: Path, recs: dict[str, dict]) -> None:
    """Atomic rewrite, one record per segment, sorted -- idempotent
    across resumes (a rerun's record replaces the old one)."""
    tmp = path.with_suffix(".jsonl.tmp")
    tmp.write_text("".join(json.dumps(recs[k]) + "\n"
                           for k in sorted(recs)))
    tmp.replace(path)


# ----------------------------------------------------------------- running

def log_tail(path: Path, n: int = LOG_TAIL_LINES) -> str:
    try:
        return "\n".join(path.read_text(errors="replace")
                         .splitlines()[-n:])
    except OSError:
        return ""


def run_segment(row: dict, driver_cmd: list[str], timeout: float,
                out: Path, env: dict) -> dict:
    """One child = one segment. New session so a timeout can kill the
    whole process group; the segment's partial workdir is inert
    (repair_multi is transactional)."""
    logdir = out / "corpus_logs"
    logdir.mkdir(parents=True, exist_ok=True)
    log = logdir / f"{row['segment']}.log"
    cmd = driver_cmd + ["--segment", row["segment"]]
    t0 = time.time()
    with open(log, "w") as lf:
        proc = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT,
                                start_new_session=True, env=env)
        try:
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                proc.kill()
            proc.wait()
            return {"segment": row["segment"], "corpus": row["corpus"],
                    "status": "timeout", "timeout_s": timeout,
                    "wall_s": round(time.time() - t0, 1),
                    "log": str(log), "log_tail": log_tail(log)}
    wall = round(time.time() - t0, 1)
    base = {"segment": row["segment"], "corpus": row["corpus"],
            "wall_s": wall, "log": str(log)}
    if rc != 0:
        return {**base, "status": "error", "returncode": rc,
                "log_tail": log_tail(log)}
    cpath = out / f"{row['segment']}_multi_certificate.json"
    if not cpath.exists():
        return {**base, "status": "error", "returncode": rc,
                "log_tail": log_tail(log),
                "reason": "child exited 0 but wrote no certificate"}
    try:
        cert = json.loads(cpath.read_text())
    except (json.JSONDecodeError, OSError) as e:
        return {**base, "status": "error",
                "reason": f"unreadable certificate: {e}"}
    return {**base, "status": "ok", "cert": str(cpath), **cert_fields(cert)}


def progress_line(rec: dict, done: int, total: int) -> str:
    res = rec.get("residual") or {}
    if rec["status"] in ("ok", "skip"):
        extra = (f"class={rec.get('class')} repaired={rec.get('repaired')} "
                 f"residual=d0:{res.get('d0')},d1:{res.get('d1')}")
    else:
        tail = rec.get("log_tail", "").splitlines()
        extra = rec.get("reason") or (tail[-1] if tail else "")
    return (f"[{done}/{total}] {rec['segment']} "
            f"{rec['status'].upper()} {extra} wall={rec.get('wall_s')}s")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="parallel corpus driver for bench/repair_multi.py")
    ap.add_argument("--corpus", help="substring filter on corpus name")
    ap.add_argument("--segments", type=Path,
                    help="file with one segment name per line")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--shard", help="i/N deterministic shard by sorted "
                                    "segment name (multi-machine)")
    ap.add_argument("--jobs", type=int,
                    default=max(1, (os.cpu_count() or 4) // 4))
    ap.add_argument("--segment-timeout", type=float, default=3600.0,
                    help="per-segment wall-clock seconds (default 3600)")
    ap.add_argument("--force", action="store_true",
                    help="rerun even when the checkpoint matches")
    ap.add_argument("--census-threads", type=int, default=None,
                    help="export WINDCHECK_CENSUS_THREADS to children "
                         "(cap: jobs * census-threads <= cores)")
    ap.add_argument("--driver-cmd", default=DRIVER_CMD,
                    help="command run per segment; '--segment <name>' is "
                         "appended (tests override this)")
    ap.add_argument("--out-root", type=Path, default=OUT,
                    help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    out = args.out_root
    out.mkdir(parents=True, exist_ok=True)
    summary_path = out / "corpus_summary.jsonl"
    commit = head_commit()
    driver_cmd = shlex.split(args.driver_cmd)
    env = dict(os.environ)
    if args.census_threads is not None:
        env["WINDCHECK_CENSUS_THREADS"] = str(args.census_threads)

    rows = select_segments(enumerate_segments(), args.corpus, args.segments,
                           args.shard, args.limit)
    if not rows:
        print("no segments selected")
        return 1
    recs = load_summary(summary_path)

    todo, done = [], 0
    for row in rows:
        action, cert = checkpoint(row, commit, args.force, out)
        if action == "skip":
            done += 1
            rec = {"segment": row["segment"], "corpus": row["corpus"],
                   "status": "skip",
                   "cert": str(out / f"{row['segment']}"
                                     f"_multi_certificate.json"),
                   **cert_fields(cert)}
            # a checkpoint skip never downgrades an existing ok record
            if recs.get(row["segment"], {}).get("status") != "ok":
                recs[row["segment"]] = rec
            print(f"[checkpoint] SKIP {row['segment']} "
                  f"(cert + input hashes + code_commit match)", flush=True)
        else:
            if action.startswith("rerun"):
                print(f"[checkpoint] {action.upper()} {row['segment']}",
                      flush=True)
            todo.append(row)
    write_summary(summary_path, recs)
    total = len(rows)
    print(f"{total} selected: {done} checkpoint-skipped, {len(todo)} to "
          f"run, jobs={args.jobs}, timeout={args.segment_timeout:.0f}s, "
          f"commit {commit[:12]}", flush=True)

    counts: dict[str, int] = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {pool.submit(run_segment, row, driver_cmd,
                            args.segment_timeout, out, env): row
                for row in todo}
        for fut in as_completed(futs):
            row = futs[fut]
            try:
                rec = fut.result()
            except Exception as e:            # driver bug: record, continue
                rec = {"segment": row["segment"], "corpus": row["corpus"],
                       "status": "error", "reason": f"driver: {e!r}"}
            rec["code_commit"] = rec.get("code_commit") or commit
            done += 1
            counts[rec["status"]] = counts.get(rec["status"], 0) + 1
            recs[rec["segment"]] = rec
            write_summary(summary_path, recs)
            print(progress_line(rec, done, total), flush=True)

    print(f"done: {json.dumps(counts)} (+{total - len(todo)} skipped); "
          f"summary {summary_path}", flush=True)
    return 0 if counts.get("error", 0) + counts.get("timeout", 0) == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
