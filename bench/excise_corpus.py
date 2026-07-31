"""W4 corpus driver: certified excision over the whole pinned corpus.

Pure orchestration. bench/excise_shadow.py --certificate does the work and
its measurement logic, thresholds and FROZEN policy are never touched from
here. This script decides only WHICH segments run, enforces the ten-minute
per-segment process limit, and guarantees round-28 Q4(2):

    EVERY entry in the base manifest ends with a terminal record.

A crash, a timeout, a duplicate alias and a clean transformation are all
equally valid outputs. The only failure mode this driver has is leaving a
segment without a record, so a child that dies still gets a certificate
written for it -- by the driver, marked as driver-written, carrying the log
tail as evidence.

    uv run python bench/excise_corpus.py --jobs 3
    uv run python bench/excise_corpus.py --corpus "Scroll 1" --limit 5
    uv run python bench/excise_corpus.py --segments picks.txt --jobs 1
    uv run python bench/excise_corpus.py --shard 0/2   # machine A
    uv run python bench/excise_corpus.py --shard 1/2   # machine B

Bases: out/corpus_bases.json (schema corpus_bases/v1). The manifest, not a
glob, decides which bytes get cut -- that is what makes a certificate
replayable. DUPLICATE ALIASES (is_canonical false) are never executed: they
get an immediate terminal certificate pointing at their canonical geometry,
so the artifact count and the unique-geometry count can both be reported
without double-counting a cut that was only performed once.

Checkpointing: a segment is skipped when a valid certificate exists AND its
recorded base hashes match the manifest's base_hashes AND its policy_hash
matches the current frozen policy AND its source-tree digest matches
the current code --
unless --force. Driver-written failure records never checkpoint-skip.

Summary: <out-root>/corpus_summary.jsonl, one record per segment, rewritten
atomically after every completion (idempotent across resumes).
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

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT / "src", ROOT / "bench"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from windcheck.excise import (FROZEN_POLICY_VERSION,            # noqa: E402
                              frozen_policy_hash)

OUT = Path("out/excised/corpus")
BASE_MANIFEST = Path("out/corpus_bases.json")
DRIVER_CMD = "uv run python bench/excise_shadow.py"
AXES = ("x", "y", "z")
LOG_TAIL_LINES = 20

# Round-28 Q2: "hard 10-minute process limit". The driver enforces it as a
# wall-clock kill, so a segment that blows the budget is a RESULT (an
# `error`/`timeout` terminal record) and never a hung corpus pass.
SEGMENT_TIMEOUT_S = 600.0

CERT_SUFFIX = "_excision_certificate.json"
RECORD_KIND = "excision certificate"


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def code_identity() -> str:
    """The content digest of the published source tree. This is the driver's
    freshness key: it changes exactly when the code that produces a result
    changes, and unlike a commit sha a public consumer can recompute it."""
    import sys as _sys
    _src = str(Path(__file__).resolve().parents[1] / "src")
    if _src not in _sys.path:
        _sys.path.insert(0, _src)
    from windcheck.provenance import source_tree_digest
    return source_tree_digest()


def cert_code_identity(cert: dict) -> str | None:
    """A certificate's code identity, accepting the pre-release spelling."""
    for key in ("source_tree_digest", "code_identity"):
        v = cert.get(key)
        if isinstance(v, str) and v.strip():
            return v
    prov = cert.get("provenance")
    if isinstance(prov, dict) and prov.get("source_tree_digest"):
        return str(prov["source_tree_digest"])
    return None


def cert_path(out: Path, segment: str) -> Path:
    return Path(out) / f"{segment}{CERT_SUFFIX}"


# ---------------------------------------------------------------- manifest

def load_manifest(path: Path) -> dict:
    """The base manifest plus its own sha256 (a certificate must be able to
    name the exact manifest revision it was cut under)."""
    path = Path(path)
    doc = json.loads(path.read_text())
    entries = list(doc.get("entries") or [])
    entries.sort(key=lambda e: e["segment"])
    return {"path": str(path), "sha256": sha(path),
            "schema": doc.get("schema"),
            "provenance": doc.get("provenance"),
            "entries": entries}


def parse_shard(spec: str) -> tuple[int, int]:
    try:
        i, n = (int(x) for x in spec.split("/"))
    except ValueError:
        raise SystemExit(f"--shard must be i/N, got {spec!r}")
    if not (n >= 1 and 0 <= i < n):
        raise SystemExit(f"--shard out of range: {spec!r}")
    return i, n


def select_segments(entries: list[dict], corpus: str | None,
                    segments_file: Path | None, shard: str | None,
                    limit: int | None) -> list[dict]:
    """Filters applied in a FIXED order: corpus -> name list -> shard ->
    limit (the same order as repair_corpus.select_segments). Entries arrive
    and stay sorted by segment name, so shard i/N is index-mod-N over the
    sorted survivors: deterministic, disjoint, and the union over i covers
    everything."""
    rows = list(entries)
    if corpus:
        rows = [r for r in rows if corpus.lower() in (r.get("corpus") or
                                                      "").lower()]
    if segments_file:
        want = {ln.strip() for ln in Path(segments_file).read_text()
                .splitlines() if ln.strip() and not ln.startswith("#")}
        rows = [r for r in rows if r["segment"] in want]
        missing = want - {r["segment"] for r in rows}
        if missing:
            print(f"WARNING: {len(missing)} requested segment(s) absent from "
                  f"the base manifest: {sorted(missing)}", flush=True)
    if shard:
        i, n = parse_shard(shard)
        rows = [r for k, r in enumerate(rows) if k % n == i]
    if limit is not None:
        rows = rows[:limit]
    return rows


# -------------------------------------------------------------- checkpoint

# Dispositions that represent a finished measurement of the segment. An
# `error` record (crash, timeout, no feasible mask) is retried on the next
# pass rather than frozen in, because nothing about the geometry was
# established by it.
TERMINAL_OK_TO_SKIP = ("transformed", "already_clean", "duplicate_alias",
                       "triangle_empty_invalid", "not_censusable",
                       "residual_transverse")


def base_hashes_match(cert: dict, entry: dict) -> bool:
    """The certificate must record the SAME base bytes the manifest pins.
    A missing or partial hash record is a mismatch, never a pass."""
    got = cert.get("base_hashes") or {}
    want = entry.get("base_hashes") or {}
    if not want:
        return False
    for ax in AXES:
        if want.get(ax) is None or got.get(ax) != want.get(ax):
            return False
    return True


def checkpoint(entry: dict, commit: str, policy_hash: str,
               force: bool = False,
               out: Path = OUT) -> tuple[str, dict | None]:
    """('skip', cert) only when a valid certificate exists AND its base
    hashes match the manifest AND its policy_hash matches the current frozen
    policy AND its source-tree digest matches the current code.
    Otherwise ('run'/'rerun_*',
    None) naming which of the three moved."""
    cpath = cert_path(out, entry["segment"])
    if not cpath.exists():
        return "run", None
    if force:
        return "rerun_forced", None
    try:
        cert = json.loads(cpath.read_text())
    except (json.JSONDecodeError, OSError):
        return "rerun_unreadable_cert", None
    if cert.get("record_kind") != RECORD_KIND:
        return "rerun_unreadable_cert", None
    if cert.get("terminal_disposition") not in TERMINAL_OK_TO_SKIP:
        # a driver-written failure record is a record, not a result to keep
        return "rerun_previous_failure", None
    if not base_hashes_match(cert, entry):
        return "rerun_base_changed", None
    if cert.get("policy_hash") != policy_hash:
        return "rerun_policy_changed", None
    if cert_code_identity(cert) != commit:
        return "rerun_code_changed", None
    return "skip", cert


# ------------------------------------------------------- terminal records

def entry_fields(entry: dict) -> dict:
    return {"segment": entry["segment"], "corpus": entry.get("corpus"),
            "base_kind": entry.get("base_kind"),
            "base_mesh": entry.get("base_mesh"),
            "base_hashes": entry.get("base_hashes"),
            "original_mesh": entry.get("original_mesh"),
            "original_hashes": entry.get("original_hashes"),
            "repair_certificate": entry.get("repair_certificate"),
            "repair_certificate_sha256":
                entry.get("repair_certificate_sha256"),
            "geometry_key": entry.get("geometry_key"),
            "original_geometry_key": entry.get("original_geometry_key"),
            "duplicate_of": entry.get("duplicate_of"),
            "is_canonical": entry.get("is_canonical"),
            "voxel_um": entry.get("voxel_um")}


def write_terminal_certificate(out: Path, entry: dict, manifest: dict,
                               commit: str, policy_hash: str,
                               disposition: str, status: str,
                               wall_seconds: float, **extra) -> Path:
    """A terminal certificate the DRIVER writes: for duplicate aliases (no
    run is defined) and for children that died without leaving one. Round-28
    Q4(2) requires all 185 entries to carry a disposition; a missing file is
    not an acceptable answer for either case."""
    rec = {
        "record_kind": RECORD_KIND,
        "written_by": "bench/excise_corpus.py",
        "driver_terminal_record": True,
        "terminal_disposition": disposition,
        "status": status,
        **entry_fields(entry),
        "base_manifest": manifest["path"],
        "base_manifest_sha256": manifest["sha256"],
        "policy_version": FROZEN_POLICY_VERSION,
        "policy_hash": policy_hash,
        "source_tree_digest": commit,
        "operational_retained_fraction": None,
        "headline_retained_fraction": None,
        "core_gate_pass": None,
        "output_mesh": None,
        "output_mesh_hashes": None,
        "wall_seconds": round(wall_seconds, 2),
        "recorded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        **extra}
    p = cert_path(out, entry["segment"])
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(rec, indent=1, default=str))
    tmp.replace(p)
    return p


DUPLICATE_NOTE = (
    "byte-identical alias of an already-cut canonical base geometry. The "
    "excision is performed ONCE, on the canonical entry; this artifact "
    "inherits that result and must NOT be counted again in the "
    "unique-geometry-weighted corpus retention (round-28 Q3).")


def duplicate_record(out: Path, entry: dict, manifest: dict, commit: str,
                     policy_hash: str) -> dict:
    """Duplicate aliases are NOT run -- there is nothing to compute. They
    get their terminal certificate immediately."""
    p = write_terminal_certificate(
        out, entry, manifest, commit, policy_hash,
        disposition="duplicate_alias", status="duplicate_alias",
        wall_seconds=0.0, note=DUPLICATE_NOTE,
        canonical_segment=entry.get("duplicate_of"))
    return summary_record(entry, "duplicate_alias", "duplicate_alias",
                          p, 0.0, cert=json.loads(p.read_text()))


# ----------------------------------------------------------------- summary

def summary_record(entry: dict, disposition: str, status: str,
                   cpath: Path | None, wall_s: float,
                   cert: dict | None = None, **extra) -> dict:
    cert = cert or {}
    return {
        "segment": entry["segment"],
        "corpus": entry.get("corpus"),
        "terminal_disposition": disposition,
        "status": status,
        "operational_retained_fraction":
            cert.get("operational_retained_fraction"),
        "headline_retained_fraction": cert.get("headline_retained_fraction"),
        "core_gate_pass": cert.get("core_gate_pass"),
        "claimed_clean": cert.get("claimed_clean"),
        "wall_seconds": (cert.get("wall_seconds")
                         if cert.get("wall_seconds") is not None
                         else round(wall_s, 2)),
        "certificate": (str(cpath) if cpath else None),
        "certificate_sha256": (sha(cpath) if cpath and Path(cpath).exists()
                               else None),
        "output_mesh": cert.get("output_mesh"),
        "base_kind": entry.get("base_kind"),
        "geometry_key": entry.get("geometry_key"),
        "is_canonical": entry.get("is_canonical"),
        "duplicate_of": entry.get("duplicate_of"),
        **extra}


def load_summary(path: Path) -> dict[str, dict]:
    recs: dict[str, dict] = {}
    if Path(path).exists():
        for line in Path(path).read_text().splitlines():
            if line.strip():
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                recs[r["segment"]] = r
    return recs


def write_summary(path: Path, recs: dict[str, dict]) -> None:
    """Atomic rewrite, one record per segment, sorted -- idempotent across
    resumes (a rerun's record replaces the old one)."""
    path = Path(path)
    tmp = path.with_suffix(".jsonl.tmp")
    tmp.write_text("".join(json.dumps(recs[k]) + "\n" for k in sorted(recs)))
    tmp.replace(path)


# ----------------------------------------------------------------- running

def log_tail(path: Path, n: int = LOG_TAIL_LINES) -> str:
    try:
        return "\n".join(Path(path).read_text(errors="replace")
                         .splitlines()[-n:])
    except OSError:
        return ""


def run_segment(entry: dict, driver_cmd: list[str], timeout: float,
                out: Path, env: dict, manifest: dict, commit: str,
                policy_hash: str) -> dict:
    """One child = one segment, in its own session so a timeout kills the
    whole process group. Returns a summary record; a child that fails to
    leave a certificate gets one written for it here."""
    seg = entry["segment"]
    logdir = Path(out) / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    log = logdir / f"{seg}.log"
    cmd = list(driver_cmd) + ["--segment", seg, "--certificate",
                              "--out-root", str(out),
                              "--base-manifest", manifest["path"]]
    cpath = cert_path(out, seg)
    t0 = time.time()
    with open(log, "w") as lf:
        lf.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
        lf.flush()
        proc = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT,
                                start_new_session=True, env=env, cwd=str(ROOT))
        try:
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                proc.kill()
            proc.wait()
            wall = time.time() - t0
            # a partial record from the killed child is not a certificate
            if cpath.exists():
                cpath.replace(cpath.with_name(cpath.name + ".killed"))
            p = write_terminal_certificate(
                out, entry, manifest, commit, policy_hash,
                disposition="error", status="timeout", wall_seconds=wall,
                timeout_s=timeout, elapsed_seconds=round(wall, 2),
                note=("exceeded the frozen policy's hard per-segment process "
                      "limit and was killed; a timeout is a RESULT, not a "
                      "driver crash"),
                log=str(log), log_tail=log_tail(log))
            return summary_record(entry, "error", "timeout", p, wall,
                                  log=str(log), timeout_s=timeout)
    wall = time.time() - t0
    if rc != 0 or not cpath.exists():
        reason = ("child exited 0 but wrote no certificate" if rc == 0
                  else f"child exited {rc}")
        p = write_terminal_certificate(
            out, entry, manifest, commit, policy_hash,
            disposition="error", status="error", wall_seconds=wall,
            returncode=rc, reason=reason, log=str(log),
            log_tail=log_tail(log))
        return summary_record(entry, "error", "error", p, wall,
                              returncode=rc, reason=reason, log=str(log),
                              log_tail=log_tail(log))
    try:
        cert = json.loads(cpath.read_text())
    except (json.JSONDecodeError, OSError) as e:
        reason = f"unreadable certificate: {e}"
        p = write_terminal_certificate(
            out, entry, manifest, commit, policy_hash,
            disposition="error", status="error", wall_seconds=wall,
            reason=reason, log=str(log), log_tail=log_tail(log))
        return summary_record(entry, "error", "error", p, wall,
                              reason=reason, log=str(log))
    disp = cert.get("terminal_disposition")
    if disp is None:
        reason = "certificate carries no terminal_disposition"
        p = write_terminal_certificate(
            out, entry, manifest, commit, policy_hash,
            disposition="error", status="error", wall_seconds=wall,
            reason=reason, log=str(log), log_tail=log_tail(log))
        return summary_record(entry, "error", "error", p, wall,
                              reason=reason, log=str(log))
    return summary_record(entry, disp, cert.get("status") or disp, cpath,
                          wall, cert=cert, log=str(log))


def progress_line(rec: dict, done: int, total: int) -> str:
    if rec["terminal_disposition"] in ("transformed", "residual_transverse"):
        extra = (f"headline={rec.get('headline_retained_fraction')} "
                 f"operational={rec.get('operational_retained_fraction')} "
                 f"core_gate={rec.get('core_gate_pass')}")
    elif rec["terminal_disposition"] == "duplicate_alias":
        extra = f"alias_of={rec.get('duplicate_of')}"
    else:
        tail = (rec.get("log_tail") or "").splitlines()
        extra = rec.get("reason") or (tail[-1] if tail else "")
    return (f"[{done}/{total}] {rec['segment']} "
            f"{rec['terminal_disposition'].upper()} ({rec['status']}) "
            f"{extra} wall={rec.get('wall_seconds')}s")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="corpus driver for certified excision "
                    "(bench/excise_shadow.py --certificate)")
    ap.add_argument("--base-manifest", type=Path, default=BASE_MANIFEST,
                    help="round-28 Q3 base manifest (default "
                         f"{BASE_MANIFEST})")
    ap.add_argument("--out-root", type=Path, default=OUT,
                    help=f"certificates, meshes and logs (default {OUT})")
    ap.add_argument("--corpus", help="substring filter on corpus name")
    ap.add_argument("--segments", type=Path,
                    help="file with one segment name per line")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--shard", help="i/N deterministic shard by sorted "
                                    "segment name (multi-machine)")
    ap.add_argument("--jobs", type=int,
                    default=max(1, (os.cpu_count() or 4) // 4))
    ap.add_argument("--segment-timeout", type=float,
                    default=SEGMENT_TIMEOUT_S,
                    help="HARD per-segment wall-clock limit in seconds "
                         f"(default {SEGMENT_TIMEOUT_S:.0f}, the frozen "
                         "policy's ten-minute process gate)")
    ap.add_argument("--force", action="store_true",
                    help="rerun even when the checkpoint matches")
    ap.add_argument("--census-threads", type=int, default=None,
                    help="export WINDCHECK_CENSUS_THREADS to children "
                         "(cap: jobs * census-threads <= cores)")
    ap.add_argument("--driver-cmd", default=DRIVER_CMD,
                    help="command run per segment; '--segment <name> "
                         "--certificate --out-root <root> --base-manifest "
                         "<path>' is appended (tests override this)")
    args = ap.parse_args(argv)

    out = Path(args.out_root)
    out.mkdir(parents=True, exist_ok=True)
    summary_path = out / "corpus_summary.jsonl"
    commit = code_identity()
    policy_hash = frozen_policy_hash()
    driver_cmd = shlex.split(args.driver_cmd)
    env = dict(os.environ)
    if args.census_threads is not None:
        env["WINDCHECK_CENSUS_THREADS"] = str(args.census_threads)

    manifest = load_manifest(args.base_manifest)
    rows = select_segments(manifest["entries"], args.corpus, args.segments,
                           args.shard, args.limit)
    if not rows:
        print("no segments selected")
        return 1
    recs = load_summary(summary_path)

    todo: list[dict] = []
    done = 0
    n_dup = 0
    for entry in rows:
        if not entry.get("is_canonical", True):
            rec = duplicate_record(out, entry, manifest, commit, policy_hash)
            recs[rec["segment"]] = rec
            done += 1
            n_dup += 1
            print(f"[duplicate] ALIAS {entry['segment']} -> "
                  f"{entry.get('duplicate_of')} (not run)", flush=True)
            continue
        action, cert = checkpoint(entry, commit, policy_hash, args.force, out)
        if action == "skip":
            done += 1
            recs[entry["segment"]] = summary_record(
                entry, cert.get("terminal_disposition"),
                cert.get("status") or "skip",
                cert_path(out, entry["segment"]), 0.0, cert=cert,
                checkpoint="skip")
            print(f"[checkpoint] SKIP {entry['segment']} "
                  f"(cert + base hashes + policy_hash + source-tree digest match)",
                  flush=True)
        else:
            if action != "run":
                print(f"[checkpoint] {action.upper()} {entry['segment']}",
                      flush=True)
            todo.append(entry)
    write_summary(summary_path, recs)
    total = len(rows)
    print(f"{total} selected: {n_dup} duplicate alias(es), "
          f"{done - n_dup} checkpoint-skipped, {len(todo)} to run, "
          f"jobs={args.jobs}, timeout={args.segment_timeout:.0f}s, "
          f"policy {FROZEN_POLICY_VERSION} {policy_hash[:12]}, "
          f"commit {commit[:12]}", flush=True)

    counts: dict[str, int] = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {pool.submit(run_segment, entry, driver_cmd,
                            args.segment_timeout, out, env, manifest,
                            commit, policy_hash): entry for entry in todo}
        for fut in as_completed(futs):
            entry = futs[fut]
            try:
                rec = fut.result()
            except Exception as e:        # a driver bug is also a RESULT
                p = write_terminal_certificate(
                    out, entry, manifest, commit, policy_hash,
                    disposition="error", status="error", wall_seconds=0.0,
                    reason=f"driver: {e!r}")
                rec = summary_record(entry, "error", "error", p, 0.0,
                                     reason=f"driver: {e!r}")
            done += 1
            counts[rec["terminal_disposition"]] = counts.get(
                rec["terminal_disposition"], 0) + 1
            recs[rec["segment"]] = rec
            write_summary(summary_path, recs)
            print(progress_line(rec, done, total), flush=True)

    missing = [e["segment"] for e in rows if e["segment"] not in recs]
    if missing:
        print(f"BUG: {len(missing)} selected segment(s) left no record: "
              f"{missing[:10]}", flush=True)
    print(f"done: {json.dumps(counts)} (+{total - len(todo)} not run); "
          f"summary {summary_path}", flush=True)
    return 0 if not counts.get("error") and not missing else 2


if __name__ == "__main__":
    raise SystemExit(main())
