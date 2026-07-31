"""Experiment ledger: every attempt recorded and recoverable.

An ATTEMPT is any run whose outcome we might want to defend, repeat, or
undo: a pilot, a corpus pass, a kernel change validated against goldens,
a cloud run. Two commands:

    uv run python bench/ledger.py open  --id W1-shadow-monster \
        --goal "measure global-excision area cost on the monster" \
        --cmd "uv run python bench/excise_shadow.py --segment ..." \
        [--inputs path ...]

    uv run python bench/ledger.py close --id W1-shadow-monster \
        --verdict kept|reverted|inconclusive \
        --result "one line of what happened" \
        [--outputs path ...] [--notes "..."]

`open` refuses on a dirty tree (an attempt must be attributable to a
commit), records HEAD + a git tag `attempt/<id>`, hashes the declared
inputs, and appends a JSON record to out/ledger/attempts.jsonl plus a
human row in notes/LEDGER.md. `close` hashes the declared outputs,
records the verdict and the exact restore recipe.

Recovery for any attempt: `git checkout attempt/<id>` reproduces the
code; the record carries input hashes (so you know the data was the
same), output hashes (so you can tell whether a rerun matched), and the
command. Artifacts under out/ are gitignored, so `close` also copies
small declared outputs into out/ledger/<id>/ as a snapshot.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

LEDGER_DIR = Path("out/ledger")
JSONL = LEDGER_DIR / "attempts.jsonl"
MD = Path("notes/LEDGER.md")
SNAPSHOT_MAX_MB = 32.0


def git(*a: str) -> str:
    return subprocess.run(["git", *a], capture_output=True,
                          text=True).stdout.strip()


def digest(p: Path) -> str:
    """sha256 of a file, or of a directory's sorted (name, content)."""
    h = hashlib.sha256()
    if p.is_dir():
        for f in sorted(p.rglob("*")):
            if f.is_file():
                h.update(str(f.relative_to(p)).encode())
                h.update(f.read_bytes())
    elif p.exists():
        h.update(p.read_bytes())
    else:
        return "MISSING"
    return h.hexdigest()


def hash_all(paths) -> dict:
    return {str(p): digest(Path(p)) for p in paths or []}


def records() -> list[dict]:
    if not JSONL.exists():
        return []
    return [json.loads(x) for x in JSONL.read_text().splitlines() if x]


def write(rec: dict) -> None:
    LEDGER_DIR.mkdir(parents=True, exist_ok=True)
    with JSONL.open("a") as f:
        f.write(json.dumps(rec) + "\n")


def md_row(rec: dict) -> None:
    MD.parent.mkdir(parents=True, exist_ok=True)
    if not MD.exists():
        MD.write_text(
            "# Experiment ledger\n\n"
            "Every attempt from 2026-07-30 onward. Machine-readable\n"
            "records with hashes: `out/ledger/attempts.jsonl`. Code for\n"
            "any attempt: `git checkout attempt/<id>`. Written by\n"
            "`bench/ledger.py` — do not hand-edit rows.\n\n"
            "| id | opened | commit | goal | verdict | result |\n"
            "|---|---|---|---|---|---|\n")
    row = (f"| `{rec['id']}` | {rec['opened'][:16]} | `{rec['commit'][:8]}` "
           f"| {rec['goal']} | {rec.get('verdict', '_open_')} "
           f"| {rec.get('result', '')} |\n")
    MD.write_text(MD.read_text() + row)


def cmd_open(a) -> int:
    if git("status", "--porcelain"):
        print("REFUSED: working tree is dirty — commit first so the "
              "attempt is attributable to an exact code state.")
        return 1
    if any(r["id"] == a.id for r in records()):
        print(f"REFUSED: attempt id {a.id} already exists.")
        return 1
    commit = git("rev-parse", "HEAD")
    subprocess.run(["git", "tag", "-f", f"attempt/{a.id}", commit],
                   capture_output=True)
    rec = {"id": a.id, "opened": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                               time.gmtime()),
           "commit": commit, "tag": f"attempt/{a.id}", "goal": a.goal,
           "command": a.cmd, "inputs": hash_all(a.inputs),
           "uv_lock_sha256": digest(Path("uv.lock")),
           "engine_selfcross_sha256": digest(Path("engines/selfcross")),
           "state": "open"}
    write(rec)
    md_row(rec)
    print(f"OPENED {a.id} at {commit[:8]} (tag attempt/{a.id})\n  goal: "
          f"{a.goal}\n  cmd:  {a.cmd}")
    return 0


def cmd_close(a) -> int:
    recs = records()
    opens = [r for r in recs if r["id"] == a.id and r.get("state") == "open"]
    if not opens:
        print(f"REFUSED: no open attempt with id {a.id}.")
        return 1
    rec = dict(opens[-1])
    snap = LEDGER_DIR / a.id
    kept = []
    for p in a.outputs or []:
        src = Path(p)
        if not src.exists():
            continue
        mb = sum(f.stat().st_size for f in src.rglob("*") if f.is_file()) / 1e6 \
            if src.is_dir() else src.stat().st_size / 1e6
        if mb <= SNAPSHOT_MAX_MB:
            snap.mkdir(parents=True, exist_ok=True)
            dst = snap / src.name
            if dst.exists():
                shutil.rmtree(dst) if dst.is_dir() else dst.unlink()
            shutil.copytree(src, dst) if src.is_dir() else shutil.copy(src, dst)
            kept.append(str(dst))
    rec.update({
        "state": "closed",
        "closed": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "verdict": a.verdict, "result": a.result, "notes": a.notes or "",
        "outputs": hash_all(a.outputs),
        "snapshots": kept,
        "commit_at_close": git("rev-parse", "HEAD"),
        "recovery": (f"git checkout attempt/{a.id}  # code as run; inputs "
                     f"and outputs verified by the hashes in this record; "
                     f"rerun: {rec['command']}"),
    })
    write(rec)
    md_row(rec)
    print(f"CLOSED {a.id}: {a.verdict} — {a.result}")
    if kept:
        print("  snapshots: " + ", ".join(kept))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="sub", required=True)
    o = sub.add_parser("open")
    o.add_argument("--id", required=True)
    o.add_argument("--goal", required=True)
    o.add_argument("--cmd", required=True)
    o.add_argument("--inputs", nargs="*")
    c = sub.add_parser("close")
    c.add_argument("--id", required=True)
    c.add_argument("--verdict", required=True,
                   choices=["kept", "reverted", "inconclusive", "failed"])
    c.add_argument("--result", required=True)
    c.add_argument("--outputs", nargs="*")
    c.add_argument("--notes")
    a = ap.parse_args()
    return cmd_open(a) if a.sub == "open" else cmd_close(a)


if __name__ == "__main__":
    raise SystemExit(main())
