"""Audit a release archive's MEMBER LIST before it is uploaded.

A tarball built on macOS silently carried 1,574 AppleDouble `._*`
members and a second top-level entry into a release candidate. Nothing
inspected the member list, so nothing caught it: the payload was right
and the container was not.

Refuses on: more than one top-level entry, any `._*` or `.DS_Store`,
symlinks, hardlinks, absolute paths, or `..` traversal.
"""
from __future__ import annotations
import argparse, hashlib, json, tarfile
from pathlib import Path

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True)
    ap.add_argument("--expect-top", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    p = Path(a.archive)
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    problems, tops, n = [], set(), 0
    with tarfile.open(p, "r:gz") as tf:
        for m in tf.getmembers():
            n += 1
            parts = Path(m.name).parts
            tops.add(parts[0] if parts else m.name)
            base = Path(m.name).name
            if base.startswith("._"):
                problems.append(f"AppleDouble member: {m.name}")
            if base == ".DS_Store":
                problems.append(f".DS_Store member: {m.name}")
            if m.issym() or m.islnk():
                problems.append(f"link member: {m.name}")
            if m.name.startswith("/"):
                problems.append(f"absolute path: {m.name}")
            if ".." in parts:
                problems.append(f"parent traversal: {m.name}")
    if len(tops) != 1:
        problems.append(f"expected exactly one top-level entry, got {sorted(tops)}")
    if a.expect_top and tops != {a.expect_top}:
        problems.append(f"top-level entry is {sorted(tops)}, expected {a.expect_top!r}")
    out = {"archive": str(p), "sha256": h.hexdigest(), "bytes": p.stat().st_size,
           "n_members": n, "top_level": sorted(tops),
           "n_problems": len(problems), "problems": problems[:20],
           "verdict": "CLEAN" if not problems else "REFUSED"}
    if a.out:
        Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in
                      ("sha256", "bytes", "n_members", "top_level",
                       "n_problems", "verdict")}, indent=1))
    for x in problems[:10]:
        print("  ", x)
    return 0 if not problems else 1
if __name__ == "__main__":
    raise SystemExit(main())
