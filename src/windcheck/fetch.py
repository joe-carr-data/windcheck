"""Fetch tifxyz surface maps and pin them in a manifest.

Nothing under `data/` is committed, so `data/MANIFEST.json` is the provenance
record: exact S3 keys, byte counts and SHA-256 for every file a result depends
on. Any downstream number can be traced back to pinned bytes, and a corrupted
or silently-updated upstream file shows up as a hash mismatch rather than as a
mysteriously different result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from . import S3_BUCKET

# Only the two volume-space variants. The `intermediate/` flattened and
# normalized grids are outputs of the flattening step, not the traced geometry,
# so they are not part of the invariant check.
INCLUDE = "*/mesh/*.tifxyz/*"


def download(sample: str, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["aws", "s3", "cp", f"s3://{S3_BUCKET}/{sample}/segments/", str(dest),
         "--no-sign-request", "--recursive", "--exclude", "*", "--include", INCLUDE,
         "--only-show-errors"],
        check=True,
    )


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def write_manifest(sample: str, root: Path, out: Path) -> dict:
    """Hash every downloaded file and write the manifest."""
    files = sorted(p for p in root.rglob("*") if p.is_file())
    entries = []
    total = 0
    for p in files:
        rel = p.relative_to(root).as_posix()
        size = p.stat().st_size
        total += size
        entries.append({
            "path": rel,
            "s3_key": f"{sample}/segments/{rel}",
            "bytes": size,
            "sha256": sha256(p),
        })

    manifest = {
        "bucket": S3_BUCKET,
        "sample": sample,
        "include_pattern": INCLUDE,
        "n_files": len(entries),
        "total_bytes": total,
        "files": entries,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=1) + "\n")
    return manifest


def verify(root: Path, manifest_path: Path) -> list[str]:
    """Re-hash and return a list of problems (empty list == clean)."""
    manifest = json.loads(manifest_path.read_text())
    problems = []
    for e in manifest["files"]:
        p = root / e["path"]
        if not p.exists():
            problems.append(f"missing: {e['path']}")
        elif p.stat().st_size != e["bytes"]:
            problems.append(f"size mismatch: {e['path']}")
        elif sha256(p) != e["sha256"]:
            problems.append(f"hash mismatch: {e['path']}")
    return problems


def main(argv: list[str] | None = None) -> int:
    root_dir = Path(__file__).resolve().parents[2]
    p = argparse.ArgumentParser(prog="windcheck.fetch")
    p.add_argument("--sample", default="PHerc0172")
    p.add_argument("--skip-download", action="store_true")
    p.add_argument("--verify", action="store_true", help="check against an existing manifest")
    a = p.parse_args(argv)

    data = root_dir / "data" / "scroll5_tifxyz"
    man = root_dir / "data" / "MANIFEST.json"

    if a.verify:
        problems = verify(data, man)
        print("\n".join(problems) if problems else "MANIFEST VERIFIED: all files match")
        return 1 if problems else 0

    if not a.skip_download:
        download(a.sample, data)
    m = write_manifest(a.sample, data, man)
    print(f"manifest: {m['n_files']} files, {m['total_bytes'] / 1e6:.1f} MB -> {man}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
