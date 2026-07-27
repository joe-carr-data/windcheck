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

# Only the volume-space meshes. The `intermediate/` flattened and normalized
# grids are outputs of the flattening step, not the traced geometry, so they are
# not part of the invariant check.
#
# The pattern is narrowed per corpus to the ONE published volume the audit
# measures. Without that, a scroll published at four resolutions pulls all four
# -- for Scroll 1 that is 133 MB per channel per segment against 2 MB for the
# volume actually used. Narrowing also makes the manifest describe exactly the
# files a result depends on, rather than a superset.
INCLUDE = "*/mesh/*{volume}*.tifxyz/*"

# Local directory for each published sample. All four share the S3 layout
# `<sample>/segments/`, verified against the bucket.
#
# The headline is measured across all four, so all four have to be pinned: a
# result that rests on corpora nobody else can check is not reproducible, it is
# just asserted. PHerc0172 keeps the unsuffixed manifest name because the
# already-published tool refers to it.
# sample -> (local directory, the volume id the audit measures)
CORPORA = {
    "PHercParis4": ("scroll1_tifxyz", "20230205180739"),
    "PHerc0172": ("scroll5_tifxyz", "20241024131839"),
    "PHerc0814": ("PHerc0814_tifxyz", "20250804134230"),
    "PHerc0139": ("PHerc0139_tifxyz", "20250728140407"),
    "PHerc1667": ("PHerc1667_tifxyz", "20231117161658"),
}


def paths_for(sample: str, root_dir: Path) -> tuple[Path, Path]:
    """(data directory, manifest path) for a sample."""
    if sample not in CORPORA:
        raise SystemExit(f"unknown sample {sample!r}; known: "
                         f"{', '.join(sorted(CORPORA))}")
    data = root_dir / "data" / CORPORA[sample][0]
    name = "MANIFEST.json" if sample == "PHerc0172" else f"MANIFEST-{sample}.json"
    return data, root_dir / "data" / name


def volume_of(sample: str) -> str:
    return CORPORA[sample][1]


def include_for(sample: str) -> str:
    return INCLUDE.format(volume=CORPORA[sample][1])


def download(sample: str, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["aws", "s3", "cp", f"s3://{S3_BUCKET}/{sample}/segments/", str(dest),
         "--no-sign-request", "--recursive", "--exclude", "*",
         "--include", include_for(sample), "--only-show-errors"],
        check=True,
    )


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def write_manifest(sample: str, root: Path, out: Path) -> dict:
    """Hash the files this corpus's result depends on, and write the manifest.

    Restricted to the audited volume, so the manifest is exactly the set
    `download` fetches -- a manifest listing files the fetcher would not bring
    down is not a reproduction recipe.
    """
    vol = CORPORA[sample][1]
    files = sorted(p for p in root.rglob("*")
                   if p.is_file() and "/mesh/" in p.as_posix()
                   and vol in p.as_posix() and ".tifxyz/" in p.as_posix())
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
        "include_pattern": include_for(sample),
        "volume": vol,
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
    p.add_argument("--sample", default="PHerc0172",
                   help=f"one of {', '.join(sorted(CORPORA))}, or 'all'")
    p.add_argument("--skip-download", action="store_true")
    p.add_argument("--verify", action="store_true", help="check against an existing manifest")
    a = p.parse_args(argv)

    samples = sorted(CORPORA) if a.sample == "all" else [a.sample]
    failed = 0
    for sample in samples:
        data, man = paths_for(sample, root_dir)

        if a.verify:
            if not man.exists():
                print(f"{sample}: NO MANIFEST at {man}")
                failed += 1
                continue
            problems = verify(data, man)
            if problems:
                failed += 1
                print(f"{sample}: {len(problems)} problem(s)")
                print("\n".join(problems[:20]))
            else:
                print(f"{sample}: VERIFIED, all files match")
            continue

        if not a.skip_download:
            download(sample, data)
        m = write_manifest(sample, data, man)
        print(f"{sample}: {m['n_files']} files, "
              f"{m['total_bytes'] / 1e9:.2f} GB -> {man}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
