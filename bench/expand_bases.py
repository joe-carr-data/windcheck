"""BASE MANIFEST for the EXPANSION corpus (the 99 traces beyond the pin).

The August target says "every published Herculaneum surface trace". The
round-28 base manifest (out/corpus_bases.json) pins the 185 segments of
the five original samples; bench/corpus_expand.py later censused the other
nine samples' 99 traces (out/expand/expand.jsonl, published as
results/corpus/expansion_records.jsonl). This writes the SAME kind of
manifest for those 99, so bench/excise_corpus.py can run the identical
frozen transform pass over them:

    uv run python bench/expand_bases.py
    uv run python bench/excise_corpus.py \
        --base-manifest out/expand_bases.json --out-root out/excised/expand

It is an ADDITIVE SECOND MANIFEST: out/corpus_bases.json is not read,
not restated and not touched, so the two corpora report separately or
together without either being recomputed.

Differences from the pinned manifest, all deliberate:

  * base_kind is "original" for every entry -- no displacement-repair
    campaign ever ran on these traces, so there is no repaired mesh to
    verify and no certificate to cite.
  * There is no published per-sample MANIFEST-*.json for these samples;
    the enumeration authority is the census record file itself, whose
    sha256 every entry carries (census_records / census_records_sha256).
  * voxel_um is the TRACE'S OWN claim, derived from its meta.json as
    sqrt(area_cm2 / area_vx2) * 1e4, recorded with the derivation named.
    Where meta.json carries no area fields the value is null -- NEVER a
    guessed default: the excision operates in voxel units and only
    physical-unit reporting depends on it. (Cross-check: for PHercMANBp
    and PHerc0332 the derived 2.399 um equals the voxel size in the
    sample's single published volume name.)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
for _p in (REPO_ROOT / "src", REPO_ROOT / "bench"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from windcheck.manifest import mesh_manifest                    # noqa: E402
from corpus_bases import (geometry_key, head_commit,            # noqa: E402
                          repo_path, sha_file)

SCHEMA = "corpus_bases/v1"
AXES = ("x", "y", "z")
RECORDS = REPO_ROOT / "out" / "expand" / "expand.jsonl"
DATA_ROOT = REPO_ROOT / "data" / "expand"
DEFAULT_OUT = REPO_ROOT / "out" / "expand_bases.json"

VOXEL_NOTE = ("derived from the trace's own meta.json as "
              "sqrt(area_cm2/area_vx2)*1e4; null when meta.json carries "
              "no area fields (physical-unit fields in the certificate "
              "are then null; the transform itself is voxel-space)")


def trace_voxel_um(mesh_dir: Path) -> float | None:
    meta = mesh_dir / "meta.json"
    if not meta.is_file():
        return None
    try:
        doc = json.loads(meta.read_text())
    except (OSError, ValueError):
        return None
    cm2, vx2 = doc.get("area_cm2"), doc.get("area_vx2")
    if not (isinstance(cm2, (int, float)) and isinstance(vx2, (int, float))
            and cm2 > 0 and vx2 > 0):
        return None
    return math.sqrt(cm2 / vx2) * 1e4


def build(records: Path = RECORDS, data_root: Path = DATA_ROOT) -> dict:
    rows = [json.loads(l) for l in records.read_text().splitlines() if l]
    names = [r["segment"] for r in rows]
    if len(set(names)) != len(names):
        raise SystemExit("segment names are not unique across samples; the "
                         "manifest keys on segment name")

    entries = []
    for r in sorted(rows, key=lambda r: r["segment"]):
        mesh = data_root / r["sample"] / r["segment"]
        if not mesh.is_dir():
            raise SystemExit(f"{r['sample']}/{r['segment']}: no local mesh "
                             f"at {mesh} -- refetch with corpus_expand.py")
        hashes = mesh_manifest(mesh)
        for ax in AXES:
            if hashes.get(ax) is None:
                raise SystemExit(f"{r['segment']}: missing {ax}.tif")
        planes = {a: hashes.get(a) for a in AXES}
        planes["mask"] = hashes.get("mask")
        entries.append({
            "segment": r["segment"],
            "corpus": r["sample"],
            "volume": r.get("volume") or None,
            "voxel_um": trace_voxel_um(mesh),
            "voxel_um_note": VOXEL_NOTE,
            "original_mesh": repo_path(mesh),
            "original_hashes": planes,
            "published_manifest_hashes": None,
            "published_manifest_source": None,
            "published_manifest_agrees_with_disk": None,
            "s3_prefix": r.get("s3_prefix"),
            "census_records": repo_path(records),
            "base_kind": "original",
            "base_mesh": repo_path(mesh),
            "base_hashes": planes,
            "repair_certificate": None,
            "repair_certificate_sha256": None,
            "repair_cert_output_hashes_verified": None,
            "repair_verification_note": ("expansion corpus: no displacement-"
                                         "repair campaign ran on these "
                                         "traces; base is the original "
                                         "fetched mesh"),
            "geometry_key": geometry_key(planes),
            "original_geometry_key": geometry_key(planes),
            "duplicate_of": None,
            "is_canonical": True,
        })

    # Duplicate aliasing, same rule as the pinned manifest: identical
    # coordinate planes alias to the lexicographically smallest name.
    by_key: dict[str, list[dict]] = {}
    for e in entries:
        by_key.setdefault(e["geometry_key"], []).append(e)
    n_dup = 0
    for group in by_key.values():
        if len(group) < 2:
            continue
        canon = min(group, key=lambda e: e["segment"])
        for e in group:
            if e is not canon:
                e["is_canonical"] = False
                e["duplicate_of"] = canon["segment"]
                n_dup += 1

    records_sha = sha_file(records)
    for e in entries:
        e["census_records_sha256"] = records_sha

    return {
        "schema": SCHEMA,
        "corpus_note": ("EXPANSION corpus base manifest: the nine samples "
                        "and 99 traces bench/corpus_expand.py added beyond "
                        "the 185 pinned segments. Additive to "
                        "out/corpus_bases.json, which is unchanged."),
        "generated_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "code_commit": head_commit(),
        "n_segments": len(entries),
        "n_canonical": sum(1 for e in entries if e["is_canonical"]),
        "n_duplicate_aliases": n_dup,
        "n_displacement_repaired": 0,
        "n_original": len(entries),
        "entries": entries,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--records", type=Path, default=RECORDS)
    args = ap.parse_args(argv)
    doc = build(args.records)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(doc, indent=1) + "\n")
    samples = sorted({e["corpus"] for e in doc["entries"]})
    n_vox = sum(1 for e in doc["entries"] if e["voxel_um"] is not None)
    print(f"{args.out}: {doc['n_segments']} entries "
          f"({doc['n_canonical']} canonical, "
          f"{doc['n_duplicate_aliases']} duplicate aliases) across "
          f"{len(samples)} samples; voxel_um known for {n_vox}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
