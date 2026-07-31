"""Exact nonlocal self-intersection census of the published corpus.

The first measurement of whether published scroll meshes actually cross
themselves, as opposed to merely coming close. Proximity is ambiguous because
wraps in a crushed scroll genuinely lie microns apart. A transverse crossing is
not ambiguous: an embedded surface cannot pass through itself at any packing
density, which the test suite pins with folds 0.5 vx apart reporting nothing.

Run on everything, at both quad diagonals, and report per segment. A crossing
that appears under one diagonal and not the other is tessellation-sensitive and
is reported separately rather than counted.

    uv run python bench/crossing_census.py --root data/scroll5_tifxyz
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np

from windcheck import atlas, tifxyz

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "engines" / "selfcross"
# Per-run output directory. This was hardcoded to out/crossing/, so successive
# runs over different corpora wrote their CSVs into the SAME place and later
# analyses read a mixture of four scrolls. It produced a labelled-segment
# separation of 7,238 cells (~1,100 mm) that no segment in that corpus can
# support. Always keep runs separate.
WORK = Path("out/crossing")


_MESH_RE = re.compile(r"-on-(\d+)-(\d+\.?\d*)um\.tifxyz$")


def _volume_of(name: str) -> str | None:
    m = _MESH_RE.search(name)
    return m.group(1) if m else None


def _resolution_of(name: str) -> float | None:
    m = _MESH_RE.search(name)
    return float(m.group(2)) if m else None


def census_one(path: Path, name: str, exclude: int, cell: float,
               threads: int, maxedge: float, work: Path = WORK) -> dict | None:
    threads = int(os.environ.get("WINDCHECK_CENSUS_THREADS", threads))
    work.mkdir(parents=True, exist_ok=True)
    surf = tifxyz.read(path)
    if surf.valid.sum() < 5000:
        return None
    abin = work / "c_atlas.bin"
    atlas.write_atlas([_E(path)], abin)

    # Record the parameters AND the exact mesh with the result. A summary
    # without the parameters cannot be reproduced; a summary without the mesh
    # identity cannot be safely joined to. The CSVs below hold (v, u) indices
    # into *this* grid, and a later analysis that globs a different resolution
    # of the same segment reinterprets them against the wrong one. That happened
    # on PHerc1667 across all 20 segments and produced a plausible, entirely
    # fictional result -- so the mesh travels with the numbers now.
    row: dict = {"segment": name, "grid": list(surf.shape),
                 "valid_cells": int(surf.valid.sum()),
                 "mesh": {"path": str(path), "name": path.name,
                          "volume": _volume_of(path.name),
                          "resolution_um": _resolution_of(path.name)},
                 "params": {"exclude": exclude, "cell": cell,
                            "maxedge": maxedge, "touch_tol": 1e-3}}
    for diag in (0, 1):
        out = work / f"{name[:40]}_d{diag}.csv"
        r = subprocess.run(
            [str(ENGINE), str(abin), str(out), str(threads), str(cell),
             str(exclude), str(diag), str(maxedge)],
            capture_output=True, text=True, check=True)
        j = json.loads(r.stdout.strip().splitlines()[-1])
        row[f"d{diag}"] = {k: j[k] for k in
                           ("triangles", "quads_dropped", "pairs_tested",
                            "transverse", "coplanar", "grazing")}
        row[f"csv_d{diag}"] = str(out)
    # NOT a setwise comparison -- kept only as a coarse screen. Equal counts do
    # not mean the same regions crossed. crossing_analyse.py does it properly.
    row["transverse_min_of_diagonals"] = min(row["d0"]["transverse"],
                                             row["d1"]["transverse"])
    return row


class _E:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.winding = None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/scroll5_tifxyz"))
    ap.add_argument("--volume", default="20241024131839",
                    help="substring of the volume id / resolution to select")
    ap.add_argument("--exclude", type=int, default=1,
                    help="quads within this Chebyshev grid distance share a\n                         vertex; 1 is shared-vertex exclusion")
    ap.add_argument("--cell", type=float, default=40.0)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--maxedge", type=float, default=60.0,
                    help="drop quads with any edge longer than this (voxels); "
                         "grid pitch is ~20, discontinuities reach 1848")
    ap.add_argument("--json", type=Path, default=Path("out/crossing/census.json"))
    ap.add_argument("--work", type=Path, default=None,
                    help="per-run CSV directory; defaults to the --json parent "
                         "so separate corpora cannot contaminate each other")
    a = ap.parse_args()

    rows = []
    segs = sorted(d for d in a.root.iterdir() if d.is_dir())
    print(f"{'segment':46s} {'tris':>9s} {'d0 tvs':>7s} {'d1 tvs':>7s} "
          f"{'copl':>7s} {'graze':>7s}")
    print("-" * 92)
    for d in segs:
        m = sorted(d.glob(f"mesh/*{a.volume}*.tifxyz"))
        if not m:
            continue
        try:
            row = census_one(m[0], d.name, a.exclude, a.cell, a.threads,
                             a.maxedge, work=(a.work or a.json.parent))
        except subprocess.CalledProcessError as e:
            print(f"{d.name[:44]:46s}  ENGINE FAILED: {e.stderr.strip()[:40]}")
            continue
        if row is None:
            continue
        rows.append(row)
        print(f"{d.name[:44]:46s} {row['d0']['triangles']:9,d} "
              f"{row['d0']['transverse']:7,d} {row['d1']['transverse']:7,d} "
              f"{row['d0']['coplanar']:7,d} {row['d0']['grazing']:7,d}",
              flush=True)

    if not rows:
        print("no segments analysed")
        return

    tv0 = np.array([r["d0"]["transverse"] for r in rows])
    tv1 = np.array([r["d1"]["transverse"] for r in rows])
    both = np.array([r["transverse_min_of_diagonals"] for r in rows])
    print("-" * 92)
    print(f"{len(rows)} segments analysed, "
          f"{int((both > 0).sum())} with transverse crossings under BOTH diagonals")
    print(f"  total transverse: diagonal 0 = {tv0.sum():,}, "
          f"diagonal 1 = {tv1.sum():,}")
    print(f"  coplanar total: {sum(r['d0']['coplanar'] for r in rows):,}")
    print(f"  grazing total:  {sum(r['d0']['grazing'] for r in rows):,}")
    # Diagonal agreement is a SETWISE question and is answered in
    # crossing_analyse.py. Equal counts do not mean the same regions crossed.

    a.json.parent.mkdir(parents=True, exist_ok=True)
    a.json.write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {a.json}")

    # Deliberately no verdict line here. Two reasons it used to be wrong:
    # `min(d0, d1)` does not show the same regions crossed under both
    # tessellations, and triangle-pair counts are not defect counts. Run
    # bench/crossing_analyse.py, which clusters pairs into regions and events,
    # matches diagonals setwise, and reports the separation distribution instead
    # of a count past one cutoff.
    print("\nRaw pair counts above are NOT the finding. Run:")
    print("    uv run python bench/crossing_analyse.py")


if __name__ == "__main__":
    main()
