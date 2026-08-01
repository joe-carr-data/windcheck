"""Diagnose the patches the sweep flagged, and test the proposed remedy.

For each flagged patch: where do the crossings sit relative to the patch
boundary, and does eroding the valid mask remove them? The second half
matters more than the first. The dataset's author already suspected the
boundary ("the very last vertex on the boundary sometimes curls a bit
inwards or outwards; you can do a boundary erosion to take care of
this"), so the useful contribution is not repeating the hypothesis but
measuring it and running the fix.

Erosion is applied to the validity mask only. No coordinate is moved, so
a patch that comes back clean is clean because cells were dropped, not
because geometry was massaged -- the same discipline as the excision
operator on the pinned corpus.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import tifffile  # noqa: E402
from scipy import ndimage  # noqa: E402

from windcheck import pipeline, tifxyz  # noqa: E402

BASE = ("https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/"
        "verified_patches")


def census(mesh: Path, work: Path, tag: str) -> tuple[int, list[dict]]:
    """Both-diagonal transverse count plus the participating quads."""
    total, rows = 0, []
    for d in (0, 1):
        csvp, counts = pipeline.run_engine(mesh, tag, work, d, {"threads": 1})
        total += int(counts["transverse"])
        with csvp.open() as fh:
            for r in csv.DictReader(fh):
                if r["verdict"] == "transverse":
                    rows.append({"diagonal": d, **r})
    return total, rows


def diagnose(name: str, data: Path, work: Path, max_erosion: int) -> dict:
    mesh = data / f"{name}"
    s = tifxyz.read(mesh)
    base_total, rows = census(mesh, work, name[:36])

    # Depth of every valid cell, so "near the boundary" is measured against
    # how deep this patch actually is rather than asserted.
    depth = ndimage.distance_transform_edt(s.valid)
    inside = depth[s.valid]
    quads = sorted({(int(r[k[0]]), int(r[k[1]]))
                    for r in rows for k in (("v1", "u1"), ("v2", "u2"))})
    quad_depth = [float(min(depth[v, u], depth[v, u + 1],
                            depth[v + 1, u], depth[v + 1, u + 1]))
                  for v, u in quads]

    out = {
        "patch": name,
        "source": f"{BASE}/{name}",
        "grid": [int(x) for x in s.shape],
        "n_valid": s.n_valid,
        "transverse_both_diagonals": base_total,
        "participating_quads": [list(q) for q in quads],
        "quad_depth_cells": quad_depth,
        "median_depth_cells": float(np.median(inside)),
        "p90_depth_cells": float(np.percentile(inside, 90)),
        "max_penetration_vx": max(float(r["penetration"]) for r in rows),
        "grid_separation": [[abs(int(r["v1"]) - int(r["v2"])),
                             abs(int(r["u1"]) - int(r["u2"]))] for r in rows],
        "erosion": [],
    }

    for k in range(1, max_erosion + 1):
        V = ndimage.binary_erosion(s.valid, np.ones((3, 3), bool), iterations=k)
        dst = data / f"_eroded{k}_{name[:26]}.tifxyz"
        if dst.exists():
            shutil.rmtree(dst)
        dst.mkdir(parents=True)
        for i, ax in enumerate("xyz"):
            a = np.asarray(s.points[..., i], np.float32).copy()
            a[~V] = -1.0
            tifffile.imwrite(dst / f"{ax}.tif", a)
        tifxyz.write_meta(mesh, dst, s.points, V)
        tot, _ = census(dst, work, f"er{k}_{name[:30]}")
        out["erosion"].append({
            "cells": k, "transverse_both_diagonals": tot,
            "valid_retained": int(V.sum()),
            "fraction_retained": float(V.sum() / max(s.valid.sum(), 1)),
            "clean": tot == 0})
        shutil.rmtree(dst, ignore_errors=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit", default="out/patches/audit_full.jsonl")
    ap.add_argument("--data", default="out/patches/flagged")
    ap.add_argument("--work", default="out/patches/flagged/work")
    ap.add_argument("--out", default="out/patches/flagged/diagnosis.json")
    ap.add_argument("--max-erosion", type=int, default=3)
    a = ap.parse_args()

    flagged = [json.loads(l) for l in Path(a.audit).read_text().splitlines()
               if l.strip()]
    flagged = [r for r in flagged if r.get("status") == "ok" and not r["clean"]]
    print(f"{len(flagged)} flagged patches")

    data, work = Path(a.data), Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    results = []
    for r in flagged:
        name = r["patch"]
        if not (data / name).is_dir():
            print(f"  fetching {name}")
            import urllib.request
            (data / name).mkdir(parents=True, exist_ok=True)
            for f in ("meta.json", "x.tif", "y.tif", "z.tif"):
                with urllib.request.urlopen(f"{BASE}/{name}/{f}",
                                            timeout=120) as resp:
                    (data / name / f).write_bytes(resp.read())
        d = diagnose(name, data, work, a.max_erosion)
        d["sha256"] = r.get("sha256")
        results.append(d)
        first_clean = next((e["cells"] for e in d["erosion"] if e["clean"]),
                           None)
        print(f"  {name}")
        print(f"    transverse {d['transverse_both_diagonals']}, "
              f"quads at depth {[round(x, 1) for x in d['quad_depth_cells']]} "
              f"(patch median {d['median_depth_cells']:.1f})")
        print(f"    clean after erosion of {first_clean} cell(s)"
              if first_clean else "    NOT cleaned by erosion")

    depths = [x for d in results for x in d["quad_depth_cells"]]
    cleaned = [d for d in results
               if any(e["clean"] for e in d["erosion"])]
    print(f"\nparticipating quads at depth <= 1 cell: "
          f"{sum(1 for x in depths if x <= 1.0)}/{len(depths)}")
    print(f"patches cleaned by boundary erosion alone: {len(cleaned)}/{len(results)}")
    Path(a.out).write_text(json.dumps(
        {"n_flagged": len(results), "results": results}, indent=1) + "\n")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
