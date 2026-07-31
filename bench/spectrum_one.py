"""Rerun the intrinsic spectrum for ONE segment, optionally with regions.

The golden spectra (out/spectrum_final_d{0,1}.json) store each event's
numbers but not its region quad sets, so a harness that wants to match
events by REGION SIGNATURE rather than by order needs them re-emitted.
This script reruns exactly what segment_spectrum does for one segment --
same graph, same grouping, same sort, same fields -- and, with
--emit-regions, attaches the parity-normalised quad sets to every event.
The event loop is a line-for-line mirror of segment_spectrum (which does
not expose regions); bench/geodesic_regress.py verifies the mirrored
output against the golden file on every segment, so any drift between
the two is caught, not assumed away.

    uv run python bench/spectrum_one.py "Scroll 5" <segment> \
        --diagonal 0 --emit-regions --out out/one.json
"""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np

from windcheck import tifxyz
from windcheck.check import load_pairs
from windcheck.intrinsic import (MAXEDGE_DEFAULT, SurfaceGraph,
                                 event_separation, oriented_events)

# bench/spectrum_census.py's table, verbatim
CORPORA = [
    ("Scroll 1", "data/scroll1_tifxyz", "20230205180739", "out/crossing_s1"),
    ("Scroll 5", "data/scroll5_tifxyz", "20241024131839", "out/crossing"),
    ("PHerc0139", "data/PHerc0139_tifxyz", "20250728140407", "out/crossing_0139"),
    ("PHerc0814", "data/PHerc0814_tifxyz", "20250804134230", "out/crossing_0814"),
    ("PHerc1667", "data/PHerc1667_tifxyz", "20231117161658", "out/crossing_1667"),
]
RES = re.compile(r"-(\d+\.?\d*)um\.tifxyz$")


def find_segment(corpus: str, segment: str, diagonal: int):
    """Mesh path, pairs CSV path and voxel size, discovered exactly as
    bench/spectrum_census.py discovers them."""
    for name, root, volume, work in CORPORA:
        if name != corpus:
            continue
        d = Path(root) / segment
        if not d.is_dir():
            raise FileNotFoundError(f"no segment dir {d}")
        m = sorted(d.glob(f"mesh/*{volume}*.tifxyz"))
        if not m or not (m[0] / "x.tif").exists():
            raise FileNotFoundError(f"no volume-matched mesh under {d}")
        res = RES.search(m[0].name)
        if not res:
            raise ValueError(f"no resolution in mesh name {m[0].name}")
        csv = Path(work) / f"{segment[:40]}_d{diagonal}.csv"
        return m[0], csv, float(res.group(1))
    raise KeyError(f"unknown corpus {corpus!r}")


def spectrum_with_regions(P, V, rec, voxel_um, diagonal,
                          emit_regions=False, maxedge=MAXEDGE_DEFAULT):
    """segment_spectrum's body, uncapped, with optional region emission."""
    if len(rec) == 0:
        return []
    g = SurfaceGraph(P, V, diagonal, maxedge)
    vv, uu = np.nonzero(g.idx >= 0)
    X = np.empty((g.n, 3))
    X[g.idx[vv, uu]] = P[vv, uu]

    evs = sorted(oriented_events(rec), key=lambda e: -len(e["rows"]))
    mm = voxel_um / 1000.0
    out = []
    for ev in evs:
        arr = rec[ev["rows"]]
        du = np.abs(arr["u1"].astype(int) - arr["u2"].astype(int))
        row = {
            "n_pairs": int(len(arr)),
            "ambiguous": ev["ambiguous"],
            "self_touching": ev["self_touching"],
            "du_max": int(du.max()),
            "median_intersection_length_vx": float(np.nanmedian(arr["pen"])),
        }
        if ev["ambiguous"]:
            row.update({"separation_mm": None, "same_component": None,
                        "endpoint_exact": None, "distance_exact": None})
        else:
            r = event_separation(g, rec, ev, X)
            row.update({
                "separation_mm": (round(r["separation_vx"] * mm, 4)
                                  if r["separation_vx"] is not None else None),
                "same_component": r["same_component"],
                "endpoint_exact": r["endpoint_exact"],
                "distance_exact": r["distance_exact"],
            })
        if emit_regions:
            row["region_a"] = sorted(map(list, ev["region_a"]))
            row["region_b"] = sorted(map(list, ev["region_b"]))
        out.append(row)
    return out


def run_one(corpus: str, segment: str, diagonal: int,
            emit_regions: bool = False) -> dict:
    mesh, csv, vx = find_segment(corpus, segment, diagonal)
    s = tifxyz.read(mesh)
    rec = load_pairs(csv)
    nv, nu = s.points.shape[:2]
    if len(rec):
        keep = ((rec["v1"] < nv - 1) & (rec["v2"] < nv - 1)
                & (rec["u1"] < nu - 1) & (rec["u2"] < nu - 1))
        rec = rec[keep]
    t = time.time()
    events = spectrum_with_regions(s.points, s.valid, rec, vx, diagonal,
                                   emit_regions=emit_regions)
    return {"corpus": corpus, "segment": segment, "diagonal": diagonal,
            "voxel_um": vx, "pairs": int(len(rec)),
            "elapsed_s": round(time.time() - t, 3), "events": events}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("corpus")
    ap.add_argument("segment")
    ap.add_argument("--diagonal", type=int, default=0, choices=(0, 1))
    ap.add_argument("--emit-regions", action="store_true")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    row = run_one(a.corpus, a.segment, a.diagonal, a.emit_regions)
    text = json.dumps(row)
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(text)
        print(f"{row['pairs']} pairs, {len(row['events'])} events "
              f"in {row['elapsed_s']}s -> {a.out}")
    else:
        print(text)


if __name__ == "__main__":
    main()
