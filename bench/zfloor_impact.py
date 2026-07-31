"""How much of the audit rests on cells the pipeline would have discarded?

`QuadSurface.cpp` invalidates every tifxyz point with `z <= 0` before it
applies the mask. This reader does not, so our valid set is a SUPERSET of
the one the pipeline sees. That direction is safe for a cleanliness claim
-- we cleaned more than was required -- but it is not safe for a DEFECT
COUNT, because a crossing whose quads live entirely in the discarded
region is a crossing the pipeline never had.

This measures the exposure directly: for every audited segment, how many
reported transverse rows involve a quad with at least one `z <= 0`
corner, and how many segments would change verdict if those rows were
dropped.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import tifffile

DEFAULT_CERTS = Path("results/certificates")
CROSSING_DIRS = [Path(p) for p in (
    "out/crossing", "out/crossing_0139", "out/crossing_0814",
    "out/crossing_1667", "out/crossing_s1")]


def quad_has_zfloor(z: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """(V-1, U-1) bool: quad has a corner we call valid but z <= 0."""
    hit = valid & (z <= 0)
    return (hit[:-1, :-1] | hit[:-1, 1:] | hit[1:, :-1] | hit[1:, 1:])


def load_mesh(path: Path):
    x = tifffile.imread(path / "x.tif").astype(np.float32)
    y = tifffile.imread(path / "y.tif").astype(np.float32)
    z = tifffile.imread(path / "z.tif").astype(np.float32)
    valid = ~((x == -1) & (y == -1) & (z == -1))
    valid &= np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    for name in ("mask.tif", "mask.png"):
        m = path / name
        if m.exists() and m.suffix == ".tif":
            mask = np.asarray(tifffile.imread(m))
            if mask.shape == valid.shape:
                valid &= mask.astype(bool)
    return z, valid


def csvs_for(stem: str) -> list[Path]:
    out = []
    for d in CROSSING_DIRS:
        out.extend(sorted(d.glob(f"{stem}_d*.csv")))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out/zfloor_impact.json")
    ap.add_argument("--certs", default=str(DEFAULT_CERTS),
                    help="directory of per-segment audit certificates")
    a = ap.parse_args()
    certs = Path(a.certs)
    if not certs.is_dir():
        raise SystemExit(f"no certificate directory at {certs}")

    rows = []
    for cert_path in sorted(certs.glob("*_certificate.json")):
        cert = json.loads(cert_path.read_text())
        mesh = Path(cert["mesh"]["path"])
        stem = cert_path.name[: -len("_certificate.json")]
        csvs = csvs_for(stem)
        if not mesh.is_dir() or not csvs:
            continue

        z, valid = load_mesh(mesh)
        bad_quad = quad_has_zfloor(z, valid)
        n_zfloor_cells = int((valid & (z <= 0)).sum())

        total = tainted = 0
        for c in csvs:
            with c.open() as fh:
                for r in csv.DictReader(fh):
                    if r["verdict"] != "transverse":
                        continue
                    total += 1
                    v1, u1, v2, u2 = (int(r["v1"]), int(r["u1"]),
                                      int(r["v2"]), int(r["u2"]))
                    for v, u in ((v1, u1), (v2, u2)):
                        if 0 <= v < bad_quad.shape[0] and 0 <= u < bad_quad.shape[1] \
                                and bad_quad[v, u]:
                            tainted += 1
                            break

        rows.append({
            "segment": stem,
            "zfloor_cells": n_zfloor_cells,
            "transverse_rows": total,
            "rows_touching_zfloor": tainted,
            "fraction": (tainted / total) if total else 0.0,
            "verdict_would_flip": bool(total > 0 and tainted == total),
        })
        if n_zfloor_cells:
            print(f"{stem[:52]:52s} cells={n_zfloor_cells:8d} "
                  f"rows={total:7d} tainted={tainted:7d} "
                  f"{'FLIPS' if total and tainted == total else ''}")

    with_cross = [r for r in rows if r["transverse_rows"] > 0]
    flips = [r for r in rows if r["verdict_would_flip"]]
    tot_rows = sum(r["transverse_rows"] for r in rows)
    tot_taint = sum(r["rows_touching_zfloor"] for r in rows)

    print(f"\nsegments measured                       {len(rows)}")
    print(f"segments with any transverse row        {len(with_cross)}")
    print(f"segments with z<=0 cells we call valid  "
          f"{sum(1 for r in rows if r['zfloor_cells'])}")
    print(f"transverse rows total                   {tot_rows}")
    print(f"rows touching a z<=0 quad               {tot_taint} "
          f"({100.0 * tot_taint / tot_rows:.4f}%)" if tot_rows else "")
    print(f"segments that would lose their crossing verdict  {len(flips)}")
    for r in flips:
        print(f"   {r['segment']}")

    Path(a.out).write_text(json.dumps(
        {"segments": rows,
         "totals": {"segments": len(rows), "transverse_rows": tot_rows,
                    "rows_touching_zfloor": tot_taint,
                    "verdicts_that_would_flip": len(flips)}}, indent=1) + "\n")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
