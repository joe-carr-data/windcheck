"""Are the crossings in bad mesh, or in good mesh?

This is the last thing that can invalidate the census. If crossings cluster in
quads that are twisted, stretched or nearly degenerate, then what we found is an
inadequate mesh representation rather than a defect in the traced surface. Both
are real, but they are different findings with different consumers.

Per participating quad, with corners p00 p10 p01 p11:

  twist       |p00 - p10 - p01 + p11|, the bilinear cross-term. This is the
              quantity a discrete Laplacian is blind to -- a saddle z=uv has zero
              Laplacian and nonzero twist -- and it is exactly what the
              volume-cartographer maintainer described when he said the quadmesh
              can "bilinearly interpolate to even closer positions". Our earlier
              Laplacian check reported no effect; that check could not have found
              one.
  planarity   RMS distance of the four corners from their best-fit plane
  aspect      longest edge / shortest edge
  max_edge    longest edge, in voxels; the grid pitch is ~20

Each metric is compared against the same metric over the whole surface, so the
question is always "is this region unusual for this mesh", never "is this number
big".
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crossing_analyse import components, load          # noqa: E402

from windcheck import tifxyz                            # noqa: E402


def quad_metrics(P: np.ndarray, V: np.ndarray):
    """Per-quad twist, planarity, aspect and max edge over the whole grid."""
    p00 = P[:-1, :-1]
    p10 = P[1:, :-1]
    p01 = P[:-1, 1:]
    p11 = P[1:, 1:]
    ok = V[:-1, :-1] & V[1:, :-1] & V[:-1, 1:] & V[1:, 1:]

    twist = np.linalg.norm(p00 - p10 - p01 + p11, axis=-1)

    e = np.stack([np.linalg.norm(p01 - p00, axis=-1),
                  np.linalg.norm(p11 - p01, axis=-1),
                  np.linalg.norm(p10 - p11, axis=-1),
                  np.linalg.norm(p00 - p10, axis=-1)], axis=-1)
    max_edge = e.max(-1)
    aspect = max_edge / np.maximum(e.min(-1), 1e-6)

    # planarity: the four corners' RMS deviation from their own best-fit plane.
    # For a bilinear quad this is governed by the twist term, so it is reported
    # as a cross-check rather than as independent evidence.
    ctr = (p00 + p10 + p01 + p11) / 4.0
    n = np.cross(p11 - p00, p10 - p01)
    ln = np.linalg.norm(n, axis=-1, keepdims=True)
    nn = n / np.maximum(ln, 1e-9)
    planar = np.zeros_like(twist)
    for c in (p00, p10, p01, p11):
        planar += np.abs(((c - ctr) * nn).sum(-1)) ** 2
    planar = np.sqrt(planar / 4.0)

    for a in (twist, max_edge, aspect, planar):
        a[~ok] = np.nan
    return {"twist": twist, "max_edge": max_edge, "aspect": aspect,
            "planarity": planar}


def describe(vals: np.ndarray) -> dict:
    v = vals[np.isfinite(vals)]
    if len(v) == 0:
        return {"n": 0}
    return {"n": int(len(v)), "median": float(np.median(v)),
            "p90": float(np.percentile(v, 90)),
            "p99": float(np.percentile(v, 99)), "max": float(v.max())}


def analyse(seg_dir: Path, csv: Path, volume: str, cut: int) -> dict | None:
    m = sorted(seg_dir.glob(f"mesh/*{volume}*.tifxyz"))
    if not m:
        return None
    s = tifxyz.read(m[0])
    met = quad_metrics(s.points, s.valid)
    rec = load(csv)
    if len(rec) == 0:
        return None

    sep = np.maximum(np.abs(rec["v1"] - rec["v2"]), np.abs(rec["u1"] - rec["u2"]))
    out: dict = {"segment": seg_dir.name, "grid": list(s.shape)}

    for tag, sel in (("local", sep <= cut), ("nonlocal", sep > cut)):
        sub = rec[sel]
        if len(sub) == 0:
            out[tag] = {"pairs": 0}
            continue
        # Deduplicate quad origins. Without this a quad taking part in many
        # crossing pairs is weighted by its pair count, so the statistic is
        # pair-weighted while being described as per-quad.
        pairs = np.stack([np.concatenate([sub["v1"], sub["v2"]]),
                          np.concatenate([sub["u1"], sub["u2"]])], axis=1)
        uniq = np.unique(pairs, axis=0)
        vs, us = uniq[:, 0], uniq[:, 1]
        inb = (vs < met["twist"].shape[0]) & (us < met["twist"].shape[1])
        vs, us = vs[inb], us[inb]
        lab, nreg = components(sub)
        ent = {"pairs": int(len(sub)), "regions": nreg,
               "penetration_median": float(np.nanmedian(sub["pen"])),
               "angle_median": float(np.nanmedian(sub["ang"]))}
        for k, arr in met.items():
            ent[k] = describe(arr[vs, us])
        out[tag] = ent

    out["whole_surface"] = {k: describe(v.ravel()) for k, v in met.items()}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/scroll5_tifxyz"))
    ap.add_argument("--volume", default="20241024131839")
    ap.add_argument("--dir", type=Path, default=Path("out/crossing"))
    ap.add_argument("--cut", type=int, default=200)
    ap.add_argument("--json", type=Path, default=Path("out/crossing/quality.json"))
    a = ap.parse_args()

    rows = []
    for d in sorted(a.root.iterdir()):
        if not d.is_dir():
            continue
        csv = a.dir / f"{d.name[:40]}_d0.csv"
        r = analyse(d, csv, a.volume, a.cut)
        if r:
            rows.append(r)

    def agg(rs, tag, key):
        vals = [r[tag][key]["median"] for r in rs
                if r.get(tag, {}).get("pairs") and r[tag].get(key, {}).get("n")]
        return float(np.median(vals)) if vals else float("nan")

    lab = [r for r in rows if "-w0" in r["segment"]]
    ag = [r for r in rows if "auto_grown" in r["segment"]]

    print(f"Quality of quads participating in crossings, vs the whole surface.")
    print(f"separation cut = {a.cut} cells\n")
    hdr = f"{'population':22s} {'twist':>9s} {'planarity':>10s} {'aspect':>8s} {'max_edge':>9s}"
    print(hdr); print("-" * len(hdr))
    for name, rs in (("labelled (44)", lab), ("multi-wrap (9)", ag)):
        for tag in ("whole_surface", "local", "nonlocal"):
            if tag == "whole_surface":
                vals = [r["whole_surface"][k]["median"] for r in rs for k in ("twist",)]
                row = [float(np.median([r["whole_surface"][k]["median"] for r in rs]))
                       for k in ("twist", "planarity", "aspect", "max_edge")]
            else:
                row = [agg(rs, tag, k)
                       for k in ("twist", "planarity", "aspect", "max_edge")]
            label = f"{name} {tag}"
            print(f"{label:22s} " + " ".join(
                f"{v:9.3f}" if np.isfinite(v) else f"{'--':>9s}" for v in row))
        print()

    a.json.write_text(json.dumps(rows, indent=2))
    print(f"wrote {a.json}")


if __name__ == "__main__":
    main()
