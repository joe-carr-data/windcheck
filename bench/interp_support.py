"""Test C: are the flags supported by the raw samples, or made by our interpolation?

Gate pre-registered in notes/PREREG-routeA.md (addendum).

The detector triangulates a 20 vx quadmesh and measures point-to-*interpolated
surface*. The volume-cartographer maintainers point out that in tightly packed
regions with some wobble, interpolating that mesh can reach positions closer
than the samples support. If so, the detector manufactures its own positives.

`d_vertex >= d_interp` by construction and 20 vx sampling inflates `d_vertex`
regardless, so the absolute value is uninformative. The question is whether
flagged cells *separate* from unflagged ones on the raw-sample distance.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from windcheck import atlas, selfgap, tifxyz


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """P(a random positive scores lower than a random negative), rank-based.

    Flagged cells should have *smaller* d_vertex, so the discriminative
    direction is inverted relative to the usual convention.
    """
    both = np.concatenate([pos, neg])
    r = np.argsort(np.argsort(both)).astype(float) + 1
    rp = r[: len(pos)].sum()
    a = (rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))
    return float(1.0 - a)


def vertex_distance(pts: np.ndarray, U: np.ndarray, exclude_u: int,
                    block: int = 64) -> np.ndarray:
    """Nearest same-trace vertex at least `exclude_u` columns away in u.

    Blocked so the u-exclusion can be applied by construction rather than by
    filtering k-nearest results, which would silently truncate when the local
    neighbourhood is larger than k.
    """
    out = np.full(len(pts), np.inf, dtype=np.float64)
    order = np.argsort(U)
    Us, Ps = U[order], pts[order]
    lo_u, hi_u = int(U.min()), int(U.max())
    for start in range(lo_u, hi_u + 1, block):
        stop = start + block
        sel = np.nonzero((U >= start) & (U < stop))[0]
        if len(sel) == 0:
            continue
        # admissible partners: far enough in u from EVERY query in the block
        left = np.searchsorted(Us, start - exclude_u, "left")
        right = np.searchsorted(Us, stop - 1 + exclude_u, "right")
        cand = np.concatenate([Ps[:left], Ps[right:]])
        if len(cand) < 8:
            continue
        d, _ = cKDTree(cand).query(pts[sel], k=1)
        out[sel] = d
    return out


def wobble(surf, stride: int) -> np.ndarray:
    """Discrete Laplacian magnitude of vertex positions: local non-planarity."""
    P = surf.points
    V = surf.valid
    lap = np.full(P.shape[:2], np.nan)
    core = V[1:-1, 1:-1] & V[:-2, 1:-1] & V[2:, 1:-1] & V[1:-1, :-2] & V[1:-1, 2:]
    d = (P[:-2, 1:-1] + P[2:, 1:-1] + P[1:-1, :-2] + P[1:-1, 2:]
         - 4.0 * P[1:-1, 1:-1])
    mag = np.linalg.norm(d, axis=-1)
    mag[~core] = np.nan
    lap[1:-1, 1:-1] = mag
    return lap[::stride, ::stride]


def run(seg_dir: Path, name: str, stride: int, threads: int) -> dict | None:
    m = sorted(seg_dir.glob("mesh/*-on-20241024131839-*.tifxyz"))
    if not m:
        return None
    path = m[0]
    work = Path("out/interp")
    work.mkdir(parents=True, exist_ok=True)

    ex = selfgap.estimate_exclude_u(path, stride, work, threads)
    if ex is None:
        return None
    surf = tifxyz.read(path)
    v, u = np.nonzero(surf.valid[::stride, ::stride])
    V, U = v * stride, u * stride
    if len(V) < 5000:
        return None
    pts = surf.points[V, U]

    atlas.write_atlas([selfgap._Entry(path)], work / "a.bin")
    atlas.write_queries_grouped(pts, U, work / "q.bin")
    res = atlas.run_engine(work / "a.bin", work / "q.bin", work / "r.bin",
                           threads=threads, exclude_u=ex)
    d_interp = res["d1"].astype(np.float64)
    flagged = np.isfinite(d_interp) & (d_interp < selfgap.FLAG_TH)
    if flagged.sum() < 200:
        return None

    d_vertex = vertex_distance(pts, U, ex)
    wob = wobble(surf, stride)[v, u]

    fin = np.isfinite(d_vertex)
    pos = d_vertex[flagged & fin]
    neg = d_vertex[(~flagged) & fin]
    a = auc(pos, neg)

    wf = wob[flagged & np.isfinite(wob)]
    wn = wob[(~flagged) & np.isfinite(wob)]

    out = {
        "trace": name,
        "exclude_u": int(ex),
        "n_queries": int(len(pts)),
        "n_flagged": int(flagged.sum()),
        "d_interp_median_flagged": float(np.median(d_interp[flagged])),
        "d_vertex_median_flagged": float(np.median(pos)),
        "d_vertex_median_unflagged": float(np.median(neg)),
        "d_vertex_p05_unflagged": float(np.percentile(neg, 5)),
        "auc": a,
        "wobble_median_flagged": float(np.median(wf)) if len(wf) else float("nan"),
        "wobble_median_unflagged": float(np.median(wn)) if len(wn) else float("nan"),
        "verdict": ("supported" if a >= 0.90 else
                    "partial" if a >= 0.70 else "NOT SUPPORTED"),
    }
    print(f"{name[:40]:42s} n={out['n_flagged']:6d}  "
          f"d_int {out['d_interp_median_flagged']:5.2f}  "
          f"d_vtx flag {out['d_vertex_median_flagged']:6.2f} vs "
          f"unflag {out['d_vertex_median_unflagged']:6.2f}  "
          f"AUC {a:.3f}  wobble {out['wobble_median_flagged']:.2f}/"
          f"{out['wobble_median_unflagged']:.2f}  {out['verdict']}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/scroll5_tifxyz"))
    ap.add_argument("--stride", type=int, default=3)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--only", default="auto_grown")
    ap.add_argument("--json", type=Path, default=Path("out/interp/testC.json"))
    a = ap.parse_args()

    rows = []
    for d in sorted(a.root.iterdir()):
        if not d.is_dir() or a.only not in d.name:
            continue
        r = run(d, d.name, a.stride, a.threads)
        if r:
            rows.append(r)

    if rows:
        aucs = np.array([r["auc"] for r in rows])
        print(f"\nAUC across {len(rows)} traces: median {np.median(aucs):.3f}  "
              f"range {aucs.min():.3f}-{aucs.max():.3f}")
        print("pre-registered: >=0.90 supported, 0.70-0.90 partial, "
              "<0.70 NOT SUPPORTED")
        a.json.parent.mkdir(parents=True, exist_ok=True)
        a.json.write_text(json.dumps(rows, indent=2))
