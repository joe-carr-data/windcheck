"""Is there a twist magnitude above which a quad cannot be triangulated cleanly?

The triangulability result says ~90% of detected crossings survive every diagonal
choice. If that is governed by quad twist, then there should be a threshold in
twist-relative-to-quad-size above which no planar split works — and a threshold
is something the mesher can act on, where a description is not.

Twist is scale-dependent, so it is normalised by mean edge length: a twist of
2 voxels means something different on a 20-voxel quad than on a 200-voxel one.

Reports the unavoidable fraction binned by normalised twist, and the local vs
nonlocal classes separately, since they had very different twist profiles.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crossing_analyse import load                          # noqa: E402
from triangulability import quad_tris, tri_tri_cross       # noqa: E402

from windcheck import tifxyz                               # noqa: E402


def norm_twist(P: np.ndarray, v: int, u: int) -> float:
    p00, p10 = P[v, u], P[v + 1, u]
    p01, p11 = P[v, u + 1], P[v + 1, u + 1]
    tw = np.linalg.norm(p00 - p10 - p01 + p11)
    e = np.mean([np.linalg.norm(p01 - p00), np.linalg.norm(p11 - p01),
                 np.linalg.norm(p10 - p11), np.linalg.norm(p00 - p10)])
    return float(tw / max(e, 1e-9))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/scroll5_tifxyz"))
    ap.add_argument("--volume", default="20241024131839")
    ap.add_argument("--dir", type=Path, default=Path("out/crossing"))
    ap.add_argument("--cut", type=int, default=200)
    ap.add_argument("--sample", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--json", type=Path,
                    default=Path("out/crossing/twist_threshold.json"))
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    recs = []
    for d in sorted(a.root.iterdir()):
        if not d.is_dir():
            continue
        csv = a.dir / f"{d.name[:40]}_d0.csv"
        m = sorted(d.glob(f"mesh/*{a.volume}*.tifxyz"))
        if not csv.exists() or not m:
            continue
        s = tifxyz.read(m[0])
        P, V = s.points.astype(np.float64), s.valid
        rec = load(csv)
        if len(rec) == 0:
            continue
        sep = np.maximum(np.abs(rec["v1"] - rec["v2"]),
                         np.abs(rec["u1"] - rec["u2"]))
        for tag, sel in (("local", sep <= a.cut), ("nonlocal", sep > a.cut)):
            sub = rec[sel]
            if len(sub) == 0:
                continue
            for i in rng.permutation(len(sub))[:a.sample]:
                r = sub[i]
                v1, u1 = int(r["v1"]), int(r["u1"])
                v2, u2 = int(r["v2"]), int(r["u2"])
                if (v1 + 1 >= V.shape[0] or v2 + 1 >= V.shape[0]
                        or u1 + 1 >= V.shape[1] or u2 + 1 >= V.shape[1]):
                    continue
                hits = sum(
                    1 for da in (0, 1) for db in (0, 1)
                    if any(tri_tri_cross(A, B)
                           for A in quad_tris(P, v1, u1, da)
                           for B in quad_tris(P, v2, u2, db)))
                recs.append({"cls": tag,
                             "twist": max(norm_twist(P, v1, u1),
                                          norm_twist(P, v2, u2)),
                             "hits": hits,
                             "pen": float(r["pen"])})

    arr = np.array([(r["twist"], r["hits"], r["cls"] == "local") for r in recs],
                   dtype=[("tw", "f8"), ("h", "i4"), ("loc", "?")])
    print(f"{len(arr):,} sampled crossing pairs\n")
    edges = [0, 0.02, 0.05, 0.10, 0.20, 0.40, 1.0, 1e9]
    print(f"{'normalised twist':>20s} {'n':>7s} {'unavoidable':>12s} "
          f"{'local n':>8s} {'nonloc n':>9s}")
    print("-" * 62)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (arr["tw"] >= lo) & (arr["tw"] < hi)
        if m.sum() == 0:
            continue
        una = (arr["h"][m] == 4).mean()
        lab = f"{lo:.2f}-{hi:.2f}" if hi < 1e8 else f">{lo:.2f}"
        print(f"{lab:>20s} {int(m.sum()):7,d} {una:11.1%} "
              f"{int((m & arr['loc']).sum()):8,d} "
              f"{int((m & ~arr['loc']).sum()):9,d}")

    for tag, m in (("local", arr["loc"]), ("nonlocal", ~arr["loc"])):
        if m.sum():
            print(f"\n{tag}: median normalised twist "
                  f"{np.median(arr['tw'][m]):.4f}, "
                  f"unavoidable {(arr['h'][m] == 4).mean():.1%}")
    a.json.write_text(json.dumps(recs, indent=2))
    print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
