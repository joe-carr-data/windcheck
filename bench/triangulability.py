"""Can a twisted quad pair be triangulated at all without self-intersecting?

The local crossings survive both *global* diagonal choices, but that is not the
same as "no planar triangulation exists". Each quad has two possible splits, so a
pair of quads has four combinations, and a mixed assignment might avoid a
crossing that both uniform assignments hit.

This tests all four per pair. If every combination crosses, then no choice of
diagonals realises those two quads as disjoint triangles: the quad grid at that
location is not realisable as an embedded triangulated surface, which is a
statement about the FORMAT rather than about this particular mesh.

An independent triangle-triangle implementation is used deliberately -- it shares
no code with engines/selfcross.cpp, so agreement is a cross-check rather than a
tautology.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crossing_analyse import load                        # noqa: E402

from windcheck import tifxyz                             # noqa: E402


def tri_tri_cross(A: np.ndarray, B: np.ndarray, tol: float = 1e-3) -> bool:
    """True if the two triangle interiors penetrate. Independent of the C++."""
    def plane(T):
        n = np.cross(T[1] - T[0], T[2] - T[0])
        ln = np.linalg.norm(n)
        return (n / ln, ln) if ln > 1e-12 else (None, 0.0)

    n1, l1 = plane(A)
    n2, l2 = plane(B)
    if n1 is None or n2 is None:
        return False
    d1 = np.array([np.dot(n2, p - B[0]) for p in A])
    d2 = np.array([np.dot(n1, p - A[0]) for p in B])
    e = 1e-9
    if (d1 > e).all() or (d1 < -e).all():
        return False
    if (d2 > e).all() or (d2 < -e).all():
        return False
    if np.abs(d1).max() < e:            # coplanar: not a transverse crossing
        return False

    D = np.cross(n1, n2)
    ax = int(np.argmax(np.abs(D)))

    def interval(T, d):
        apex = -1
        for i in range(3):
            j, k = (i + 1) % 3, (i + 2) % 3
            if (d[i] > 0 and d[j] <= 0 and d[k] <= 0) or \
               (d[i] < 0 and d[j] >= 0 and d[k] >= 0):
                apex = i
                break
        if apex < 0:
            return None
        j, k = (apex + 1) % 3, (apex + 2) % 3
        pa = T[apex][ax]
        t0 = pa + (T[j][ax] - pa) * (d[apex] / (d[apex] - d[j]))
        t1 = pa + (T[k][ax] - pa) * (d[apex] / (d[apex] - d[k]))
        return (t0, t1) if t0 <= t1 else (t1, t0)

    ia, ib = interval(A, d1), interval(B, d2)
    if ia is None or ib is None:
        return False
    return (min(ia[1], ib[1]) - max(ia[0], ib[0])) > tol


def quad_tris(P, v, u, diag):
    p00, p10 = P[v, u], P[v + 1, u]
    p01, p11 = P[v, u + 1], P[v + 1, u + 1]
    if diag == 0:
        return [np.array([p00, p01, p11]), np.array([p00, p11, p10])]
    return [np.array([p00, p01, p10]), np.array([p01, p11, p10])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/scroll5_tifxyz"))
    ap.add_argument("--volume", default="20241024131839")
    ap.add_argument("--dir", type=Path, default=Path("out/crossing"))
    ap.add_argument("--cut", type=int, default=200)
    ap.add_argument("--sample", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    tally = {"local": [0, 0, 0], "nonlocal": [0, 0, 0]}   # all4, some, none

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
            idx = rng.permutation(len(sub))[:a.sample]
            for i in idx:
                r = sub[i]
                v1, u1, v2, u2 = (int(r["v1"]), int(r["u1"]),
                                  int(r["v2"]), int(r["u2"]))
                if (v1 + 1 >= V.shape[0] or v2 + 1 >= V.shape[0]
                        or u1 + 1 >= V.shape[1] or u2 + 1 >= V.shape[1]):
                    continue
                hits = 0
                for da in (0, 1):
                    for db in (0, 1):
                        if any(tri_tri_cross(A, B)
                               for A in quad_tris(P, v1, u1, da)
                               for B in quad_tris(P, v2, u2, db)):
                            hits += 1
                if hits == 4:
                    tally[tag][0] += 1
                elif hits > 0:
                    tally[tag][1] += 1
                else:
                    tally[tag][2] += 1

    print("Of the four possible per-quad diagonal assignments, how many cross?\n")
    print(f"{'class':10s} {'all 4':>8s} {'some':>8s} {'none':>8s} "
          f"{'unavoidable':>12s}")
    print("-" * 52)
    for tag in ("local", "nonlocal"):
        a4, so, no = tally[tag]
        tot = a4 + so + no
        if tot == 0:
            continue
        print(f"{tag:10s} {a4:8,d} {so:8,d} {no:8,d} {a4/tot:11.1%}")
    print("\n'all 4' means no choice of diagonals separates the two quads:")
    print("the grid is not realisable as an embedded triangulated surface there.")
    print("'none' means the C++ and this independent implementation disagree,")
    print("which would be a cross-check failure and must be zero-ish.")


if __name__ == "__main__":
    main()
