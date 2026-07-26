"""Check the predicate against FCL, an established collision library.

The earlier Python cross-check shared Moller's interval-overlap structure and was
sampled only from our own positives, so it could confirm plumbing but not catch a
systematic error in the algorithm or a false negative. FCL is a mature C++
collision library using BVH traversal and a different triangle test, so
disagreement means one of us is wrong rather than that we transcribed the same
formula twice.

Both directions are tested:

  POSITIVES  pairs we call TRANSVERSE. FCL should see a collision.
  NEGATIVES  nearby quad pairs we did NOT report. FCL should mostly see none,
             and every case where it does is a candidate false negative of ours.

One expected systematic difference: FCL reports contact, we require strictly
positive penetration above a touch tolerance. So pairs that merely touch should
show up as "FCL yes, we no", and that is correct behaviour rather than a bug.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import trimesh
from scipy.spatial import cKDTree
from trimesh.collision import CollisionManager

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crossing_analyse import load                            # noqa: E402

from windcheck import tifxyz                                 # noqa: E402


def quad_mesh(P: np.ndarray, v: int, u: int) -> trimesh.Trimesh:
    """The two triangles of one quad, diagonal 0, matching the engine."""
    verts = np.array([P[v, u], P[v, u + 1], P[v + 1, u + 1], P[v + 1, u]])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    return trimesh.Trimesh(vertices=verts, faces=faces, process=False)


def collides(P, a, b) -> bool:
    m = CollisionManager()
    m.add_object("a", quad_mesh(P, a[0], a[1]))
    return bool(m.in_collision_single(quad_mesh(P, b[0], b[1])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seg", default="20251115002740-auto_grown_20251115002740308_0_flatboi")
    ap.add_argument("--root", type=Path, default=Path("data/scroll5_tifxyz"))
    ap.add_argument("--volume", default="20241024131839")
    ap.add_argument("--dir", type=Path, default=Path("out/crossing"))
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    d = a.root / a.seg
    path = sorted(d.glob(f"mesh/*{a.volume}*.tifxyz"))[0]
    s = tifxyz.read(path)
    P, V = s.points.astype(np.float64), s.valid
    rec = load(a.dir / f"{a.seg[:40]}_d0.csv")
    print(f"{a.seg[-16:]}: {len(rec):,} transverse pairs on record")

    rng = np.random.default_rng(a.seed)
    reported = {((int(r["v1"]), int(r["u1"])), (int(r["v2"]), int(r["u2"])))
                for r in rec}
    reported |= {(b, aa) for aa, b in reported}

    # --- positives -------------------------------------------------------
    idx = rng.permutation(len(rec))[:a.n]
    pos_agree = pos_total = 0
    for i in idx:
        r = rec[i]
        q1 = (int(r["v1"]), int(r["u1"]))
        q2 = (int(r["v2"]), int(r["u2"]))
        if (q1[0] + 1 >= V.shape[0] or q2[0] + 1 >= V.shape[0]
                or q1[1] + 1 >= V.shape[1] or q2[1] + 1 >= V.shape[1]):
            continue
        pos_total += 1
        if collides(P, q1, q2):
            pos_agree += 1

    # --- negatives: nearby quads we did NOT report ------------------------
    vv, uu = np.nonzero(V[:-1, :-1])
    keep = rng.permutation(len(vv))[:40000]
    vv, uu = vv[keep], uu[keep]
    cent = 0.25 * (P[vv, uu] + P[vv, uu + 1] + P[vv + 1, uu] + P[vv + 1, uu + 1])
    tree = cKDTree(cent)
    pairs = tree.query_pairs(r=25.0, output_type="ndarray")
    rng.shuffle(pairs)

    neg_agree = neg_total = 0
    fcl_yes_we_no = []
    for i, j in pairs:
        if neg_total >= a.n:
            break
        q1 = (int(vv[i]), int(uu[i]))
        q2 = (int(vv[j]), int(uu[j]))
        if abs(q1[0] - q2[0]) <= 1 and abs(q1[1] - q2[1]) <= 1:
            continue                       # adjacency, excluded by both
        if (q1, q2) in reported:
            continue                       # a known positive
        neg_total += 1
        if not collides(P, q1, q2):
            neg_agree += 1
        else:
            fcl_yes_we_no.append((q1, q2))

    print(f"\nPOSITIVES (we say transverse):  FCL agrees {pos_agree}/{pos_total}"
          f"  = {pos_agree/max(pos_total,1):.1%}")
    print(f"NEGATIVES (nearby, unreported): FCL agrees {neg_agree}/{neg_total}"
          f"  = {neg_agree/max(neg_total,1):.1%}")
    if fcl_yes_we_no:
        print(f"\n{len(fcl_yes_we_no)} cases where FCL sees contact and we do not.")
        print("Expected for touching-but-not-penetrating pairs, since FCL reports")
        print("contact while we require positive penetration. Sample:")
        for q1, q2 in fcl_yes_we_no[:5]:
            print(f"    {q1} vs {q2}")


if __name__ == "__main__":
    main()
