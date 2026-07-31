"""Equivalence pins for the broad-phase cache in local_field_contacts.

`local_field_contacts(..., broadphase=BroadphaseCache(P_in, V, maxedge))`
must return exactly the contacts of the uncached full-rebuild path: the
cache stores the BASE mesh's retained-quad mask and per-quad corner AABBs,
and per call only quads incident to a moved vertex are recomputed (under
the deformed coordinates B) and patched over the base arrays. Retention
and AABBs are per-quad functions of the four corner coordinates, so this
patched state equals the full recompute from B bitwise.

Contact lists are compared as SORTED lists: the tuple entries are fully
ordered, and sorting makes the pin independent of set-iteration order.
(In practice both paths build moved_q identically and emit the same
order; sorting is belt-and-braces.)
"""
from __future__ import annotations

import numpy as np

from windcheck.check import PAIR_DTYPE
from windcheck.intrinsic import SurfaceGraph, oriented_events, retained_quads
from windcheck.repair import (BroadphaseCache, apply_field,
                              local_field_contacts)


def planted():
    """The certified-repair fixture: quad (0,0) crosses quad (0,5)."""
    yz = [(0, -0.5), (1, -0.5), (2, -0.5), (3, -0.5),
          (3.5, 2.0), (1.5, 2.0), (1.5, -2.0)]
    P = np.array([[[v, y, z] for (y, z) in yz] for v in (0.0, 1.0)])
    V = np.ones((2, 7), dtype=bool)
    rec = np.array([(0, 1, 0, 5, 1.0, 10.0)], dtype=PAIR_DTYPE)
    ev = oriented_events(rec)[0]
    return P, V, ev


def grid20():
    """20x20 grid with validity-mask holes and mild z relief."""
    rng = np.random.default_rng(3)
    nv = nu = 20
    P = np.zeros((nv, nu, 3))
    vv, uu = np.meshgrid(np.arange(nv, dtype=float),
                         np.arange(nu, dtype=float), indexing="ij")
    P[..., 0], P[..., 1] = vv, uu
    P[..., 2] = 0.25 * rng.standard_normal((nv, nu))
    V = np.ones((nv, nu), dtype=bool)
    V[5:8, 5:9] = False                       # interior hole
    V[14, 2] = False                          # lone invalid vertex
    ev = {"region_a": {(2, 2), (2, 3)}, "region_b": {(11, 12)}}
    return P, V, ev


def both(A, B, V, ev, maxedge):
    cold = local_field_contacts(A, B, V, ev, maxedge)
    warm = local_field_contacts(A, B, V, ev, maxedge,
                                broadphase=BroadphaseCache(A, V, maxedge))
    return sorted(cold), sorted(warm)


def quantize(P):
    return P.astype(np.float32).astype(np.float64)


def test_harmonic_deformation_matches_on_planted_grid():
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    P2, _ = apply_field(g, P, V, ev, np.array([0.0, 0.0, 1.0]), 0.7,
                        "symmetric", 6.0)
    B = quantize(P2)
    cold, warm = both(P, B, V, ev, 10.0)
    assert cold == warm


def test_harmonic_deformation_matches_on_20x20_grid_with_holes():
    P, V, ev = grid20()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=4.0)
    bpc = BroadphaseCache(P, V, 4.0)
    for d, t in ((np.array([0.0, 0.0, 1.0]), 1.5),
                 (np.array([0.3, -0.2, 0.9]), 0.8)):
        P2, _ = apply_field(g, P, V, ev, d, t, "symmetric", 8.0)
        B = quantize(P2)
        cold = sorted(local_field_contacts(P, B, V, ev, 4.0))
        # one cache object reused across different deformations of one base
        warm = sorted(local_field_contacts(P, B, V, ev, 4.0,
                                           broadphase=bpc))
        assert cold == warm


def test_retention_change_newly_dropped_matches_full_recompute():
    P, V, ev = grid20()
    maxedge = 2.5
    B = quantize(P.copy())
    B[10, 3, 2] += 3.0                        # edge > maxedge: quads dropped
    Q0 = retained_quads(np.asarray(P, np.float64), V, maxedge)
    Q1 = retained_quads(B, V, maxedge)
    assert (Q0 != Q1).any()                   # retention really changed
    cold, warm = both(P, B, V, ev, maxedge)
    assert cold == warm


def test_retention_change_newly_retained_matches_full_recompute():
    P, V, ev = grid20()
    maxedge = 2.5
    A = P.copy()
    A[10, 3, 2] += 3.0                        # base has a dropped quad
    B = quantize(A.copy())
    B[10, 3, 2] -= 3.0                        # deformation re-retains it
    Q0 = retained_quads(np.asarray(A, np.float64), V, maxedge)
    Q1 = retained_quads(np.asarray(B, np.float64), V, maxedge)
    assert (~Q0 & Q1).any()                   # a quad became retained
    cold, warm = both(A, B, V, ev, maxedge)
    assert cold == warm


def test_actual_contact_reported_identically():
    P, V, ev = planted()
    # deform quad (0,0) through the vertical sheet quad (0,5) (plane
    # y=1.5, z in [-2,2]): its corners move to z=0, y straddling 1.5 --
    # a genuine crossing at Chebyshev grid distance 5
    B = quantize(P.copy())
    B[0, 0], B[1, 0] = [0.0, 1.2, 0.0], [1.0, 1.2, 0.0]
    B[0, 1], B[1, 1] = [0.0, 1.8, 0.0], [1.0, 1.8, 0.0]
    cold, warm = both(P, B, V, ev, 10.0)
    assert cold == warm
    assert len(cold) > 0                      # the contact is really there


def bigrid(nv=8, nu=8, spacing=30.0, seed=7):
    """Planted grid whose quads straddle multiple 40.0-voxel census cells
    (spacing 30 > cell/2, coords up to 210: cells 0..5 per axis)."""
    rng = np.random.default_rng(seed)
    P = np.zeros((nv, nu, 3))
    vv, uu = np.meshgrid(np.arange(nv, dtype=float),
                         np.arange(nu, dtype=float), indexing="ij")
    P[..., 0], P[..., 1] = spacing * vv, spacing * uu
    P[..., 2] = 2.0 * rng.standard_normal((nv, nu))
    V = np.ones((nv, nu), dtype=bool)
    ev = {"region_a": {(1, 1)}, "region_b": {(6, 6)}}
    return P, V, ev


def test_bucket_candidates_equal_fullscan_on_base_aabbs():
    """query() ∪ exact-AABB-filter == full scan over the base AABBs, for
    boxes that sit inside one cell and boxes straddling cell boundaries."""
    P, V, ev = bigrid()
    bp = BroadphaseCache(P, V, 100.0)
    n = len(bp.flat0)
    assert n > 0 and len(bp.buckets) > 1      # index really spans cells
    rng = np.random.default_rng(11)
    boxes = [  # (blo, bhi): cell-interior, boundary-straddling, cell-sized
        (np.array([5.0, 5.0, -3.0]), np.array([12.0, 12.0, 3.0])),
        (np.array([38.0, 38.0, -1.0]), np.array([42.0, 42.0, 1.0])),
        (np.array([75.0, 115.0, -5.0]), np.array([125.0, 165.0, 5.0])),
    ]
    for _ in range(20):
        c = rng.uniform(-10.0, 220.0, 3) * np.array([1, 1, 0.02])
        h = rng.uniform(1.0, 45.0, 3)
        boxes.append((c - h, c + h))
    for blo, bhi in boxes:
        full = np.nonzero(np.all(bp.qhi0 >= blo, axis=1)
                          & np.all(bp.qlo0 <= bhi, axis=1))[0]
        cand = bp.query(blo, bhi)
        assert set(full) <= set(cand)         # conservative prefilter
        kept = cand[np.all(bp.qhi0[cand] >= blo, axis=1)
                    & np.all(bp.qlo0[cand] <= bhi, axis=1)]
        assert kept.tolist() == full.tolist() # exact after the AABB test


def test_moved_quad_aabb_crossing_cell_boundaries_matches():
    """A moved quad whose swept AABB straddles 40.0-cell boundaries in all
    axes, planted through a distant obstacle: bucket path == full scan."""
    P, V, ev = bigrid()
    B = P.copy()
    # carry quad (1,1) across cell boundaries onto quad (5,5)'s interior
    # (base corners of (5,5) at 150..180: cells 3-4) and pierce its surface
    for (v, u) in ((1, 1), (1, 2), (2, 1), (2, 2)):
        B[v, u] = [155.0 + 12.0 * (v - 1), 155.0 + 12.0 * (u - 1),
                   -6.0 + 12.0 * (v - 1)]
    B = quantize(B)
    maxedge = 400.0                           # keep transition quads
    Q1 = retained_quads(np.asarray(B, np.float64), V, maxedge)
    assert Q1[1, 1] and Q1[5, 5]
    # the swept AABB of (1,1) spans multiple census cells in every axis
    c0 = np.array([P[1, 1], P[2, 1], P[1, 2], P[2, 2]])
    c1 = np.array([B[1, 1], B[2, 1], B[1, 2], B[2, 2]])
    blo = np.minimum(c0.min(0), c1.min(0))
    bhi = np.maximum(c0.max(0), c1.max(0))
    assert (np.floor(bhi / 40.0) - np.floor(blo / 40.0) >= 1).all()
    cold, warm = both(P, B, V, ev, maxedge)
    assert cold == warm
    assert len(cold) > 0


def test_moved_obstacle_leaving_base_cell_matches():
    """A moved OBSTACLE whose AABB left its base cells entirely: it can
    only be found via the always-candidate moved set, never the buckets."""
    P, V, ev = bigrid()
    B = P.copy()
    # obstacle quad (6,6) (base corners 180..210, cells 4-5) relocates
    # wholesale next to quad (1,1) (corners 30..60) and tilts through the
    # plane z~0 that quad (1,1) spans -- its new AABB shares no cell with
    # its base AABB
    for (v, u) in ((6, 6), (6, 7), (7, 6), (7, 7)):
        B[v, u] = [38.0 + 10.0 * (v - 6), 38.0 + 10.0 * (u - 6),
                   -8.0 + 16.0 * (v - 6)]
    # nudge quad (1,1)'s corners so it is itself a moved quad too
    for (v, u) in ((1, 1), (1, 2), (2, 1), (2, 2)):
        B[v, u, 2] = 0.0
    B = quantize(B)
    maxedge = 400.0
    Q1 = retained_quads(np.asarray(B, np.float64), V, maxedge)
    assert Q1[1, 1] and Q1[6, 6]
    bp = BroadphaseCache(P, V, maxedge)
    # the relocated obstacle's current AABB overlaps no bucket that its
    # base AABB was registered in
    new_lo = np.array([B[6, 6], B[7, 6], B[6, 7], B[7, 7]]).min(0)
    new_hi = np.array([B[6, 6], B[7, 6], B[6, 7], B[7, 7]]).max(0)
    i66 = int(np.nonzero(bp.flat0 == 6 * (V.shape[1] - 1) + 6)[0][0])
    assert i66 not in set(bp.query(new_lo, new_hi).tolist())
    cold, warm = both(P, B, V, ev, maxedge)
    assert cold == warm
    assert len(cold) > 0


def test_deformed_branch_contact_on_20x20():
    P, V, ev = grid20()
    # carry quad (2,2) over quad (10,10) and tilt it through that quad's
    # z~0 surface: a crossing between quads at Chebyshev distance 8
    B = P.copy()
    for (v, u) in ((2, 2), (2, 3), (3, 2), (3, 3)):
        B[v, u] = [10.2 + 0.6 * (v - 2), 10.2 + 0.6 * (u - 2),
                   -1.0 + 2.0 * (v - 2)]
    B = quantize(B)
    maxedge = 20.0                            # keep transition quads retained
    Q1 = retained_quads(np.asarray(B, np.float64), V, maxedge)
    assert Q1[2, 2] and Q1[10, 10]            # both quads survive
    cold, warm = both(P, B, V, ev, maxedge)
    assert cold == warm
    assert len(cold) > 0
