"""Analytic tests for the intrinsic-separation engine.

Every distance asserted here is computed by hand on a synthetic surface first.
The engine replaces two retracted attempts (rounds 11 and 12); each retraction
gets a pinning test: the graph must be the retained-triangle complex (a
one-row slit DISCONNECTS it, a single-vertex chain is not surface, a
maxedge-dropped quad contributes nothing), grouping must preserve branch
identity (straight/swapped parity, conflicts flagged ambiguous), and event
endpoints are the recomputed triangle-intersection points attached
barycentrically -- pinned by a planted transverse crossing whose separation
is known exactly.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from windcheck.check import PAIR_DTYPE
from windcheck.intrinsic import (SurfaceGraph, _tri_tri_segment,
                                 oriented_events, event_separation,
                                 segment_spectrum)

SQ2 = math.sqrt(2.0)


def plane(n: int = 7, spacing: float = 1.0):
    v, u = np.mgrid[0:n, 0:n].astype(np.float64)
    P = np.stack([u * spacing, v * spacing, np.zeros_like(u)], axis=-1)
    V = np.ones((n, n), dtype=bool)
    return P, V


def pairs(*quads) -> np.ndarray:
    return np.array([(v1, u1, v2, u2, 1.0, 10.0) for v1, u1, v2, u2 in quads],
                    dtype=PAIR_DTYPE)


def graph_X(g: SurfaceGraph, P: np.ndarray) -> np.ndarray:
    vv, uu = np.nonzero(g.idx >= 0)
    X = np.empty((g.n, 3))
    X[g.idx[vv, uu]] = P[vv, uu]
    return X


# ------------------------------------------------------------ graph = complex
def test_row_distance_is_exact():
    P, V = plane(7)
    g = SurfaceGraph(P, V, diagonal=0)
    assert g.vertex_distance((0, 0), (0, 6)) == pytest.approx(6.0)
    assert g.vertex_distance((3, 1), (5, 1)) == pytest.approx(2.0)


def test_diagonal_orientation_matters():
    """d0 runs down-right, d1 up-right; a corner-to-corner path must feel it."""
    P, V = plane(7)
    g0 = SurfaceGraph(P, V, diagonal=0)
    g1 = SurfaceGraph(P, V, diagonal=1)
    assert g0.vertex_distance((0, 0), (5, 5)) == pytest.approx(5 * SQ2)
    assert g1.vertex_distance((0, 0), (5, 5)) == pytest.approx(10.0)
    assert g1.vertex_distance((5, 0), (0, 5)) == pytest.approx(5 * SQ2)
    assert g0.vertex_distance((5, 0), (0, 5)) == pytest.approx(10.0)


def test_one_row_slit_disconnects():
    """A one-vertex-wide opening supports no triangle: the surface is CUT.

    The round-12 correction: the old vertex-grid graph walked through this
    opening (14 + 2*sqrt(2)); the censused triangle complex does not contain
    it, so the two sides are different components.
    """
    P, V = plane(7)
    V[0:6, 3] = False
    for diag in (0, 1):
        g = SurfaceGraph(P, V, diagonal=diag)
        a, b = g.vertex(0, 0), g.vertex(0, 6)
        assert a >= 0 and b >= 0
        assert g.comp[a] != g.comp[b]


def test_two_row_passage_forces_the_detour():
    """A two-row opening admits one retained quad column: detour is exact.

    Hand computation, d0: (0,0) -diag- (1,1) -diag- (2,2), down 3 to (5,2),
    across (5,3),(5,4), up 5 to (0,4), right 2 -> 2*sqrt(2)+3+2+5+2 = 12+2√2.
    d1 mirrors it (down 5, across, up-right diagonals): same total.
    """
    P, V = plane(7)
    V[0:5, 3] = False
    for diag in (0, 1):
        g = SurfaceGraph(P, V, diagonal=diag)
        assert g.vertex_distance((0, 0), (0, 6)) == pytest.approx(12 + 2 * SQ2)


def test_single_vertex_chain_is_not_surface():
    """A chain of valid vertices with no quad around it joins nothing."""
    P, V = plane(5)
    V[:, 2] = False
    V[2, 2] = True           # lone bridge vertex between the two blocks
    g = SurfaceGraph(P, V)
    assert g.idx[2, 2] == -1                 # not in any retained quad
    assert g.ncomp == 2


def test_maxedge_drops_the_quad():
    """A quad with a stretched edge is dropped exactly as selfcross drops it."""
    P, V = plane(4)
    P[1, 1, 2] = 5.0                          # spike: edges from (1,1) > 2
    g = SurfaceGraph(P, V, maxedge=2.0)
    assert g.idx[1, 1] == -1                  # vertex belongs to no quad
    assert g.idx[0, 0] == -1                  # its only quad was dropped
    assert g.ncomp == 1                       # the rest stays one sheet


# --------------------------------------------------------------- grouping
def test_oriented_grouping_straight_and_swapped():
    rec = pairs((0, 0, 10, 0),
                (1, 0, 11, 0),      # straight continuation
                (12, 0, 2, 0),      # same event, recorded swapped
                (0, 20, 10, 20))    # unrelated event
    evs = sorted(oriented_events(rec), key=lambda e: -len(e["rows"]))
    assert len(evs) == 2
    big = evs[0]
    assert len(big["rows"]) == 3 and not big["ambiguous"]
    assert big["region_a"] == {(0, 0), (1, 0), (2, 0)}
    assert big["region_b"] == {(10, 0), (11, 0), (12, 0)}
    assert not big["self_touching"]


def test_orientation_conflict_is_ambiguous():
    """A pair adjacent to another on BOTH assignments cannot be oriented."""
    rec = pairs((0, 0, 2, 0), (1, 0, 1, 1))
    evs = oriented_events(rec)
    assert len(evs) == 1
    assert evs[0]["ambiguous"]
    assert evs[0]["self_touching"]


def test_ambiguous_events_are_reported_not_measured():
    P, V = plane(8)
    rec = pairs((0, 0, 2, 0), (1, 0, 1, 1))
    spec = segment_spectrum(P, V, rec, voxel_um=1000.0)
    assert len(spec) == 1
    assert spec[0]["ambiguous"] and spec[0]["separation_mm"] is None


# ------------------------------------------------------------- endpoints
def test_tri_tri_segment_exact():
    """Flat triangle vs vertical wall: segment is y=1, z=0, x in [0,3]."""
    T1 = np.array([[0, 0, 0], [4, 0, 0], [0, 4, 0]], float)
    T2 = np.array([[0, 1, -1], [4, 1, -1], [0, 1, 3]], float)
    seg = _tri_tri_segment(T1, T2)
    assert seg is not None
    got = sorted(map(tuple, np.round(seg, 9)))
    assert got[0] == pytest.approx((0.0, 1.0, 0.0))
    assert got[1] == pytest.approx((3.0, 1.0, 0.0))
    # parallel planes: no segment
    assert _tri_tri_segment(T1, T1 + [0, 0, 5.0]) is None


def test_planted_crossing_measured_barycentrically():
    """A strip that pierces itself; the walk locus-to-locus is known exactly.

    Two rows (x = v), columns trace a curve in (y, z): flat run at z=-0.5
    through y=0..3, up to (3.5, 2), back to (1.5, 2), then straight down to
    (1.5, -2) -- the last column drop is a wall at y=1.5 crossing the flat
    run at (y, z) = (1.5, -0.5). Crossing pair: base quad (0,1) vs wall quad
    (0,5). Shortest walk from the crossing locus via the base sheet to the
    wall sheet: 0.5 (locus to corner y=2) + 1 + sqrt(6.5) + 2 + 2.5 (corner
    at (1.5,2) back down to the locus) = 6 + sqrt(6.5).
    """
    yz = [(0, -0.5), (1, -0.5), (2, -0.5), (3, -0.5),
          (3.5, 2.0), (1.5, 2.0), (1.5, -2.0)]
    P = np.array([[[v, y, z] for (y, z) in yz] for v in (0.0, 1.0)])
    V = np.ones((2, 7), dtype=bool)
    rec = pairs((0, 1, 0, 5))
    spec = segment_spectrum(P, V, rec, voxel_um=1000.0)
    assert len(spec) == 1
    e = spec[0]
    assert e["endpoint_exact"] is True
    assert e["same_component"] is True
    assert e["separation_mm"] == pytest.approx(6 + math.sqrt(6.5), abs=1e-3)


def test_seeded_distance_matches_field_and_fallback():
    P, V = plane(7)
    V[0:5, 3] = False                      # forces the 12+2sqrt(2) detour
    g = SurfaceGraph(P, V, diagonal=0)
    a = [(g.vertex(0, 0), 0.0)]
    b = [(g.vertex(0, 6), 0.0)]
    exact = 12 + 2 * SQ2
    d, ok = g.seeded_distance(a, b)
    assert ok and d == pytest.approx(exact)
    # budget exhaustion must stay exact even with UNEQUAL offsets -- the
    # round-13 defect: the old fallback mixed one source's path with
    # another source's smaller offset. Seed a decoy nearer the target in
    # graph terms but carrying a large offset; correct answer must pick the
    # (0,0) seed's full detour, not combine decoy path with 0.0 offset.
    a2 = [(g.vertex(0, 0), 0.0), (g.vertex(0, 5), 100.0)]
    want = min(exact, 100.0 + 1.0)
    d2, ok2 = g.seeded_distance(a2, b, max_pops=1)
    assert ok2 and d2 == pytest.approx(want)
    d3, ok3 = g.seeded_distance(a2, b)
    assert ok3 and d3 == pytest.approx(want)


def test_disconnection_is_an_event_property():
    P, V = plane(10)
    V[:, 4:6] = False                      # two islands
    g = SurfaceGraph(P, V)
    X = graph_X(g, P)
    assert g.ncomp == 2
    rec = pairs((0, 0, 0, 7))
    ev = oriented_events(rec)[0]
    r = event_separation(g, rec, ev, X)
    assert r["separation_vx"] is None and r["same_component"] is False
    rec2 = pairs((0, 0, 0, 2))
    ev2 = oriented_events(rec2)[0]
    r2 = event_separation(g, rec2, ev2, X)
    # flat plane: no true intersection, corner fallback; quads (0,0),(0,2)
    # nearest corners u=1 and u=2 -> exactly 1
    assert r2["same_component"] is True
    assert r2["separation_vx"] == pytest.approx(1.0)


def test_spectrum_cap_keeps_largest_events():
    P, V = plane(14)
    rec = pairs((0, 0, 0, 8), (5, 0, 5, 8), (6, 1, 6, 9))
    spec = segment_spectrum(P, V, rec, voxel_um=1000.0, max_events=1)
    assert len(spec) == 1 and spec[0]["n_pairs"] == 2
