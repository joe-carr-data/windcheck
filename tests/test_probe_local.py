"""Equivalence pins for the probe-locality optimization.

`event_pairs_local` must agree with the SurfaceGraph-based
`event_crossing_pairs` on crossing counts, cleared states, and
dropped-quad errors; `apply_field` with a prebuilt neutral graph and a
field cache must produce byte-identical coordinates to the uncached
path. These pins license replacing the O(grid) per-probe rebuild with
the O(event) local predicate inside certified_repair / search_repair.
"""
from __future__ import annotations

import numpy as np
import pytest

from windcheck.check import PAIR_DTYPE
from windcheck.intrinsic import SurfaceGraph, oriented_events
from windcheck.repair import (apply_field, event_crossing_pairs,
                              event_pairs_local)


def _X(g, P):
    vv, uu = np.nonzero(g.idx >= 0)
    X = np.empty((g.n, 3))
    X[g.idx[vv, uu]] = P[vv, uu]
    return X


def planted():
    """The certified-repair fixture: quad (0,0) crosses quad (0,5)."""
    yz = [(0, -0.5), (1, -0.5), (2, -0.5), (3, -0.5),
          (3.5, 2.0), (1.5, 2.0), (1.5, -2.0)]
    P = np.array([[[v, y, z] for (y, z) in yz] for v in (0.0, 1.0)])
    V = np.ones((2, 7), dtype=bool)
    rec = np.array([(0, 1, 0, 5, 1.0, 10.0)], dtype=PAIR_DTYPE)
    ev = oriented_events(rec)[0]
    return P, V, ev


def graph_count(P, V, ev, maxedge, diag):
    g = SurfaceGraph(P, V, diagonal=diag, maxedge=maxedge)
    return event_crossing_pairs(g, _X(g, P), ev)


def test_counts_match_on_planted_crossing():
    P, V, ev = planted()
    for diag in (0, 1):
        assert event_pairs_local(P, V, ev, 10.0, diag) \
            == graph_count(P, V, ev, 10.0, diag)
    assert event_pairs_local(P, V, ev, 10.0, 0) > 0


def test_counts_match_on_cleared_state():
    P, V, ev = planted()
    P2 = P.copy()
    for v, u in ev["region_a"]:                  # lift branch A clear
        for c in ((v, u), (v + 1, u), (v, u + 1), (v + 1, u + 1)):
            P2[c] = P2[c] + np.array([0.0, 0.0, 5.0])
    for diag in (0, 1):
        assert event_pairs_local(P2, V, ev, 100.0, diag) \
            == graph_count(P2, V, ev, 100.0, diag) == 0


def test_dropped_quad_raises_in_both():
    P, V, ev = planted()
    (v0, u0) = sorted(ev["region_a"])[0]
    V2 = V.copy()
    V2[v0, u0] = False                            # invalidate a corner
    with pytest.raises(ValueError):
        graph_count(P, V2, ev, 10.0, 0)
    with pytest.raises(ValueError):
        event_pairs_local(P, V2, ev, 10.0, 0)


def test_maxedge_drop_raises_in_both():
    P, V, ev = planted()
    (v0, u0) = sorted(ev["region_b"])[0]
    P2 = P.copy()
    P2[v0, u0, 2] += 50.0                         # stretch past maxedge
    with pytest.raises(ValueError):
        graph_count(P2, V, ev, 10.0, 0)
    with pytest.raises(ValueError):
        event_pairs_local(P2, V, ev, 10.0, 0)


def test_apply_field_cache_is_bit_identical():
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    gf = SurfaceGraph(np.asarray(P, np.float64), V, diagonal=-1,
                      maxedge=10.0)
    d = np.array([0.0, 0.0, 1.0])
    for mode in ("pinned", "symmetric", 0.3):
        cache: dict = {}
        base, rep0 = apply_field(g, P, V, ev, d, 0.7, mode, 6.0)
        warm1, rep1 = apply_field(g, P, V, ev, d, 0.7, mode, 6.0,
                                  gf=gf, field_cache=cache)
        assert cache                              # field was stored
        warm2, rep2 = apply_field(g, P, V, ev, d, 0.7, mode, 6.0,
                                  gf=gf, field_cache=cache)
        assert np.array_equal(base, warm1)
        assert np.array_equal(base, warm2)
        assert rep0 == rep1 == rep2
        # a different t reuses the cached field, still bit-identical
        b3, _ = apply_field(g, P, V, ev, d, 1.3, mode, 6.0)
        w3, _ = apply_field(g, P, V, ev, d, 1.3, mode, 6.0,
                            gf=gf, field_cache=cache)
        assert np.array_equal(b3, w3)
