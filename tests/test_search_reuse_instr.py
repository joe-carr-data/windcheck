"""Round-24 pins: apply_field patch-not-copy workspace, search_repair
one-pass extra candidates, and the kernel profiler.

1. `apply_field(..., workspace=ws)` must be BYTE-identical to the
   fresh-copy path across a sequence of different (direction, lam, t)
   probes -- the buffer resets only previously-moved vertices and applies
   the identical in-place arithmetic.
2. `search_repair(..., extra_candidates=k)` must return the plain call's
   candidates first, in the plain call's exact order (selection AND
   ranking), with the extras strictly appended -- the compatibility
   guarantee that lets transact() re-search at most once without changing
   any golden candidate sequence.
3. The kernel profiler is observation-only: counters populate, reset
   clears, and instrumented functions return unchanged values.
"""
from __future__ import annotations

import numpy as np

from windcheck.check import PAIR_DTYPE
from windcheck.intrinsic import SurfaceGraph, oriented_events
from windcheck.repair import (KERNEL_PROFILE, apply_field,
                              event_pairs_local, kernel_profile_reset,
                              kernel_profile_snapshot, search_repair)


def planted():
    """The certified-repair fixture: quad (0,0) crosses quad (0,5)."""
    yz = [(0, -0.5), (1, -0.5), (2, -0.5), (3, -0.5),
          (3.5, 2.0), (1.5, 2.0), (1.5, -2.0)]
    P = np.array([[[v, y, z] for (y, z) in yz] for v in (0.0, 1.0)])
    V = np.ones((2, 7), dtype=bool)
    rec = np.array([(0, 1, 0, 5, 1.0, 10.0)], dtype=PAIR_DTYPE)
    ev = oriented_events(rec)[0]
    return P, V, ev


def test_apply_field_workspace_bit_identical_over_probe_sequence():
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    gf = SurfaceGraph(np.asarray(P, np.float64), V, diagonal=-1,
                      maxedge=10.0)
    fcache: dict = {}
    ws: dict = {}
    rng = np.random.default_rng(3)
    probes = [(np.array([0.0, 0.0, 1.0]), "symmetric", 0.7),
              (np.array([0.0, 0.0, 1.0]), "symmetric", 1.3),
              (np.array([0.0, 1.0, 0.5]), 0.3, 0.9),
              (np.array([0.0, -1.0, 1.0]), "pinned", 0.4),
              (np.array([0.0, 0.2, -1.0]), 0.7, 1.1),
              (np.array([0.0, 0.0, 1.0]), "symmetric", 0.7)]
    probes += [(rng.normal(size=3) + np.array([0.0, 0.0, 2.0]),
                float(rng.uniform(0.1, 0.9)), float(rng.uniform(0.2, 1.8)))
               for _ in range(6)]
    for d, mode, t in probes:
        fresh, rep_f = apply_field(g, P, V, ev, d, t, mode, 6.0,
                                   gf=gf, field_cache=fcache)
        buf, rep_b = apply_field(g, P, V, ev, d, t, mode, 6.0,
                                 gf=gf, field_cache=fcache, workspace=ws)
        assert fresh.tobytes() == buf.tobytes()      # byte-identical
        assert rep_f == rep_b
    assert "base" in ws and ws["touched"]            # buffer really engaged


def test_workspace_resets_vertices_the_next_probe_does_not_move():
    """A symmetric probe (partner moves) followed by a pinned probe
    (partner fixed at 0): vertices moved only by the first must be
    restored to base bitwise by the buffer reset."""
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    gf = SurfaceGraph(np.asarray(P, np.float64), V, diagonal=-1,
                      maxedge=10.0)
    ws: dict = {}
    d = np.array([0.0, 0.0, 1.0])
    apply_field(g, P, V, ev, d, 1.0, "symmetric", 8.0, gf=gf, workspace=ws)
    sym_touched = set(ws["touched"])
    buf, _ = apply_field(g, P, V, ev, d, 1.0, "pinned", 8.0,
                         gf=gf, workspace=ws)
    fresh, _ = apply_field(g, P, V, ev, d, 1.0, "pinned", 8.0, gf=gf)
    assert set(ws["touched"]) < sym_touched      # partner no longer moves
    assert fresh.tobytes() == buf.tobytes()


def test_search_repair_extra_candidates_preserves_head_order():
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    plain = search_repair(g, P, V, ev, budget_point_vx=3.0, support_vx=6.0)
    more = search_repair(g, P, V, ev, budget_point_vx=3.0, support_vx=6.0,
                         extra_candidates=6)
    assert len(plain) == 6                       # fixture saturates the scan
    assert len(more) > len(plain)                # extras really collected
    for a, b in zip(plain, more):                # head: same candidates,
        assert a["direction"] == b["direction"]  # same order
        assert a["lam"] == b["lam"]
        assert a["t_rel"] == b["t_rel"]
        assert a["max_point_vx"] == b["max_point_vx"]
        assert np.array_equal(a["P32"], b["P32"])
        assert a["field_report"] == b["field_report"]
    ranks = [c["max_point_vx"] for c in more[len(plain):]]
    assert ranks == sorted(ranks)                # tail ranked best-first
    assert all(not c.get("extra") for c in more[:len(plain)])
    assert all(c.get("extra") for c in more[len(plain):])
    assert all("extra" not in c for c in plain)  # plain call dicts unchanged


def test_search_repair_extra_zero_is_the_plain_scan():
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    a = search_repair(g, P, V, ev, budget_point_vx=3.0, support_vx=6.0)
    b = search_repair(g, P, V, ev, budget_point_vx=3.0, support_vx=6.0,
                      extra_candidates=0)
    assert [(c["direction"], c["lam"], c["t_rel"]) for c in a] \
        == [(c["direction"], c["lam"], c["t_rel"]) for c in b]


def test_kernel_profile_counts_and_reset():
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    kernel_profile_reset()
    n0 = event_pairs_local(P, V, ev, 10.0, 0)
    search_repair(g, P, V, ev, budget_point_vx=3.0, support_vx=6.0)
    snap = kernel_profile_snapshot()
    assert snap["calls"]["search_repair"] == 1
    assert snap["calls"]["event_pairs_local"] >= 1
    assert snap["calls"]["apply_field"] >= 1
    assert snap["calls"]["harmonic_field"] >= 1
    assert snap["counts"]["field_cache_hits"] >= 1
    assert snap["calls"]["local_field_contacts"] >= 1
    assert snap["counts"]["tri_tri_tests"] > 0
    assert snap["counts"]["broadphase_candidates"] > 0
    assert snap["calls"]["pair_intervals_batched_crossing"] > 0
    # round-25 hardening: clocks nest, must not be summed; the snapshot
    # says so explicitly and names the key inclusive_seconds
    assert snap["clocks_are_inclusive"] is True
    assert "seconds" not in snap
    assert snap["inclusive_seconds"]["search_repair"] > 0
    kernel_profile_reset()
    assert kernel_profile_snapshot() == {"calls": {}, "inclusive_seconds": {},
                                         "counts": {},
                                         "clocks_are_inclusive": True}
    assert not KERNEL_PROFILE["calls"]
    # instrumented function still returns the unchanged value
    assert event_pairs_local(P, V, ev, 10.0, 0) == n0
