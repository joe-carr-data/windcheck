"""Planted controls for rigid directional clearance.

Every asserted value is hand-computed from the construction. The epsilon
control is the round-13 mandate: a patch planted to penetrate by a known
epsilon must report clearance exactly epsilon (symmetric motion halves it,
which is the CALLER's convention to apply, not this module's).
"""
from __future__ import annotations

import numpy as np
import pytest

from windcheck.clearance import (clear_time, event_rigid_clearance,
                                 pair_interval)

BASE = np.array([[-2.0, -2.0, 0.0], [4.0, -2.0, 0.0], [0.0, 4.0, 0.0]])


def vertical(zmin: float, zmax: float) -> np.ndarray:
    """Vertical triangle through the base plane near the origin."""
    return np.array([[0.0, 0.0, zmin], [0.0, 0.0, zmax], [0.0, 1.0, 0.5 * (zmin + zmax)]])


def test_pair_interval_hand_computed():
    """A spans z in [-1, 2]: intersects while t in [-2, 1] along +z."""
    A = vertical(-1.0, 2.0)
    lo, hi = pair_interval(A, BASE, np.array([0.0, 0.0, 1.0]))
    assert lo == pytest.approx(-2.0, abs=1e-6)
    assert hi == pytest.approx(1.0, abs=1e-6)
    # reversed direction mirrors the interval
    lo2, hi2 = pair_interval(A, BASE, np.array([0.0, 0.0, -1.0]))
    assert (lo2, hi2) == (pytest.approx(-1.0, abs=1e-6),
                          pytest.approx(2.0, abs=1e-6))


def test_pair_interval_disjoint_is_none_or_positive():
    """A floating above the plane: no intersection at t=0."""
    A = vertical(0.5, 2.0)
    iv = pair_interval(A, BASE, np.array([0.0, 0.0, 1.0]))
    assert iv is not None
    lo, hi = iv
    assert lo < hi < 0.0 + 1e-9 or lo > 0.0   # interval excludes t=0
    assert not (lo <= 0.0 <= hi)


def test_planted_epsilon_control():
    """Penetration by exactly eps -> clearance along +z exactly eps."""
    eps = 0.05
    A = vertical(-eps, 1.0)
    r = event_rigid_clearance([A], [BASE],
                              [np.array([0.0, 0.0, 1.0])])
    assert r is not None
    assert r["t"] == pytest.approx(eps, abs=1e-4)


def test_direction_matters_and_best_is_chosen():
    """Deep along +z (t=1) but shallow along -z (t=2): +z wins; a sideways
    direction parallel to the base plane never separates the pair (the
    footprint stays inside the big base triangle for small t; the LP interval
    is wide), so the minimiser must be +z with t=1."""
    A = vertical(-1.0, 2.0)
    r = event_rigid_clearance(
        [A], [BASE],
        [np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, -1.0])])
    assert r is not None
    assert r["t"] == pytest.approx(1.0, abs=1e-4)
    assert r["direction"][2] == pytest.approx(1.0)


def test_neighbourhood_blocking_advances_t():
    """Clearing the crossing at t=1 lands inside a blocker spanning
    t in [0.8, 1.5]: the accepted translation must be 1.5, not 1."""
    assert clear_time(1.0, [(0.8, 1.5)]) == pytest.approx(1.5, abs=1e-6)
    assert clear_time(1.0, [(2.0, 3.0)]) == pytest.approx(1.0)
    # chained blockers: [0.9,1.4] then [1.35,1.8] -> 1.8
    assert clear_time(1.0, [(0.9, 1.4), (1.35, 1.8)]) == pytest.approx(1.8, abs=1e-6)


def test_event_clearance_with_blocking_triangle():
    """Full event path: crossing needs t=1 along +z, a static ceiling
    triangle occupies the band z in [1.0, 1.5] over the origin -- the moving
    patch (z span 3) intersects it for t in [-2.5, 2.5+eps...]; hand
    computation: moving A spans z in [-1,2]+t; ceiling at z in [1.0,1.5];
    overlap while -1+t <= 1.5 and 2+t >= 1.0 -> t in [-1, 2.5]. Required
    t=1 sits inside -> advance to 2.5."""
    A = vertical(-1.0, 2.0)
    ceiling = np.array([[-2.0, -2.0, 1.25], [4.0, -2.0, 1.25],
                        [0.0, 4.0, 1.25]])
    r = event_rigid_clearance([A], [BASE], [np.array([0.0, 0.0, 1.0])],
                              neighbourhood=[(A, ceiling)])
    assert r is not None
    # ceiling is a plane at z=1.25: A intersects it while -1+t<=1.25<=2+t
    # -> t in [-0.75, 2.25]; required 1 -> advanced to 2.25
    assert r["t"] == pytest.approx(2.25, abs=1e-4)


def test_swept_neighbourhood_sees_grid_distant_blocker():
    """Round-14 negative control 4: a shelf 7 grid columns away but 0.8 vx
    above the moving patch must enter the obstacle set when the budget can
    reach it (t_max=2) and stay out when it cannot (t_max=0.1)."""
    import numpy as np
    from windcheck.intrinsic import SurfaceGraph
    from windcheck.repair import region_triangles, swept_neighbourhood

    # strip: flat base at z=-0.5 (u=0..3), excursion, then a shelf at z=0.3
    # hovering over the base run -- grid-distant, spatially adjacent
    yz = [(0, -0.5), (1, -0.5), (2, -0.5), (3, -0.5),
          (6, -0.5), (6, 4.0), (0.5, 4.0), (0.5, 0.3), (2.5, 0.3)]
    P = np.array([[[v, y, z] for (y, z) in yz] for v in (0.0, 1.0)])
    V = np.ones((2, 9), dtype=bool)
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    moving = [(0, 1)]                       # base quad under the shelf
    X = _X(g, P)
    mov = region_triangles(g, X, moving)
    far = swept_neighbourhood(g, X, mov, moving, [], t_max=2.0)
    near = swept_neighbourhood(g, X, mov, moving, [], t_max=0.1)
    def has_shelf(tris):
        return any(np.any(np.isclose(T[:, 2], 0.3)) for T in tris)
    assert has_shelf(far)                   # budget reaches the shelf
    assert not has_shelf(near)              # tight budget excludes it


def _X(g, P):
    import numpy as np
    vv, uu = np.nonzero(g.idx >= 0)
    X = np.empty((g.n, 3))
    X[g.idx[vv, uu]] = P[vv, uu]
    return X




def test_certified_repair_removes_planted_crossing_both_modes():
    """Round-15 end-to-end: clearance -> harmonic field -> float32-certified
    removal. Pinned mode must leave every partner vertex unmoved; symmetric
    must split the motion; both must clear the crossing on the QUANTIZED
    mesh with clean gates."""
    import numpy as np
    from windcheck.intrinsic import SurfaceGraph, oriented_events
    from windcheck.repair import (certified_repair, event_clearance,
                                  event_crossing_pairs)
    from windcheck.check import PAIR_DTYPE

    yz = [(0, -0.5), (1, -0.5), (2, -0.5), (3, -0.5),
          (3.5, 2.0), (1.5, 2.0), (1.5, -2.0)]
    P = np.array([[[v, y, z] for (y, z) in yz] for v in (0.0, 1.0)])
    V = np.ones((2, 7), dtype=bool)
    rec = np.array([(0, 1, 0, 5, 1.0, 10.0)], dtype=PAIR_DTYPE)
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    X = _X(g, P)
    ev = oriented_events(rec)[0]
    assert event_crossing_pairs(g, X, ev) > 0          # planted and present
    r = event_clearance(g, X, ev, 1000.0, t_max_vx=3.0)
    assert r is not None
    region = ev if r["side"] == "a" else \
        dict(ev, region_a=ev["region_b"], region_b=ev["region_a"])
    for mode in ("pinned", "symmetric"):
        out = certified_repair(g, P, V, region, r["direction"], r["t_vx"],
                               mode=mode, budget_vx=3.0, support_vx=6.0)
        assert out is not None, mode
        P32, rep = out
        g2 = SurfaceGraph(P32, V, diagonal=0, maxedge=10.0)
        assert event_crossing_pairs(g2, _X(g2, P32), region) == 0
        assert rep["triangle_inversions_d0"] == 0
        assert rep["triangle_inversions_d1"] == 0
        assert rep["quads_newly_dropped"] == 0
        assert rep["additional_displacement_to_local_clearance_vx"] >= 0.0
        assert rep["bracket_clean_vx"] >= rep["bracket_dirty_vx"]
        assert rep["bracket_clean_vx"] == rep["applied_relative_vx"]
        if mode == "pinned":
            assert rep["field_partner_nonzero"] == 0
            assert rep["quantized_partner_moved"] == 0
        else:
            assert rep["quantized_partner_moved"] > 0
            assert rep["max_disp_vx"] <= 0.51 * rep["applied_relative_vx"]


def test_harmonic_field_is_intrinsic_not_grid():
    """Round-15 hole 2 pin: an island two grid columns from the core but
    with NO surface path to it must receive zero weight."""
    import numpy as np
    from windcheck.intrinsic import SurfaceGraph
    from windcheck.repair import harmonic_field

    P, V = np.zeros((6, 8, 3)), np.ones((6, 8), bool)
    v, u = np.mgrid[0:6, 0:8]
    P[..., 0], P[..., 1] = u, v
    V[:, 3] = False                       # cuts the sheet: island at u>=4
    g = SurfaceGraph(P, V)
    core = {g.vertex(2, 1), g.vertex(3, 1)}
    w = harmonic_field(g, core, set(), 1.0, 0.0, support_vx=20.0)
    island = [g.vertex(2, 5), g.vertex(3, 6)]
    assert all(w.get(n, 0.0) == 0.0 for n in island)
    near = g.vertex(2, 2)                 # connected: positive weight
    assert w.get(near, 0.0) > 0.0


def test_float32_input_at_scroll_magnitude():
    """Round-16 blocker 2 pin: real tifxyz is float32 near coordinate 5e3.
    Deformation must be float64 internally (P.copy() on float32 input kept
    all arithmetic in float32) and still certify after quantization."""
    import numpy as np
    from windcheck.intrinsic import SurfaceGraph, oriented_events
    from windcheck.repair import (apply_field, certified_repair,
                                  event_crossing_pairs)
    from windcheck.check import PAIR_DTYPE

    yz = [(0, -0.5), (1, -0.5), (2, -0.5), (3, -0.5),
          (3.5, 2.0), (1.5, 2.0), (1.5, -2.0)]
    P = np.array([[[v, y, z] for (y, z) in yz] for v in (0.0, 1.0)])
    P += np.array([5000.0, 5000.0, 5000.0])          # scroll magnitude
    P32in = P.astype(np.float32)                     # the real input dtype
    V = np.ones((2, 7), dtype=bool)
    rec = np.array([(0, 1, 0, 5, 1.0, 10.0)], dtype=PAIR_DTYPE)
    g = SurfaceGraph(P32in, V, diagonal=0, maxedge=10.0)
    ev = oriented_events(rec)[0]
    P2, _ = apply_field(g, P32in, V, ev, [0.0, 0.0, 1.0], 0.5)
    assert P2.dtype == np.float64                    # float64 pipeline
    out = certified_repair(g, P32in, V, ev, [0.0, 0.0, 1.0], 1.5,
                           mode="symmetric", budget_vx=3.0, support_vx=6.0)
    assert out is not None
    P32, rep = out
    g2 = SurfaceGraph(P32, V, diagonal=0, maxedge=10.0)
    assert event_crossing_pairs(g2, _X(g2, P32), ev) == 0
    assert rep["numerical_resolution_vx"] >= 4 * 4.8e-4  # ULP at 5e3


def test_merge_events_orients_the_union():
    """d0/d1 regions of one event unify with consistent orientation."""
    from windcheck.repair import merge_events
    ev0 = {"region_a": {(504, 130), (505, 130)},
           "region_b": {(504, 132), (505, 132)}}
    aligned = {"region_a": {(503, 130)}, "region_b": {(503, 132)}}
    m = merge_events(ev0, aligned)
    assert (503, 130) in m["region_a"] and (503, 132) in m["region_b"]
    swapped = {"region_a": {(503, 132)}, "region_b": {(503, 130)}}
    m2 = merge_events(ev0, swapped)
    assert (503, 130) in m2["region_a"] and (503, 132) in m2["region_b"]


def test_pair_min_exit_planted_epsilon():
    """Planted penetration eps: the cheapest exit anywhere is eps (up)."""
    from windcheck.clearance import pair_min_exit, event_min_exit
    eps = 0.05
    A = vertical(-eps, 1.0)
    assert pair_min_exit(A, BASE) == pytest.approx(eps, abs=1e-6)
    assert event_min_exit([A], [BASE]) == pytest.approx(eps, abs=1e-6)
    # separated pair: 0 (nothing to certify)
    assert pair_min_exit(vertical(0.5, 2.0), BASE) == 0.0


def test_match_events_outcomes():
    """Matcher yields matched / unmatched / ambiguous explicitly."""
    from windcheck.repair import match_events
    A = {"region_a": {(10, 10)}, "region_b": {(10, 40)}}
    A1s = {"region_a": {(11, 10)}, "region_b": {(11, 40)}}   # matches A
    B = {"region_a": {(50, 5)}, "region_b": {(50, 90)}}       # d0 only
    C1 = {"region_a": {(80, 8)}, "region_b": {(80, 60)}}      # d1 only
    r = match_events([A, B], [A1s, C1])
    assert r["matched"] and r["matched"][0][:2] == (0, 0)
    assert r["unmatched_d0"] == [1] and r["unmatched_d1"] == [1]
    assert not r["ambiguous"]
    # two near-identical d1 partners -> ambiguous, not guessed
    r2 = match_events([A], [A1s, {"region_a": {(9, 10)},
                                  "region_b": {(9, 40)}}])
    assert r2["ambiguous"] and not r2["matched"]


def test_match_events_reciprocal_ownership():
    """Round-23 blocker 2: a d1 event courted equally by two d0 events
    must land in exactly ONE ambiguous group -- never both matched and
    rival -- and every event index belongs to exactly one bucket."""
    from windcheck.repair import match_events
    X = {"region_a": {(10, 10)}, "region_b": {(10, 40)}}
    Z = {"region_a": {(10, 11)}, "region_b": {(10, 41)}}
    Y = {"region_a": {(11, 10)}, "region_b": {(11, 40)}}  # near X and Z
    r = match_events([X, Z], [Y])
    assert not r["matched"]
    assert len(r["ambiguous"]) == 1
    g0, g1 = r["ambiguous"][0]
    assert set(g0) == {0, 1} and g1 == [0]
    assert r["unmatched_d0"] == [] and r["unmatched_d1"] == []
