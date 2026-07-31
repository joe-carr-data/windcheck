"""Round-25 pins: exact remaining per-vertex budget intervals (round-24
Part C, approved spec) and the apply_field workspace base-identity token.

The seven reviewer-mandated interval cases:
  1. partial cancellation (a second repair that REDUCES |r| admits more
     than the scalar 1-|r| remaining budget would);
  2. orthogonal second displacement;
  3. vertex exactly on the budget boundary moving outward -> upper bound 0;
  4. asymmetric split + zero-weight vertices impose no constraint;
  5. float32 quantization at the analytic boundary (no ULP shrink; the
     emitted-float32 gate is the honest rejector -- recorded policy);
  6. empty intersection before any probe (candidate skipped, labelled);
  7. disjoint supports compose independently.

Plus: bitwise-default guarantees (P_orig omitted == today's behaviour),
certified_repair budget_point_vx plumbing, and the workspace token.
"""
from __future__ import annotations

import numpy as np
import pytest

from windcheck.check import PAIR_DTYPE
from windcheck.intrinsic import SurfaceGraph, oriented_events
from windcheck.repair import (BUDGET_ULP_ALLOWANCE_VX, apply_field,
                              budget_interval_from_arrays, certified_repair,
                              event_clearance, field_residual_arrays,
                              field_weights, remaining_budget_interval,
                              search_repair, split_hi_weight)


def planted():
    """The certified-repair fixture: quad (0,0) crosses quad (0,5)."""
    yz = [(0, -0.5), (1, -0.5), (2, -0.5), (3, -0.5),
          (3.5, 2.0), (1.5, 2.0), (1.5, -2.0)]
    P = np.array([[[v, y, z] for (y, z) in yz] for v in (0.0, 1.0)])
    V = np.ones((2, 7), dtype=bool)
    rec = np.array([(0, 1, 0, 5, 1.0, 10.0)], dtype=PAIR_DTYPE)
    ev = oriented_events(rec)[0]
    return P, V, ev


# ---------------------------------------------------------------------------
# 1-3: the per-vertex quadratic w^2 t^2 + 2 w (r.d) t + (|r|^2 - B^2) <= 0


def test_partial_cancellation_admits_more_than_scalar_remaining():
    B, w = 1.0, 0.5
    d = np.array([1.0, 0.0, 0.0])
    R = np.array([[-0.5, 0.0, 0.0]])          # moving +d REDUCES |r|
    iv = budget_interval_from_arrays(np.array([w]), R, d, B, w)
    assert iv is not None
    assert iv[0] == 0.0
    # displacement is (0.5 t - 0.5) x-hat: within budget until t = 3
    assert np.isclose(iv[1], 3.0)
    naive = (B - 0.5) / w                     # scalar 1-|r| remaining budget
    assert iv[1] > naive                      # cancellation admits MORE


def test_orthogonal_second_displacement():
    B, w = 1.0, 0.5
    d = np.array([1.0, 0.0, 0.0])
    R = np.array([[0.0, 0.6, 0.0]])           # residual orthogonal to d
    iv = budget_interval_from_arrays(np.array([w]), R, d, B, w)
    assert iv is not None and iv[0] == 0.0
    # |r + w t d|^2 = 0.36 + (0.5 t)^2 <= 1  ->  t <= sqrt(0.64)/0.5
    assert np.isclose(iv[1], np.sqrt(B * B - 0.36) / w)
    assert iv[1] < B / w                      # tighter than the naive cap


def test_boundary_vertex_moving_outward_upper_bound_zero():
    B, w = 1.0, 0.5
    d = np.array([1.0, 0.0, 0.0])
    R = np.array([[1.0, 0.0, 0.0]])           # ON the boundary, r || d
    iv = budget_interval_from_arrays(np.array([w]), R, d, B, w)
    assert iv == (0.0, 0.0)                   # any t > 0 leaves the budget


# ---------------------------------------------------------------------------
# 4: asymmetric split + zero-weight vertices


def test_asymmetric_split_and_zero_weight_vertices():
    P, V, ev = planted()
    gf = SurfaceGraph(np.asarray(P, np.float64), V, diagonal=-1,
                      maxedge=10.0)
    n_core, n_res, n_zero = gf.vertex(0, 0), gf.vertex(0, 1), gf.vertex(1, 6)
    w = {n_core: 0.3, n_res: -0.7, n_zero: 0.0}
    d = np.array([1.0, 0.0, 0.0])
    P_orig = np.asarray(P, np.float64).copy()
    P_cur = P_orig.copy()
    P_cur[0, 1] += np.array([0.0, 0.9, 0.0])  # orthogonal residual, w=-0.7
    P_cur[1, 6] += np.array([10.0, 0.0, 0.0])  # HUGE residual, w = 0
    iv = remaining_budget_interval(w, gf, P_cur, P_orig, d, 1.0, 0.3)
    # the zero-weight vertex is 10 vx outside the budget: were it counted,
    # the intersection would be empty -- it must impose no constraint
    assert iv is not None
    # binding constraint: the |w|=0.7 vertex with orthogonal residual 0.9
    assert np.isclose(iv[1], np.sqrt(1.0 - 0.81) / 0.7)
    # all-zero residuals: bitwise the historical scalar cap for lam=0.3
    iv0 = remaining_budget_interval(w, gf, P_orig, P_orig, d, 1.0, 0.3)
    assert iv0 == (0.0, 1.0 / max(0.3, 1.0 - 0.3))
    assert split_hi_weight(0.3) == max(0.3, 1.0 - 0.3)


# ---------------------------------------------------------------------------
# 5: float32 quantization at the analytic boundary (recorded policy)


def test_float32_quantization_at_analytic_boundary():
    # POLICY (recorded): the analytic root is NOT shrunk by an ULP
    # allowance; the emitted float32 cumulative gate rejects honestly.
    assert BUDGET_ULP_ALLOWANCE_VX == 0.0
    B, w = 1.0, 0.5
    x0 = 5000.0                               # float32-representable base
    disps = []
    for th in np.linspace(0.05, np.pi / 2 - 0.05, 64):
        d = np.array([np.cos(th), np.sin(th), 0.0])
        iv = budget_interval_from_arrays(np.array([w]), None, d, B, w)
        assert iv == (0.0, 2.0)               # analytic bound, unshrunk
        # at t = iv[1] the vertex moves w*t*d = exactly B along d; float32
        # emission rounds each component to the local grid, so the emitted
        # displacement can land an ULP EITHER side of B
        moved = np.float64(np.float32(x0 + w * iv[1] * d))
        disps.append(float(np.linalg.norm(moved - x0)))
    ulp = float(np.spacing(np.float32(5001.0)))
    # never more than quantization noise beyond the budget ...
    assert all(dv <= B + 2.0 * ulp for dv in disps)
    # ... but SOME bases do land beyond B: the analytic interval admitted
    # them and only the emitted-float32 gate rejects -- honestly
    assert any(dv > B for dv in disps)


# ---------------------------------------------------------------------------
# 6: empty intersection before any probe -> candidate skipped, labelled


def test_empty_intersection_candidate_skipped_and_labelled():
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    # every vertex already sits 1.5 vx (> B = 1) from its original: for
    # most directions the quadratic has no real root (disc < 0), for the
    # rest the per-vertex intersection over differing weights is empty
    P_orig = np.asarray(P, np.float64).copy()
    P_orig[..., 1] -= 1.5
    log: list = []
    out = search_repair(g, P, V, ev, budget_point_vx=1.0, support_vx=6.0,
                        P_orig=P_orig, skipped_log=log)
    assert out == []                          # nothing admissible was probed
    assert len(log) > 0                       # skips are LABELLED
    assert all(e["reason"] == "empty_remaining_budget_interval" for e in log)
    assert all("direction" in e and "lam" in e for e in log)
    # zero-residual run: historical silent skips are NOT labelled
    log2: list = []
    search_repair(g, P, V, ev, budget_point_vx=1.0, support_vx=6.0,
                  P_orig=P, skipped_log=log2)
    assert log2 == []


# ---------------------------------------------------------------------------
# 7: disjoint supports compose independently


def test_disjoint_supports_compose_independently():
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    gf = SurfaceGraph(np.asarray(P, np.float64), V, diagonal=-1,
                      maxedge=10.0)
    support = 1.5
    w = field_weights(gf, ev, 0.5, support)
    outside = [(v, u) for v in range(2) for u in range(7)
               if w.get(gf.vertex(v, u), 0.0) == 0.0]
    assert outside, "fixture must have vertices outside the field support"
    # an earlier repair moved ONLY vertices this event's field never
    # touches (disjoint supports)
    P_orig = np.asarray(P, np.float64).copy()
    for v, u in outside:
        P_orig[v, u] += np.array([0.9, 0.4, -0.7])
    d = np.array([0.0, 0.0, 1.0])
    iv_disjoint = remaining_budget_interval(w, gf, P, P_orig, d, 1.0, 0.5)
    iv_fresh = remaining_budget_interval(w, gf, P, None, d, 1.0, 0.5)
    assert iv_disjoint == iv_fresh == (0.0, 1.0 / 0.5)   # bitwise equal
    # search under the disjoint history == search from scratch, bitwise
    a = search_repair(g, P, V, ev, budget_point_vx=3.0, support_vx=support,
                      P_orig=P_orig)
    b = search_repair(g, P, V, ev, budget_point_vx=3.0, support_vx=support)
    assert len(a) == len(b)
    for ca, cb in zip(a, b):
        assert ca["direction"] == cb["direction"]
        assert ca["lam"] == cb["lam"]
        assert ca["t_rel"] == cb["t_rel"]
        assert np.array_equal(ca["P32"], cb["P32"])
    # contrast: a residual ON a nonzero-weight vertex DOES tighten
    inside = [(v, u) for v in range(2) for u in range(7)
              if abs(w.get(gf.vertex(v, u), 0.0)) == 0.5][0]
    P_orig2 = np.asarray(P, np.float64).copy()
    P_orig2[inside] += np.array([0.9, 0.0, 0.0])   # orthogonal to d, |r|=.9
    iv2 = remaining_budget_interval(w, gf, P, P_orig2, d, 1.0, 0.5)
    assert iv2 is not None and iv2[1] < iv_fresh[1]


# ---------------------------------------------------------------------------
# integration: bitwise-default guarantees and plumbing


def test_search_repair_explicit_p_orig_equals_default_bitwise():
    """transact now always passes P_orig; on a virgin base (r=0) the
    candidates must be BITWISE the historical ones -- the golden guard."""
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    plain = search_repair(g, P, V, ev, budget_point_vx=3.0, support_vx=6.0)
    withp = search_repair(g, P, V, ev, budget_point_vx=3.0, support_vx=6.0,
                          P_orig=P)
    assert len(plain) == len(withp) > 0
    for a, b in zip(plain, withp):
        assert a["direction"] == b["direction"]
        assert a["lam"] == b["lam"]
        assert a["t_rel"] == b["t_rel"]
        assert a["max_point_vx"] == b["max_point_vx"]
        assert np.array_equal(a["P32"], b["P32"])


def test_certified_repair_budget_point_interval():
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    gv, gu = np.nonzero(g.idx >= 0)
    X = np.empty((g.n, 3))
    X[g.idx[gv, gu]] = P[gv, gu]
    r = event_clearance(g, X, ev, 1000.0, t_max_vx=3.0)
    assert r is not None
    region = ev if r["side"] == "a" else \
        dict(ev, region_a=ev["region_b"], region_b=ev["region_a"])
    # budget_point_vx=None (default): the legacy path, still succeeds
    legacy = certified_repair(g, P, V, region, r["direction"], r["t_vx"],
                              mode="pinned", budget_vx=3.0, support_vx=6.0)
    assert legacy is not None
    assert "remaining_budget_interval_vx" not in legacy[1]
    # point budget below the rigid lower bound: empty intersection with
    # the candidate's rigid lower bound -> None before any probe
    tiny = certified_repair(g, P, V, region, r["direction"], r["t_vx"],
                            mode="pinned", budget_vx=3.0, support_vx=6.0,
                            P_orig=P, budget_point_vx=0.05 * r["t_vx"])
    assert tiny is None
    # generous budget: succeeds and records the interval + the policy
    ok = certified_repair(g, P, V, region, r["direction"], r["t_vx"],
                          mode="symmetric", budget_vx=3.0, support_vx=6.0,
                          P_orig=P, budget_point_vx=3.0)
    assert ok is not None
    rep = ok[1]
    assert rep["remaining_budget_interval_vx"] == [0.0, 3.0 / 0.5]
    assert rep["budget_interval_ulp_allowance_vx"] == 0.0


def test_zero_weight_vertices_excluded_from_residual_arrays():
    P, V, ev = planted()
    gf = SurfaceGraph(np.asarray(P, np.float64), V, diagonal=-1,
                      maxedge=10.0)
    w = {gf.vertex(0, 0): 0.5, gf.vertex(1, 6): 0.0}
    P_orig = np.asarray(P, np.float64).copy()
    P_orig[1, 6] += 100.0
    warr, R = field_residual_arrays(w, gf, P, P_orig)
    assert len(warr) == 1 and warr[0] == 0.5
    assert R.shape == (1, 3) and np.all(R == 0.0)


# ---------------------------------------------------------------------------
# workspace base-identity token (round-25 hardening)


def test_workspace_rejects_reuse_against_different_base_mesh():
    P, V, ev = planted()
    g = SurfaceGraph(P, V, diagonal=0, maxedge=10.0)
    gf = SurfaceGraph(np.asarray(P, np.float64), V, diagonal=-1,
                      maxedge=10.0)
    ws: dict = {}
    d = np.array([0.0, 0.0, 1.0])
    apply_field(g, P, V, ev, d, 0.5, "symmetric", 6.0, gf=gf, workspace=ws)
    assert "base_token" in ws
    # same shape, different content: rejected
    P2 = P.copy() + 0.25
    with pytest.raises(AssertionError, match="different base mesh"):
        apply_field(g, P2, V, ev, d, 0.5, "symmetric", 6.0, gf=gf,
                    workspace=ws)
    # same content, different object: the identity token still rejects
    # (a workspace is owned by ONE base mesh object, like gf/field_cache)
    P3 = P.copy()
    with pytest.raises(AssertionError, match="different base mesh"):
        apply_field(g, P3, V, ev, d, 0.5, "symmetric", 6.0, gf=gf,
                    workspace=ws)
    # the owning base keeps working after rejected attempts
    out, _ = apply_field(g, P, V, ev, d, 0.7, "symmetric", 6.0, gf=gf,
                         workspace=ws)
    fresh, _ = apply_field(g, P, V, ev, d, 0.7, "symmetric", 6.0, gf=gf)
    assert out.tobytes() == fresh.tobytes()
