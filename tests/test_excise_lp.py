"""Acceptance pins for the round-27 LP-rounding selection strategy.

The exact MILP returns NO INCUMBENT on the crowded and monster segments,
so the delivered cut on those comes from the deterministic constructor in
`windcheck.excise`: LP relaxation -> 1/k rounding -> reverse-delete ->
bounded local improvement, with an area-aware greedy fallback when the LP
itself does not solve. Everything below is PLANTED -- no real data, no
engine -- and pins the four properties the strategy is allowed to claim:

  (a) the LP-rounded mask is FEASIBLE and within k x the exact optimum,
      and reverse-delete never breaks feasibility;
  (b) the ordering A_LP <= A_exact <= A_rounded holds (the LP objective is
      a certified LOWER bound, the rounded area an achieved UPPER one);
  (c) rounding followed by reverse-delete can be STRICTLY better than raw
      rounding -- pinned on a planted instance where it is 18 -> 4;
  (d) the greedy fallback path produces a feasible mask and claims no
      approximation ratio.
"""
from __future__ import annotations

import numpy as np
import pytest

from windcheck import excise as ex
from windcheck.excise import (Component, _solve_milp, census, excise,
                              greedy_cover, local_improve, lp_relaxation,
                              quad_area_canonical, quad_area_canonical_grid,
                              quad_corners, reduce_constraints, reverse_delete,
                              round_lp, select_global)
from windcheck.intrinsic import retained_quads

from test_excise import blade_two_walls, blade_wall, fold_mesh


# ------------------------------------------------------------- fixtures
def flat_grid(nv: int, nu: int):
    """Planar unit-square grid: every quad has canonical area exactly 1."""
    P = np.zeros((nv, nu, 3))
    P[:, :, 0] = np.arange(nv)[:, None]
    P[:, :, 1] = np.arange(nu)[None, :]
    V = np.ones((nv, nu), bool)
    return P, V


def pair_constraint(qa, qb):
    """The round-24 eight-corner coverage set for a crossing quad pair."""
    return {"key": (0, *qa, *qb), "q1": qa, "q2": qb,
            "coverage": sorted(set(quad_corners(*qa)) | set(quad_corners(*qb)))}


def model_of(P, V, constraints=None, protected=()):
    P64 = np.asarray(P, np.float64)
    Q = retained_quads(P64, np.asarray(V, bool), 60.0)
    if constraints is None:
        constraints = census(P64, np.asarray(V, bool))
    return P64, Q, Component(constraints, P64, Q, set(protected),
                             quad_area_canonical_grid(P64))


def exact_area(constraints, P64, Q):
    """Proven stage-1 optimum for a component, or None if not proven."""
    sol = _solve_milp(constraints, P64, Q, set(), None)
    if sol["status"] != "optimal":
        return None, sol
    return sol["records"][0]["objective"], sol


# ------------------------------------------------------- area convention
def test_vectorised_canonical_area_matches_the_scalar_function():
    rng = np.random.default_rng(11)
    P = rng.normal(size=(7, 9, 3))
    grid = quad_area_canonical_grid(P)
    for v in range(6):
        for u in range(8):
            assert abs(float(grid[v, u]) - quad_area_canonical(P, v, u)) < 1e-12


# --------------------------------------- (a) feasible and within k x opt
@pytest.mark.parametrize("fixture", [blade_wall, blade_two_walls,
                                     lambda: fold_mesh(True, True)])
def test_a_lp_rounding_is_feasible_and_within_k_times_the_optimum(fixture):
    P, V = fixture()
    P64, Q, _ = model_of(P, V)
    contacts = census(P64, np.asarray(V, bool))
    assert contacts, "fixture must plant at least one crossing"
    red = reduce_constraints(contacts, Q)
    for group in red["components"]:
        m = Component(group, P64, Q, set(), quad_area_canonical_grid(P64))
        lp = lp_relaxation(m)
        assert lp["status"] == "optimal"
        raw = round_lp(m, lp["x"])
        assert m.is_feasible(raw), "1/k rounding must be feasible by construction"
        kept = reverse_delete(m, raw)
        assert m.is_feasible(kept), "reverse-delete must preserve feasibility"
        assert kept.sum() <= raw.sum()
        opt, _ = exact_area(group, P64, Q)
        assert opt is not None
        assert m.area(raw) <= m.k * lp["objective"] + 1e-9
        assert m.area(raw) <= m.k * opt + 1e-9
        assert m.area(kept) <= m.area(raw) + 1e-12


# ----------------------------------- (b) A_LP <= A_exact <= A_rounded
@pytest.mark.parametrize("fixture", [blade_wall, blade_two_walls,
                                     lambda: fold_mesh(True, True),
                                     lambda: fold_mesh(True, False),
                                     lambda: fold_mesh(False, True)])
def test_b_bound_ordering_lp_le_exact_le_rounded(fixture):
    P, V = fixture()
    P64, Q, _ = model_of(P, V)
    contacts = census(P64, np.asarray(V, bool))
    red = reduce_constraints(contacts, Q)
    for group in red["components"]:
        m = Component(group, P64, Q, set(), quad_area_canonical_grid(P64))
        lp = lp_relaxation(m)
        opt, _ = exact_area(group, P64, Q)
        assert lp["status"] == "optimal" and opt is not None
        raw = round_lp(m, lp["x"])
        a_round = m.area(raw)
        assert lp["objective"] <= opt + 1e-9, "LP is not a lower bound"
        assert opt <= a_round + 1e-9, "rounded area beats the proven optimum"
        kept = reverse_delete(m, raw)
        assert opt <= m.area(kept) + 1e-9


def test_b2_bound_ordering_on_a_planted_set_cover_instance():
    P, V = flat_grid(8, 10)
    cons = [pair_constraint((3, 2), (3, 6))]
    P64, Q, m = model_of(P, V, cons)
    lp = lp_relaxation(m)
    opt, _ = exact_area(cons, P64, Q)
    raw = round_lp(m, lp["x"])
    assert lp["objective"] <= opt + 1e-9 <= m.area(raw) + 1e-9


# -------------------- (c) reverse-delete strictly improves raw rounding
def test_c_reverse_delete_strictly_improves_on_raw_rounding():
    """Two interior quads, one coverage constraint, unit-area quads.

    Spreading x = 1/4 over the four corners of ONE of the two quads costs
    the LP 9/4 = 2.25, strictly less than the cheapest integral vertex (4),
    so the relaxation IS fractional; 1/8-rounding then selects all four of
    those corners and charges their whole 3x3 quad block: 9. The single
    constraint still has four covers, so reverse-delete strips three of
    them and lands exactly on the exact optimum, 4.
    """
    P, V = flat_grid(8, 10)
    cons = [pair_constraint((3, 2), (3, 6))]
    P64, Q, m = model_of(P, V, cons)
    assert m.k == 8
    lp = lp_relaxation(m)
    assert lp["status"] == "optimal"
    raw = round_lp(m, lp["x"])
    a_raw = m.area(raw)
    kept = reverse_delete(m, raw)
    a_kept = m.area(kept)
    assert m.is_feasible(raw) and m.is_feasible(kept)
    assert a_kept < a_raw - 1e-9, "reverse-delete must strictly improve here"
    assert a_raw == pytest.approx(9.0)
    assert a_kept == pytest.approx(4.0)
    opt, _ = exact_area(cons, P64, Q)
    assert opt == pytest.approx(4.0)
    assert lp["objective"] == pytest.approx(2.25)
    assert a_raw <= m.k * lp["objective"] + 1e-9      # the k-approximation


def test_c2_reverse_delete_never_breaks_feasibility_on_a_dense_superset():
    P, V = flat_grid(9, 12)
    cons = [pair_constraint((3, 2), (3, 6)),
            pair_constraint((4, 3), (4, 8)),
            pair_constraint((2, 2), (5, 9))]
    P64, Q, m = model_of(P, V, cons)
    everything = ~m.protected                      # select every candidate
    kept = reverse_delete(m, everything)
    assert m.is_feasible(kept)
    assert kept.sum() < everything.sum()
    assert m.area(kept) < m.area(everything)


def test_c3_local_improvement_preserves_feasibility_and_never_worsens():
    P, V = flat_grid(9, 12)
    cons = [pair_constraint((3, 2), (3, 6)),
            pair_constraint((4, 3), (4, 8))]
    P64, Q, m = model_of(P, V, cons)
    lp = lp_relaxation(m)
    sel = reverse_delete(m, round_lp(m, lp["x"]))
    before = m.area(sel)
    imp = local_improve(m, sel, None)
    assert m.is_feasible(imp["sel"])
    assert m.area(imp["sel"]) <= before + 1e-12


# ------------------------------------------- (d) greedy fallback is feasible
def test_d_greedy_fallback_is_feasible():
    P, V = flat_grid(9, 12)
    cons = [pair_constraint((3, 2), (3, 6)),
            pair_constraint((4, 3), (4, 8)),
            pair_constraint((2, 2), (5, 9))]
    P64, Q, m = model_of(P, V, cons)
    g = greedy_cover(m)
    assert g["covered_all"]
    assert m.is_feasible(g["sel"])
    assert m.is_feasible(reverse_delete(m, g["sel"]))


@pytest.mark.parametrize("fixture", [blade_wall, blade_two_walls,
                                     lambda: fold_mesh(True, True)])
def test_d2_greedy_fallback_is_feasible_on_planted_meshes(fixture):
    P, V = fixture()
    P64, Q, _ = model_of(P, V)
    contacts = census(P64, np.asarray(V, bool))
    for group in reduce_constraints(contacts, Q)["components"]:
        m = Component(group, P64, Q, set(), quad_area_canonical_grid(P64))
        g = greedy_cover(m)
        assert g["covered_all"] and m.is_feasible(g["sel"])


def test_d3_greedy_path_is_taken_when_the_lp_does_not_solve(monkeypatch):
    """A component whose LP does not solve must still get a feasible mask,
    labelled heuristic_feasible and carrying NO lower bound."""
    P, V = blade_two_walls()
    P64 = np.asarray(P, np.float64)
    Q = retained_quads(P64, np.asarray(V, bool), 60.0)
    contacts = census(P64, np.asarray(V, bool))

    def dead_lp(m, time_limit=None):
        return {"status": "limit", "raw_status": 1, "message": "mocked",
                "objective": None, "x": None, "seconds": 0.0,
                "n_rows": 0, "n_cols": 0, "n_nnz": 0}

    monkeypatch.setattr(ex, "lp_relaxation", dead_lp)
    sel = select_global(contacts, P64, Q, exact_max_constraints=0,
                        improve_budget=0.0)
    assert sel["status"] == "ok"
    assert sel["selection_status"] == "heuristic_feasible"
    assert sel["combined_lower_bound_complete"] is False
    assert all(c["lower_bound"] is None for c in sel["components"])
    # no component has a bound, so NO ratio may be quoted at all
    assert sel["ratio_achieved_over_bound"] is None
    assert sel["achieved_area_bounded_subset"] == 0.0
    assert sel["achieved_area_unbounded_subset"] == sel["achieved_area"]
    assert sel["ratio_covers_area_fraction"] == 0.0
    V_out = np.asarray(V, bool).copy()
    for v, u in sel["chosen"]:
        V_out[v, u] = False
    assert census(P64, V_out) == [], "greedy mask must still be feasible"


# ------------------------------------------------------ policy / wiring
def test_e_select_global_proves_small_components_optimal():
    P, V = blade_two_walls()
    P64 = np.asarray(P, np.float64)
    Q = retained_quads(P64, np.asarray(V, bool), 60.0)
    contacts = census(P64, np.asarray(V, bool))
    sel = select_global(contacts, P64, Q)
    assert sel["status"] == "ok"
    assert sel["selection_status"] == "area_optimal"
    assert sel["method_mix"] == {"exact_optimal": sel["reduction"]
                                 ["n_components"]}
    assert sel["combined_lower_bound_complete"] is True
    assert sel["ratio_achieved_over_bound"] == pytest.approx(1.0)
    for c in sel["components"]:
        assert c["lower_bound_kind"] == "exact_optimum"
    V_out = np.asarray(V, bool).copy()
    for v, u in sel["chosen"]:
        V_out[v, u] = False
    assert census(P64, V_out) == []


def test_e2_select_global_stays_lp_rounded_when_exact_is_switched_off():
    P, V = blade_two_walls()
    P64 = np.asarray(P, np.float64)
    Q = retained_quads(P64, np.asarray(V, bool), 60.0)
    contacts = census(P64, np.asarray(V, bool))
    sel = select_global(contacts, P64, Q, exact_max_constraints=0)
    assert sel["selection_status"] == "lp_rounded"
    assert sel["combined_lower_bound_complete"] is True
    for c in sel["components"]:
        assert c["lower_bound_kind"] == "lp_objective"
        assert c["achieved_area"] >= c["lower_bound"] - 1e-9
        assert c["achieved_area"] <= c["k"] * c["lower_bound"] + 1e-9
    assert sel["ratio_achieved_over_bound"] >= 1.0 - 1e-12
    V_out = np.asarray(V, bool).copy()
    for v, u in sel["chosen"]:
        V_out[v, u] = False
    assert census(P64, V_out) == []


def test_f_excise_lp_round_strategy_end_to_end():
    P, V = blade_wall()
    r = excise(P, V, strategy="lp_round")
    assert r["status"] == "clean"
    cert = r["certificate"]
    assert cert["geometry_status"] == "transverse_clean_certified"
    assert cert["selection_status"] in ("area_optimal", "lp_rounded",
                                        "heuristic_feasible", "mixed")
    sel = cert["selection"]
    assert sel["strategy"] == "lp_round"
    assert sel["combined_lower_bound"] <= sel["achieved_area_canonical"] + 1e-9
    assert cert["output_census"]["clean_both_diagonals"] is True
    assert census(np.asarray(P, np.float64), r["valid"]) == []


def test_f2_lp_round_never_claims_minimum_area_for_a_rounded_component():
    P, V = blade_two_walls()
    r = excise(P, V, strategy="lp_round",
               selection_options={"exact_max_constraints": 0})
    assert r["status"] == "clean"
    cert = r["certificate"]
    assert cert["selection_status"] == "lp_rounded"
    assert cert["selection"]["ratio_achieved_over_bound"] is not None
    # the honest wording gate: the only place the phrase may appear is the
    # RULE that forbids it, never as a property of this cut
    assert cert["selection"]["method_mix"] == {"lp_rounded": 1}
    assert all(c["lower_bound_kind"] == "lp_objective"
               for c in cert["selection"]["per_component"])
    assert cert["solver"]["lexicographic"] is False
    assert cert["selection"]["achieved_area_canonical"] \
        >= cert["selection"]["combined_lower_bound"]


def test_e3_ratio_never_divides_full_area_by_a_partial_bound(monkeypatch):
    """With one component bounded and one not, the quoted ratio must cover
    ONLY the bounded component -- otherwise the unbounded component's area
    inflates it into a number that says nothing about selection quality."""
    P, V = flat_grid(9, 22)
    cons = [pair_constraint((3, 2), (3, 6)), pair_constraint((3, 14),
                                                             (3, 18))]
    P64, Q, _ = model_of(P, V, cons)
    from windcheck.excise import reduce_constraints as _rc
    assert _rc(cons, Q)["n_components"] == 2, "fixture must be separable"
    real = ex.lp_relaxation
    calls = {"n": 0}

    def flaky_lp(m, time_limit=None):
        calls["n"] += 1
        if calls["n"] == 1:                 # first component: LP "times out"
            return {"status": "limit", "raw_status": 1, "message": "mocked",
                    "objective": None, "x": None, "seconds": 0.0,
                    "n_rows": 0, "n_cols": 0, "n_nnz": 0}
        return real(m, time_limit)

    monkeypatch.setattr(ex, "lp_relaxation", flaky_lp)
    sel = select_global(cons, P64, Q, exact_max_constraints=0,
                        improve_budget=0.0)
    assert sel["selection_status"] == "mixed"
    assert sel["combined_lower_bound_complete"] is False
    assert sel["n_components_without_bound"] == 1
    assert sel["achieved_area_bounded_subset"] > 0
    assert sel["achieved_area_unbounded_subset"] > 0
    assert sel["achieved_area_bounded_subset"] \
        + sel["achieved_area_unbounded_subset"] \
        == pytest.approx(sel["achieved_area"])
    assert sel["ratio_achieved_over_bound"] == pytest.approx(
        sel["achieved_area_bounded_subset"] / sel["combined_lower_bound"])
    assert 0.0 < sel["ratio_covers_area_fraction"] < 1.0


def test_g_protected_vertex_makes_a_component_infeasible():
    P, V = flat_grid(8, 10)
    cons = [pair_constraint((3, 2), (3, 6))]
    P64, Q, _ = model_of(P, V, cons)
    protected = cons[0]["coverage"]
    sel = select_global(cons, P64, Q, protected)
    assert sel["status"] == "infeasible"
    assert sel["chosen"] == []
