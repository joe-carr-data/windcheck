"""Acceptance pins for the ROUND-28 FROZEN greedy-first policy.

The round-27 prototype ran the LP first and spent 600-1200 s on the giant
components for neither a mask nor a bound. The frozen policy inverts the
schedule, and the properties that make that legitimate engineering rather
than claim tuning are exactly the ones pinned here:

  (a) GREEDY-FIRST ORDERING -- every component holds a feasible incumbent
      BEFORE any solver is called, so a segment with all budgets set to
      zero still emits a complete feasible mask;
  (b) AN OPTIMIZATION FAILURE NEVER REMOVES THE INCUMBENT -- the LP, the
      exact MILP and the local improvement are each made to explode, and
      the returned mask is still feasible and still covers every
      constraint;
  (c) REPLACEMENT ONLY ON STRICTLY LOWER MEASURED AREA -- a candidate that
      is not strictly better leaves the incumbent in place;
  (d) THE 99.9%-AREA CORE -- the input-only fragmentation gate, including
      the loophole the round-27 ">=1% of segment area" rule had;
  (e) the policy hash is stable, and it changes when the policy changes;
  (f) no minimum-area claim unless EVERY component is proven optimal.

Everything is PLANTED: no real data, no engine.
"""
from __future__ import annotations

import sys

import numpy as np
import pytest

sys.path.insert(0, "bench")

from windcheck import excise as ex                              # noqa: E402
from windcheck.excise import (FROZEN_POLICY, FROZEN_POLICY_VERSION,  # noqa: E402
                              census, frozen_policy_hash,
                              quad_area_canonical_grid,
                              select_global_frozen)
from windcheck.intrinsic import retained_quads                  # noqa: E402

from test_excise import blade_two_walls, blade_wall, fold_mesh   # noqa: E402
from excise_shadow import area_core                              # noqa: E402


ZERO_BUDGETS = dict(FROZEN_POLICY,
                    lp_budget_s_per_segment=0.0,
                    improvement_budget_s_per_segment=0.0)


def setup(fixture):
    P, V = fixture()
    P64 = np.asarray(P, np.float64)
    Vb = np.asarray(V, bool)
    Q = retained_quads(P64, Vb, 60.0)
    cons = census(P64, Vb)
    return P64, Vb, Q, cons


def covered(cons, chosen) -> bool:
    sel = {tuple(map(int, v)) for v in chosen}
    return all(sel & {tuple(map(int, p)) for p in c["coverage"]}
               for c in cons)


FIXTURES = [blade_wall, blade_two_walls, lambda: fold_mesh(True, True)]


# ------------------------------------------- (a) greedy-first ordering
@pytest.mark.parametrize("fixture", FIXTURES)
def test_every_component_is_feasible_before_any_solver_runs(fixture):
    """With BOTH optional budgets at zero, no LP and no MILP may run, and
    the policy must still return a complete feasible mask."""
    P64, V, Q, cons = setup(fixture)
    calls = []
    orig_lp, orig_milp = ex.lp_relaxation, ex._solve_milp
    ex.lp_relaxation = lambda *a, **k: calls.append("lp")
    ex._solve_milp = lambda *a, **k: calls.append("milp")
    try:
        out = select_global_frozen(cons, P64, Q, (), policy=ZERO_BUDGETS,
                                   area_grid=quad_area_canonical_grid(P64))
    finally:
        ex.lp_relaxation, ex._solve_milp = orig_lp, orig_milp
    assert calls == [], "a solver ran before/despite the zero budgets"
    assert out["status"] == "ok"
    assert covered(cons, out["chosen"])
    assert set(out["method_mix"]) == {"greedy_feasible"}
    assert out["selection_status"] == "mixed"
    assert out["minimum_area_claim_admissible"] is False
    for c in out["components"]:
        assert c["greedy_incumbent_area"] >= 0.0
        assert c["greedy_construction_seconds"] >= 0.0
        assert c["achieved_area"] <= c["greedy_incumbent_area"] + 1e-12


@pytest.mark.parametrize("fixture", FIXTURES)
def test_greedy_incumbent_is_recorded_before_improvement(fixture):
    """The greedy incumbent area is a RECORD, not a leftover: it must be
    present on every component even when a later phase improved on it."""
    P64, V, Q, cons = setup(fixture)
    out = select_global_frozen(cons, P64, Q, (),
                               area_grid=quad_area_canonical_grid(P64))
    assert out["greedy_incumbent_area"] >= out["achieved_area"] - 1e-12
    assert out["improvement_over_greedy"] >= -1e-12
    for c in out["components"]:
        assert "greedy_incumbent_area" in c
        assert "greedy_construction_seconds" in c
        assert c["method"] in ("exact_optimal", "lp_improved",
                               "greedy_feasible")
        assert c["lp_attempted"] or "lp_skipped_reason" in c
        assert c["exact_attempted"] or "exact_skipped_reason" in c


# ------------------ (b) an optimization failure never removes the mask
@pytest.mark.parametrize("blow_up", ["lp", "milp", "improve", "all"])
def test_optimization_failure_cannot_remove_the_incumbent(blow_up):
    P64, V, Q, cons = setup(blade_two_walls)
    base = select_global_frozen(cons, P64, Q, (),
                                policy=ZERO_BUDGETS,
                                area_grid=quad_area_canonical_grid(P64))

    def boom(*a, **k):
        raise RuntimeError("planted solver failure")

    saved = (ex.lp_relaxation, ex._solve_milp, ex.local_improve)
    if blow_up in ("lp", "all"):
        ex.lp_relaxation = boom
    if blow_up in ("milp", "all"):
        ex._solve_milp = boom
    if blow_up in ("improve", "all"):
        ex.local_improve = boom
    try:
        out = select_global_frozen(cons, P64, Q, (),
                                   area_grid=quad_area_canonical_grid(P64))
    finally:
        ex.lp_relaxation, ex._solve_milp, ex.local_improve = saved

    assert out["status"] == "ok", "a solver failure destroyed the run"
    assert covered(cons, out["chosen"]), "the emitted mask is not feasible"
    # the incumbent survived: never worse than the greedy-only selection
    assert out["achieved_area"] <= base["achieved_area"] + 1e-9
    if blow_up in ("lp", "all"):
        assert any(c.get("lp_status") == "error" for c in out["components"])
    if blow_up in ("milp", "all"):
        assert any(c.get("exact_status") == "error"
                   for c in out["components"])


def test_an_infeasible_rounding_leaves_the_incumbent_alone():
    """A rounded candidate that fails its feasibility re-check is recorded
    and DISCARDED -- it must never be installed."""
    P64, V, Q, cons = setup(blade_two_walls)
    saved = ex.round_lp
    ex.round_lp = lambda m, x: np.zeros(len(m.verts), bool)   # infeasible
    try:
        out = select_global_frozen(cons, P64, Q, (),
                                   area_grid=quad_area_canonical_grid(P64))
    finally:
        ex.round_lp = saved
    assert out["status"] == "ok"
    assert covered(cons, out["chosen"])
    assert any(c.get("rounding_infeasible") for c in out["components"])
    assert not any(c["method"] == "lp_improved" for c in out["components"])


# --------------------- (c) replacement only on strictly lower measured area
def test_a_worse_candidate_never_replaces_the_incumbent():
    P64, V, Q, cons = setup(blade_two_walls)
    good = select_global_frozen(cons, P64, Q, (),
                                area_grid=quad_area_canonical_grid(P64))
    saved = ex.local_improve

    def worse(m, sel, deadline, *a, **k):
        s = sel.copy()
        free = np.nonzero(~s & ~m.protected)[0]
        if len(free):
            s[free[0]] = True                 # feasible but strictly worse
        return {"sel": s, "stopped": "planted", "n_swaps": 0}

    ex.local_improve = worse
    try:
        out = select_global_frozen(cons, P64, Q, (),
                                   area_grid=quad_area_canonical_grid(P64))
    finally:
        ex.local_improve = saved
    assert out["achieved_area"] <= good["achieved_area"] + 1e-12
    assert covered(cons, out["chosen"])


# ------------------------------------------------- (d) the 99.9%-area core
def test_core_is_the_smallest_prefix_reaching_999_of_input_area():
    areas = [1000.0, 500.0, 1.0, 0.4, 0.1]        # total 1501.5
    core = area_core(areas, [1.0] * 5)
    # 1000 + 500 = 1500 is 99.900...% -- 0.99900 of 1501.5 is 1500.0
    assert core["n_core_components"] == 2
    assert core["core_area"] == pytest.approx(1500.0)
    assert core["core_area_fraction"] >= 0.999
    assert core["n_tail_components"] == 3
    assert core["core_gate_pass"] is True


def test_core_closes_the_one_percent_loophole():
    """Fifty genuine 0.9%-of-area components holding 45% of a segment all
    escaped the old ">=1% of segment area" rule. They are IN the core."""
    areas = [55.0] + [0.9] * 50                   # total 100
    rmains = [1.0] + [0.5] * 50                   # every small one shattered
    core = area_core(areas, rmains)
    assert core["n_core_components"] == 51
    assert core["core_gate_pass"] is False
    assert core["n_core_components_below_gate"] == 50
    assert core["min_R_main_core"] == pytest.approx(0.5)


def test_a_genuine_tail_component_does_not_fail_the_gate():
    areas = [1e6, 1.0]                            # the speck is 1e-6 of area
    core = area_core(areas, [0.96, 0.0])
    assert core["n_core_components"] == 1
    assert core["core_gate_pass"] is True
    assert core["min_R_main_core"] == pytest.approx(0.96)
    assert core["n_tail_components"] == 1


def test_core_ordering_is_by_area_not_by_index():
    core = area_core([1.0, 900.0, 100.0], [0.2, 1.0, 1.0])
    assert core["core_component_indices"][:2] == [1, 2]
    assert core["core_gate_pass"] is True         # the 0.2 speck is the tail


def test_single_component_segment_is_entirely_core():
    core = area_core([42.0], [0.5])
    assert core["n_core_components"] == 1
    assert core["core_gate_pass"] is False


# ----------------------------------------------- (e) the policy hash
def test_policy_hash_is_stable_and_sensitive():
    h = frozen_policy_hash()
    assert h == frozen_policy_hash(FROZEN_POLICY)
    assert len(h) == 16
    assert frozen_policy_hash(dict(FROZEN_POLICY,
                                   lp_budget_s_per_segment=61.0)) != h
    assert FROZEN_POLICY["version"] == FROZEN_POLICY_VERSION
    assert FROZEN_POLICY["lp_budget_s_per_segment"] == 60.0
    assert FROZEN_POLICY["improvement_budget_s_per_segment"] == 120.0
    assert FROZEN_POLICY["lp_skip_above_reduced_constraints"] == 50000
    assert FROZEN_POLICY["process_limit_s_per_segment"] == 600.0


@pytest.mark.parametrize("fixture", FIXTURES)
def test_certificate_fields_required_by_round_28(fixture):
    P64, V, Q, cons = setup(fixture)
    out = select_global_frozen(cons, P64, Q, (),
                               area_grid=quad_area_canonical_grid(P64))
    for k in ("policy_version", "policy_hash", "greedy_incumbent_area",
              "greedy_construction_seconds", "selection_status",
              "method_mix", "achieved_area", "combined_lower_bound",
              "combined_lower_bound_complete",
              "minimum_area_claim_admissible",
              "n_components_lp_attempted", "n_components_lp_skipped",
              "n_components_exact_attempted", "n_components_exact_skipped"):
        assert k in out, k
    assert out["policy_hash"] == frozen_policy_hash()


# ------------------------ (f) no minimum-area claim unless all are proven
def test_no_minimum_area_claim_when_any_component_is_unproven():
    P64, V, Q, cons = setup(blade_two_walls)
    out = select_global_frozen(cons, P64, Q, (), policy=ZERO_BUDGETS,
                               area_grid=quad_area_canonical_grid(P64))
    assert out["minimum_area_claim_admissible"] is False
    assert out["selection_status"] != "area_optimal"


def test_lp_is_skipped_with_a_reason_above_the_threshold():
    """The 50,000-constraint skip is a SCHEDULING policy and is recorded as
    a reason on the component, never silently."""
    P64, V, Q, cons = setup(blade_two_walls)
    pol = dict(FROZEN_POLICY, lp_skip_above_reduced_constraints=0)
    out = select_global_frozen(cons, P64, Q, (), policy=pol,
                               area_grid=quad_area_canonical_grid(P64))
    assert out["n_components_lp_attempted"] == 0
    for c in out["components"]:
        assert "SCHEDULING" in c["lp_skipped_reason"]
    assert covered(cons, out["chosen"])


def test_frozen_strategy_is_reachable_from_excise():
    P, V = blade_wall()
    out = ex.excise(P, V, strategy="greedy_first")
    assert out["status"] == "clean"
    cert = out["certificate"]
    assert cert["geometry_status"] == "transverse_clean_certified"
    assert cert["selection"]["strategy"] == "greedy_first"
    assert cert["selection"]["policy_hash"] == frozen_policy_hash()
    assert cert["selection"]["policy_version"] == FROZEN_POLICY_VERSION
    assert cert["selection"]["greedy_incumbent_area"] is not None
    assert cert["output_census"]["clean_both_diagonals"] is True
