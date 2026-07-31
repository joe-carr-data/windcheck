"""Acceptance pins for certified excision (notes/CUTTER-SPEC.md, step 1).

Planted meshes only, in the tests/test_probe_local.py fixture style:
2-row "ribbon" grids whose (y, z) polyline is extruded along x, so quads
are planar walls and crossings are hand-computable 2D segment
intersections; plus one non-planar "fold" quad whose d0 and d1
tessellations cover different surfaces, for the diagonal-specific cases.

Covered acceptance cases: (a) known-optimum single crossing, (b) greedy
non-optimal vs joint MILP optimum, (c) d0-only and d1-only crossings,
(d) crossing set requiring the union of both diagonals' constraints,
(e) no-crossing identity, (f) hole + disconnected component preserved,
(g) area conservation and unique quad ownership, (h) coordinate
bit-identity, (i) shared-support refusal, (j) solver-infeasible path.

Round-24 amendment cases (see CUTTER-SPEC.md section 9): (k) planted
counterexample where the UNIQUE optimum invalidates a non-participant
fourth quad corner (8-corner coverage; participants-only is strictly
worse), (l) two-stage lexicographic solve records, (m) deterministic
mocked "best_found" incumbent emission, (n) deterministic mocked
"limit_no_incumbent" refusal, (o) maxedge-dropped quad near a crossing
(retention semantics match retained_quads exactly).
"""
from __future__ import annotations

from math import fsum

import numpy as np

from windcheck.excise import (SHARED_SUPPORT_LABEL, _solve_milp, census,
                              excise, quad_area_canonical)
from windcheck.intrinsic import retained_quads


def ribbon(yz, invalid_cols=()):
    """2-row prism: P[v,u] = (v, y_u, z_u); quads are planar walls."""
    P = np.array([[[float(v), y, z] for (y, z) in yz] for v in (0.0, 1.0)])
    V = np.ones((2, len(yz)), dtype=bool)
    for c in invalid_cols:
        V[:, c] = False
    return P, V


def blade_wall():
    """(a): blade quad (0,0), area 3, crossed once by wall quad (0,3),
    area 1; a second wall quad (0,4), area 1, does NOT cross. The unique
    minimal-area cut removes quad (0,3); under round-24 quad-level
    coverage the lexicographic stage 2 needs only ONE invalidated
    vertex, and it must be one of the two corners NOT shared with the
    kept neighbour (0,4) -- invalidating (0,4)/(1,4) would also drop
    quad (0,4) (area 2 > 1). Which of (0,3)/(1,3) is chosen is a
    solver-determined tie (spec section 8: not part of any claim)."""
    yz = [(0, 0), (3, 0),            # blade: z=0 wall, y in [0,3]
          (0, 0),                    # invalid separator column
          (0.7, -0.5), (0.7, 0.5), (0.7, 1.5)]   # wall: two stacked quads
    return ribbon(yz, invalid_cols=(2,))


def blade_two_walls():
    """(b): blade quad (0,0), area 3, crossed by wall1 (0,3) and wall2
    (0,6), each area 2. Greedy per-crossing removal takes the cheaper
    wall each time (2 + 2 = 4); the joint optimum removes the blade (3)."""
    yz = [(0, 0), (3, 0),
          (0, 0),
          (0.7, -1.0), (0.7, 1.0),
          (0, 0),
          (2.3, -1.0), (2.3, 1.0)]
    return ribbon(yz, invalid_cols=(2, 5))


def fold_mesh(att1=True, att2=True):
    """Non-planar fold quad (0,0) plus vertical attacker rectangles.

    att1 (quad (0,3), z in [0.25, 0.55]) pierces the fold ONLY under the
    d0 tessellation (fold d0 triangle in plane z=y); att2 (quad (0,6),
    z in [-0.3, 0.15]) pierces ONLY under d1 (plane z=x+y-1). Verified
    by the census keys in the tests below.
    """
    P = np.zeros((2, 8, 3))
    V = np.ones((2, 8), dtype=bool)
    V[:, 2] = False
    V[:, 5] = False
    P[0, 0] = (0, 0, 0)
    P[1, 0] = (1, 0, 0)
    P[0, 1] = (0, 1, 0)
    P[1, 1] = (1, 1, 1)                    # lifted corner: non-planar quad
    if att1:
        P[0, 3] = (0.6, 0.4, 0.25)
        P[1, 3] = (0.8, 0.4, 0.25)
        P[0, 4] = (0.6, 0.4, 0.55)
        P[1, 4] = (0.8, 0.4, 0.55)
    else:
        V[:, 3] = False
        V[:, 4] = False
    if att2:
        P[0, 6] = (0.6, 0.4, -0.3)
        P[1, 6] = (0.8, 0.4, -0.3)
        P[0, 7] = (0.6, 0.4, 0.15)
        P[1, 7] = (0.8, 0.4, 0.15)
    else:
        V[:, 6] = False
        V[:, 7] = False
    return P, V


# ------------------------------------------------------------------ cases
def test_a_known_optimum_single_crossing():
    P, V = blade_wall()
    assert len(census(P, V)) > 0
    r = excise(P, V)
    assert r["status"] == "clean"
    # hand-computed unique minimal-area cut: wall quad (0,3); blade
    # (area 3) and wall quad (0,4) survive. CHANGED in round 24: the old
    # 6-vertex triangle-participant coverage forced TWO invalidated
    # vertices [(0,3),(1,3)]; quad-level 8-corner coverage plus the true
    # lexicographic stage 2 needs exactly ONE (either non-shared corner).
    assert len(r["invalidated_vertices"]) == 1
    assert r["invalidated_vertices"][0] in {(0, 3), (1, 3)}
    assert r["removed_quads"] == [(0, 3)]
    assert r["retained_quads"] == [(0, 0), (0, 4)]
    cert = r["certificate"]
    assert cert["area"]["canonical"]["A_excised"] == 1.0
    assert cert["output_census"]["clean_both_diagonals"]
    assert cert["solver"]["status"] == "optimal"
    assert cert["solver"]["lexicographic"] is True
    # the cut has a boundary: the edge shared with the kept wall quad
    assert cert["excision"]["cut_boundary_edges"] == [[[0, 4], [1, 4]]]
    assert cert["excision"]["cut_boundary_length"] == 1.0
    assert cert["contacts"]["output_submultiset_of_input"]
    assert cert["contacts"]["missing_all_witnessed"]
    chosen = r["invalidated_vertices"][0]
    for rc in cert["removed_contacts"]:
        w = rc["witness"]
        # CHANGED in round 24: scopes are "triangle_participant" vs
        # "quad_retention" (old single "participant" scope was the
        # 6-vertex semantics); the one chosen corner witnesses contacts
        # whose crossing triangle omits it via quad retention.
        assert w is not None
        assert w["witness_scope"] in ("triangle_participant",
                                      "quad_retention")
        assert tuple(w["invalidated_vertex"]) == chosen
        assert w["destroyed_triangle"][1:3] == [0, 3]   # the wall quad
    assert cert["triangle_multisets"]["after_subset_of_before"]


def test_b_greedy_nonoptimal_milp_joint_optimum():
    P, V = blade_two_walls()
    contacts = census(P, V)
    Q = retained_quads(np.asarray(P, np.float64), V)

    def vertex_cost(cnr):
        v, u = cnr
        nv, nu = V.shape
        return fsum(quad_area_canonical(P, v + dv, u + du)
                    for dv in (-1, 0) for du in (-1, 0)
                    if 0 <= v + dv < nv - 1 and 0 <= u + du < nu - 1
                    and Q[v + dv, u + du])

    # greedy per-crossing: cheapest participating vertex, independently
    greedy_verts = {min(c["participants"], key=lambda p: (vertex_cost(p), p))
                    for c in contacts}
    greedy_removed = {(v, u) for v, u in zip(*np.nonzero(Q))
                      if {(a, b) for a in (v, v + 1) for b in (u, u + 1)}
                      & greedy_verts}
    greedy_area = fsum(quad_area_canonical(P, int(v), int(u))
                       for v, u in greedy_removed)
    assert greedy_area == 4.0            # both walls, 2 + 2

    r = excise(P, V)
    assert r["status"] == "clean"
    assert r["removed_quads"] == [(0, 0)]     # the blade: joint optimum
    a_exc = r["certificate"]["area"]["canonical"]["A_excised"]
    assert a_exc == 3.0
    assert a_exc < greedy_area


def test_c_d0_only_and_d1_only_crossings():
    P, V = fold_mesh(att1=True, att2=False)
    keys = [c["key"] for c in census(P, V)]
    assert keys and all(k[0] == 0 for k in keys)     # d0-only constraints
    r = excise(P, V)
    assert r["status"] == "clean"
    assert r["removed_quads"] == [(0, 3)]            # attacker, not fold

    P, V = fold_mesh(att1=False, att2=True)
    keys = [c["key"] for c in census(P, V)]
    assert keys and all(k[0] == 1 for k in keys)     # d1-only constraints
    r = excise(P, V)
    assert r["status"] == "clean"
    assert r["removed_quads"] == [(0, 6)]


def test_d_requires_union_of_both_diagonals():
    P, V = fold_mesh()
    keys = [c["key"] for c in census(P, V)]
    assert {k[0] for k in keys} == {0, 1}
    # resolving only the d0 constraints (cutting attacker 1) is NOT clean:
    # the d1 crossing survives -- both diagonals' constraints are required
    V_d0only = V.copy()
    V_d0only[0, 3] = False
    left = [c["key"] for c in census(P, V_d0only)]
    assert left and all(k[0] == 1 for k in left)

    r = excise(P, V)
    assert r["status"] == "clean"
    assert r["removed_quads"] == [(0, 3), (0, 6)]    # both attackers
    assert r["retained_quads"] == [(0, 0)]           # fold survives
    cert = r["certificate"]
    assert cert["output_census"]["clean_both_diagonals"]
    # non-planar fold: the two tessellations measure different areas,
    # which is exactly why the accounting is per-diagonal
    assert cert["area"]["d0"]["A_input"] != cert["area"]["d1"]["A_input"]


def test_e_no_crossing_identity():
    P, V = ribbon([(0, 0), (1, 0), (2, 0)])
    r = excise(P, V)
    assert r["status"] == "clean"
    assert np.array_equal(r["valid"], V)
    assert r["invalidated_vertices"] == []
    assert r["removed_quads"] == []
    cert = r["certificate"]
    assert cert["excision"]["invalidated_vertices"] == []
    assert cert["solver"]["status"] == "not_required"
    assert cert["solver"]["solves"] == []
    assert cert["output"]["sha256"] == cert["input"]["sha256"]
    assert cert["output_census"]["clean_both_diagonals"]
    assert cert["area"]["canonical"]["clean_recovery_fraction"] == 1.0


def holey_mesh():
    """(f): blade x wall crossing (rows 0-1) plus a DISCONNECTED flat
    patch (cols 9-13, all 4 rows, z=10) with a hole at vertex (1, 11)."""
    P = np.zeros((4, 14, 3))
    V = np.zeros((4, 14), dtype=bool)
    for v in (0, 1):
        P[v, 0] = (v, 0, 0)
        P[v, 1] = (v, 3, 0)                       # blade quad (0,0)
        P[v, 3] = (v, 0.7, -0.5)
        P[v, 4] = (v, 0.7, 0.5)                   # wall quad (0,3)
        V[v, 0] = V[v, 1] = V[v, 3] = V[v, 4] = True
    for v in range(4):
        for u in range(9, 14):
            P[v, u] = (v, u - 9.0, 10.0)          # far flat patch
            V[v, u] = True
    V[1, 11] = False                              # the hole
    return P, V


def test_f_hole_and_disconnected_component_preserved():
    P, V = holey_mesh()
    patch_quads = {(v, u) for v in range(3) for u in range(9, 13)} \
        - {(0, 10), (0, 11), (1, 10), (1, 11)}
    r = excise(P, V)
    assert r["status"] == "clean"
    assert r["removed_quads"] == [(0, 3)]         # the wall, nothing else
    assert patch_quads <= set(r["retained_quads"])
    # untouched pieces preserved exactly: mask identical away from the cut
    assert np.array_equal(r["valid"][:, 5:], V[:, 5:])
    assert np.array_equal(r["valid"][:, :3], V[:, :3])
    # HYBRID invalidation: retained coordinates bit-identical, and the ONLY
    # cells whose coordinates moved are the excised ones (stamped to -1)
    Vo = r["valid"]
    assert r["points"][Vo].tobytes() == P[Vo].tobytes()
    changed = np.any(r["points"] != P, axis=-1)
    assert np.array_equal(changed, V & ~Vo)
    cert = r["certificate"]
    assert cert["topology"]["components_before"] == 3
    assert cert["topology"]["components_after"] == 2
    assert cert["validity"]["changes_only_valid_to_invalid"]


def test_g_area_conservation_and_unique_ownership():
    P, V = blade_two_walls()
    r = excise(P, V)
    cert = r["certificate"]
    for block in (cert["area"]["d0"], cert["area"]["d1"],
                  cert["area"]["canonical"]):
        assert abs(block["A_clean"] + block["A_excised"]
                   + block["A_unresolved"] - block["A_input"]) < 1e-9
        assert block["A_unresolved"] == 0.0
    # unique ownership: removed and kept partition the input-retained set
    removed, kept = set(r["removed_quads"]), set(r["retained_quads"])
    assert not removed & kept
    Q = retained_quads(np.asarray(P, np.float64), V)
    assert removed | kept == {(int(v), int(u))
                              for v, u in zip(*np.nonzero(Q))}
    assert cert["area"]["canonical"]["clean_recovery_fraction"] == 4.0 / 7.0


def test_h_retained_coordinate_bit_identity_float32():
    """Round-25 A1: the guarantee is RETAINED-coordinate bit-identity, and
    invalidation is hybrid (mask cleared AND coordinates stamped -1)."""
    P, V = blade_wall()
    P32 = P.astype(np.float32)                    # real tifxyz dtype
    r = excise(P32, V)
    assert r["status"] == "clean"
    assert r["points"].dtype == np.float32
    Vo = r["valid"]
    assert r["points"][Vo].tobytes() == P32[Vo].tobytes()   # every RETAINED
    excised = V & ~Vo
    assert excised.any()
    assert (r["points"][excised] == -1.0).all()             # stamped
    changed = np.any(r["points"] != P32, axis=-1)
    assert np.array_equal(changed, excised)     # nothing else moved
    assert not (r["valid"] & ~V).any()             # never newly valid
    cert = r["certificate"]
    assert cert["retained_coordinate_bit_identity"]
    assert cert["coordinates_changed_only_at_excised_cells"]
    assert cert["invalidation"]["mask_cleared_at_excised"]
    assert cert["invalidation"]["coordinates_stamped_missing_at_excised"]
    assert cert["invalidation"]["missing_marker"] == -1.0
    # CHANGED in round 24: grid_shape/dtype recorded once at top level
    # (input == output by contract; the duplicate entry was removed)
    assert r["certificate"]["dtype"] == "float32"
    assert r["certificate"]["grid_shape"] == [2, 6]
    assert "dtype" not in r["certificate"]["input"]


def test_i_shared_support_refused():
    # spiral ribbon: s0 x s2, s1 x s3, s0 x s3 -- one event whose two
    # branch regions are grid-adjacent (share support at column 2)
    P, V = ribbon([(0, 0), (2, 0), (1.5, 1), (0.5, -1), (2.5, 0.8)])
    assert len(census(P, V)) > 0
    r = excise(P, V)
    assert r["status"] == "refused_shared_support"
    assert r["label"] == SHARED_SUPPORT_LABEL
    assert "points" not in r and "valid" not in r
    assert "certificate" not in r
    assert r["events"]                            # the offending regions


def test_j_infeasible_reported_without_clean_claim():
    P, V = blade_wall()
    # protect the full 8-corner coverage set: no constraint can be met
    protected = [(v, u) for v in (0, 1) for u in (0, 1, 3, 4)]
    r = excise(P, V, protected=protected)
    assert r["status"] == "infeasible"
    assert "points" not in r and "valid" not in r
    assert "certificate" not in r
    solves = r["solver"]["solves"]
    assert len(solves) == 1 and solves[0]["raw_status"] == 2
    assert solves[0]["stage"] == 1                # stage 2 never attempted


# ------------------------------------------- round-24 amendment cases
def fourth_corner_mesh():
    """(k): planted counterexample forcing a NON-PARTICIPANT fourth-corner
    optimum. A non-planar fold quad (1,1) (corners (0,0,0), (1,0,0),
    (0,1,0), lifted (1,1,1)) is pierced by a tall attacker wall, quad
    (1,5), ONLY through the fold's d0 k=1 triangle (the plane z=y): the
    attacker lives in plane y=0.4, x in [0.6,0.8], z in [0.25,20], which
    misses the fold's d0 k=0 triangle (slice x in [0,0.4]) and both d1
    triangles (z=0 and z=x+y-1 slices). Triangle participants are
    therefore fold corners {(1,1),(2,1),(2,2)} plus all four attacker
    corners; fold corner (1,2) participates in NO crossing triangle.

    Costs (canonical areas): fold ~1.3901; big planar quad (1,0)
    (area 5) hangs off fold corners (1,1),(2,1); big quad (2,1)
    (area ~7.071) hangs off (2,1),(2,2); the attacker is a 0.2 x 19.75
    wall (area 3.95), so EVERY participant vertex costs > 3.9 while the
    fourth corner (1,2) supports only the fold itself. The unique
    optimum invalidates (1,2) -- reachable only through quad-level
    8-corner coverage."""
    P = np.zeros((4, 7, 3))
    V = np.zeros((4, 7), dtype=bool)
    P[1, 1] = (0, 0, 0)
    P[2, 1] = (1, 0, 0)
    P[1, 2] = (0, 1, 0)
    P[2, 2] = (1, 1, 1)                           # fold quad (1,1)
    P[1, 0] = (0, -5, 0)
    P[2, 0] = (1, -5, 0)                          # big quad (1,0), area 5
    P[3, 1] = (6, 0, 0)
    P[3, 2] = (6, 1, 1)                           # big quad (2,1), ~7.071
    P[1, 5] = (0.6, 0.4, 0.25)
    P[2, 5] = (0.8, 0.4, 0.25)
    P[1, 6] = (0.6, 0.4, 20.0)
    P[2, 6] = (0.8, 0.4, 20.0)                    # attacker quad (1,5)
    for (v, u) in [(1, 0), (2, 0), (1, 1), (2, 1), (1, 2), (2, 2),
                   (3, 1), (3, 2), (1, 5), (2, 5), (1, 6), (2, 6)]:
        V[v, u] = True
    return P, V


def test_k_fourth_corner_quad_retention_optimum():
    P, V = fourth_corner_mesh()
    P64 = np.asarray(P, np.float64)
    contacts = census(P, V)
    assert contacts and all(c["key"][0] == 0 for c in contacts)
    assert all(c["key"][5] == 1 for c in contacts)   # fold k1 only
    participants = {p for c in contacts for p in c["participants"]}
    coverage = {p for c in contacts for p in c["coverage"]}
    assert (1, 2) in coverage and (1, 2) not in participants

    r = excise(P, V)
    assert r["status"] == "clean"
    assert r["invalidated_vertices"] == [(1, 2)]     # the fourth corner
    assert r["removed_quads"] == [(1, 1)]            # the fold only
    fold_area = quad_area_canonical(P64, 1, 1)
    cert = r["certificate"]
    assert abs(cert["area"]["canonical"]["A_excised"] - fold_area) < 1e-12
    assert cert["output_census"]["clean_both_diagonals"]
    assert cert["solver"]["status"] == "optimal"
    assert cert["solver"]["lexicographic"] is True
    for rc in cert["removed_contacts"]:
        w = rc["witness"]
        assert w["witness_scope"] == "quad_retention"
        assert w["invalidated_vertex"] == [1, 2]
        assert w["destroyed_triangle"] == [0, 1, 1, 1]

    # restricting the MILP to triangle participants is STRICTLY worse:
    # solve a second MILP whose coverage sets are the participants only
    Q = retained_quads(P64, V)
    restricted = [dict(c, coverage=c["participants"]) for c in contacts]
    full = _solve_milp(contacts, P64, Q, set(), None)
    part = _solve_milp(restricted, P64, Q, set(), None)
    assert full["records"][0]["status"] == "optimal"
    assert part["records"][0]["status"] == "optimal"
    assert full["chosen"] == [(1, 2)]
    assert all(v in participants for v in part["chosen"])
    # 1.3901 (fold via fourth corner) vs 3.95 (cheapest participant cut)
    assert full["records"][0]["objective"] \
        < part["records"][0]["objective"] - 2.5


def test_l_two_stage_lexicographic_records():
    P, V = blade_wall()
    r = excise(P, V)
    assert r["status"] == "clean"
    solver = r["certificate"]["solver"]
    solves = solver["solves"]
    assert [s["stage"] for s in solves] == [1, 2]
    s1, s2 = solves
    assert s1["purpose"] == "min_excised_area"
    assert s1["status"] == "optimal"
    assert abs(s1["objective"] - 1.0) < 1e-9      # pure area, no epsilon
    assert s1["dual_bound"] is not None and s1["mip_gap"] is not None
    assert s2["purpose"] == "min_invalidated_vertices"
    assert s2["status"] == "optimal"
    assert abs(s2["objective"] - 1.0) < 1e-9      # one vertex suffices
    assert s2["area_allowance"] > 0.0             # declared and recorded
    assert abs(s2["area_cap"] - (s1["objective"] + s2["area_allowance"])) \
        < 1e-15
    assert solver["lexicographic"] is True
    assert "stage2_allowance_rule" in solver
    assert solver["final_stage1_objective"] == s1["objective"]
    assert solver["optimality_gap"] == s1["mip_gap"]


def test_m_mocked_best_found_incumbent(monkeypatch):
    """Deterministic 'best_found': the real solve runs, then the result is
    doctored to HiGHS status 1 (time limit) with a declared gap -- no
    millisecond-timing dependence. The incumbent is the true optimum, so
    the recensus is clean and emission is allowed; the certificate must
    say best_found, record achieved area + gap + bound, skip stage 2 and
    never claim minimum-area or lexicographic optimality."""
    import windcheck.excise as ex
    real = ex._run_milp

    def doctored(c, A, lb, ub, integrality, lower, upper, options):
        res = real(c, A, lb, ub, integrality, lower, upper, options)
        res.status = 1
        res.mip_gap = 0.25
        res.mip_dual_bound = 0.75 * float(res.fun)
        return res

    monkeypatch.setattr(ex, "_run_milp", doctored)
    P, V = blade_wall()
    r = ex.excise(P, V)
    assert r["status"] == "clean"          # clean comes from the recensus
    cert = r["certificate"]
    solver = cert["solver"]
    assert solver["status"] == "best_found"
    assert solver["lexicographic"] is False
    solves = solver["solves"]
    assert [s["stage"] for s in solves] == [1]      # stage 2 skipped
    assert solves[0]["status"] == "best_found"
    assert solves[0]["objective"] is not None       # ACHIEVED area
    assert solves[0]["mip_gap"] == 0.25             # gap recorded
    assert solves[0]["dual_bound"] is not None      # bound recorded
    assert "NOT minimum-area" in solves[0]["note"]
    assert "not a proven minimum" in solver["note"]
    assert cert["output_census"]["clean_both_diagonals"]
    assert solver["optimality_gap"] == 0.25


def test_n_mocked_limit_no_incumbent(monkeypatch):
    """Deterministic 'limit_no_incumbent': the mocked solver returns
    HiGHS status 1 with no incumbent; excision must refuse -- no output
    arrays, no certificate, no clean claim."""
    import windcheck.excise as ex

    class _Res:
        status = 1
        x = None
        fun = None
        message = "time limit reached without incumbent (mocked)"

    monkeypatch.setattr(ex, "_run_milp", lambda *a, **k: _Res())
    P, V = blade_wall()
    r = ex.excise(P, V)
    assert r["status"] == "limit_no_incumbent"
    assert "points" not in r and "valid" not in r
    assert "certificate" not in r
    solves = r["solver"]["solves"]
    assert len(solves) == 1
    assert solves[0]["status"] == "limit_no_incumbent"
    assert solves[0]["raw_status"] == 1
    assert solves[0]["stage"] == 1                  # stage 2 never ran


def stretched_neighbour_mesh():
    """(o): blade x wall crossing as in blade_wall, but the wall's u+1
    neighbour quad (0,4) has side edges of length 100 > maxedge (60):
    all four of its corners are VALID yet selfcross drops it by edge
    length. Its corners (0,4),(1,4) sit inside the crossing pair's
    8-corner coverage set, so if retention semantics diverged from
    retained_quads the MILP would price (or remove) a phantom quad."""
    yz = [(0, 0), (3, 0), (0, 0),
          (0.7, -0.5), (0.7, 0.5), (0.7, 100.5)]
    return ribbon(yz, invalid_cols=(2,))


def test_o_maxedge_dropped_quad_near_crossing():
    P, V = stretched_neighbour_mesh()
    P64 = np.asarray(P, np.float64)
    Q = retained_quads(P64, V)                    # default maxedge 60
    assert V[0, 4] and V[1, 4] and V[0, 5] and V[1, 5]
    assert not Q[0, 4]                            # dropped by edge length
    assert Q[0, 0] and Q[0, 3]
    contacts = census(P, V)
    assert contacts
    assert all(q in {(0, 0), (0, 3)}
               for c in contacts for q in (c["q1"], c["q2"]))

    r = excise(P, V)
    assert r["status"] == "clean"
    assert r["removed_quads"] == [(0, 3)]
    # retention semantics match retained_quads EXACTLY: the maxedge-
    # dropped quad is neither removed nor retained -- it was never part
    # of the censused complex, so it is never priced and never counted
    owned = set(r["removed_quads"]) | set(r["retained_quads"])
    assert owned == {(int(v), int(u)) for v, u in zip(*np.nonzero(Q))}
    assert (0, 4) not in owned
    cert = r["certificate"]
    assert abs(cert["area"]["canonical"]["A_input"] - 4.0) < 1e-12
    assert abs(cert["area"]["canonical"]["A_excised"] - 1.0) < 1e-12
    # one vertex suffices; corners of the dropped quad cost the same as
    # the wall's own free corners (its area is NOT added), so any of the
    # four wall corners is an equal-area, equal-count solver tie
    assert len(r["invalidated_vertices"]) == 1
    assert r["invalidated_vertices"][0] in {(0, 3), (1, 3), (0, 4), (1, 4)}


# ------------------------- round-25/26 amendments (CUTTER-SPEC section 10)
def naive_valid(P: np.ndarray) -> np.ndarray:
    """A NAIVE consumer: implements ONLY tifxyz's x=y=z=-1 convention and
    ignores the validity mask / mask.tif sidecar entirely."""
    return ~np.all(np.asarray(P) == -1.0, axis=-1)


def as_tifxyz(P, V):
    """Mark missing cells the way real tifxyz files do, so a naive reader
    and windcheck's reader agree on the INPUT."""
    P = np.asarray(P, np.float64).copy()
    P[~V] = -1.0
    return P, V


def test_p_naive_reader_ignoring_mask_sees_no_crossing():
    """Round-25 A1. A consumer that never looks at the mask must still see
    the crossing gone -- that is the whole point of the hybrid stamps."""
    P, V = as_tifxyz(*blade_two_walls())
    # CONTROL: mask-only invalidation is exactly "hand the naive consumer
    # the input back", and the input crosses itself.
    assert np.array_equal(naive_valid(P), V)
    assert census(P, naive_valid(P))

    r = excise(P, V)
    assert r["status"] == "clean"
    V_naive = naive_valid(r["points"])            # mask.tif IGNORED
    assert not np.array_equal(V_naive, V)         # the stamps are visible
    assert np.array_equal(V_naive, r["valid"])    # both conventions agree
    assert census(r["points"], V_naive) == []     # ...and it is clean


def test_p2_naive_reader_on_a_hole_and_disconnected_mesh():
    P, V = as_tifxyz(*holey_mesh())
    r = excise(P, V)
    assert r["status"] == "clean"
    V_naive = naive_valid(r["points"])
    assert np.array_equal(V_naive, r["valid"])
    assert census(r["points"], V_naive) == []
    # the far patch and its hole are untouched under the naive reading too
    assert np.array_equal(V_naive[:, 9:], V[:, 9:])


def test_q_reduction_never_changes_the_optimal_area():
    """Round-26 Q2: dedup + dominance + decomposition are objective-
    preserving. Solve reduced and unreduced; the excised area must match."""
    for name, (P, V) in (("blade_wall", blade_wall()),
                         ("blade_two_walls", blade_two_walls()),
                         ("holey", holey_mesh()),
                         ("fold", fold_mesh()),
                         ("stretched", stretched_neighbour_mesh())):
        red = excise(P, V, reduce=True)
        raw = excise(P, V, reduce=False)
        assert red["status"] == raw["status"] == "clean", name
        a_red = red["certificate"]["area"]["canonical"]["A_excised"]
        a_raw = raw["certificate"]["area"]["canonical"]["A_excised"]
        assert abs(a_red - a_raw) < 1e-9, (name, a_red, a_raw)
        assert len(red["invalidated_vertices"]) \
            == len(raw["invalidated_vertices"]), name


def test_q2_reduction_statistics_recorded():
    P, V = blade_two_walls()
    r = excise(P, V)
    red = r["certificate"]["solver"]["constraint_reduction"]
    assert red["n_raw"] >= red["n_after_dedup"] >= red["n_after_dominance"]
    assert 0 < red["total_ratio"] <= 1.0
    assert red["n_components"] >= 1
    assert sum(red["component_sizes"]) == red["n_after_dominance"]


def test_r_dedup_and_dominance_are_exact():
    from windcheck.excise import reduce_constraints
    Q = np.ones((6, 6), bool)
    def con(key, cov):
        return {"key": key, "coverage": sorted(cov)}
    cons = [
        con("a", [(0, 0), (0, 1), (1, 0), (1, 1)]),
        con("b", [(0, 0), (0, 1), (1, 0), (1, 1)]),      # duplicate of a
        con("c", [(0, 0), (0, 1), (1, 0), (1, 1), (2, 2)]),  # superset of a
        con("d", [(3, 3), (3, 4)]),                      # independent
    ]
    red = reduce_constraints(cons, Q)
    assert red["n_raw"] == 4
    assert red["n_after_dedup"] == 3                  # b collapsed into a
    assert red["n_after_dominance"] == 2              # c dominated by a
    kept = {c["key"] for c in red["kept"]}
    assert kept == {"a", "d"}
    assert red["n_components"] == 2                   # a and d are disjoint
    assert red["dedup_ratio"] == 3 / 4
    assert red["total_ratio"] == 2 / 4


def test_r2_disconnected_components_are_solved_independently():
    """Two crossings far apart: the reducer must split them, and the split
    solve must reach the same optimum as one joint MILP."""
    P, V = holey_mesh()
    # plant a second, independent crossing inside the far patch
    P[3, 12] = (3, 3.0, 10.0)
    P[2, 12] = (2, 3.0, 10.0)
    P[0, 13] = (0, 4.0, 10.0)
    r = excise(P, V, reduce=True)
    raw = excise(P, V, reduce=False)
    assert r["status"] == raw["status"] == "clean"
    assert (r["certificate"]["area"]["canonical"]["A_excised"]
            == raw["certificate"]["area"]["canonical"]["A_excised"])


def test_s_junction_excision_when_shared_support_is_not_refused():
    """Round-26 Q1: shared support is mathematically excisable. Refused by
    default; cut and LABELLED junction_excision when asked."""
    from windcheck.excise import JUNCTION_EXCISION_LABEL
    P, V = ribbon([(0, 0), (2, 0), (1.5, 1), (0.5, -1), (2.5, 0.8)])
    assert excise(P, V)["status"] == "refused_shared_support"
    r = excise(P, V, refuse_shared_support=False)
    assert r["status"] == "clean"
    assert r["junction_excision"] is True
    cert = r["certificate"]
    assert cert["junction_excision"]["label"] == JUNCTION_EXCISION_LABEL
    assert cert["junction_excision"]["events"]
    assert cert["output_census"]["clean_both_diagonals"]
    assert census(r["points"], r["valid"]) == []


def test_t_segment_wide_constraints_cover_every_transverse_row():
    """The MILP is built from EVERY censused pair under BOTH diagonals in
    ONE solve -- not one event at a time."""
    P, V = fold_mesh()                     # d0-only and d1-only attackers
    contacts = census(P, V)
    diags = {c["key"][0] for c in contacts}
    assert diags == {0, 1}
    r = excise(P, V)
    assert r["status"] == "clean"
    cert = r["certificate"]
    assert cert["solver"]["constraint_reduction"]["n_raw"] == len(contacts)
    # one iteration: emit once, recensus once
    assert cert["solver"]["iterations"] == 1
    assert cert["output_census"]["clean_both_diagonals"]
