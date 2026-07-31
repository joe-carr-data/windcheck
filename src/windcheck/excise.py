"""Certified excision: census-clean by removing surface, never by moving it.

Implements the pre-registered contract in notes/CUTTER-SPEC.md (step 1:
planted meshes only). ~40% of crossing events are certified
rigid-infeasible: no displacement within budget can separate them. For
those the honest operation is a cut — a topology change, labelled as such:

- The output is one AGGREGATE grid with the input's dimensions; every
  RETAINED coordinate is bit-identical to the input; validity changes only
  valid -> invalid, and invalidation is HYBRID (round-25 A1): the validity
  mask is cleared AND the coordinates are stamped x=y=z=-1, so a consumer
  that honours only one of tifxyz's two missing-cell conventions still
  sees the cell as absent. Per-piece cleanliness is rejected as gameable
  (splitting colliding geometry into clean files converts
  self-intersection into inter-piece intersection).
- The cut is chosen by a vertex-mask MILP (scipy/HiGHS): binary x_v
  (vertex invalidated) and y_q (quad removed). tifxyz validity is
  QUAD-level: invalidating any corner of a quad drops the whole quad and
  both of its triangles, so the coverage set of a crossing triangle pair
  is all four corners of quad1 UNION all four corners of quad2 (normally
  8 vertices) -- a non-participant fourth corner can be the cheaper cut.
  y_q >= x_v for every quad incident to v. The objective is truly
  lexicographic in two stages: stage 1 minimizes excised area ONLY; iff
  stage 1 is proven optimal, stage 2 minimizes the invalidated-vertex
  count subject to area <= stage-1 optimum + a declared, recorded
  allowance (STAGE2_ALLOWANCE_RULE). A time-limited stage-1 incumbent
  skips stage 2 and is never called minimum-area or lexicographically
  optimized.
- Shared-support (self-touching) events are a THIRD class the cutter
  REFUSES (SHARED_SUPPORT_LABEL): a vertex mask there destroys both
  branches at once and cannot be certified as a branch separation.
- The clean claim only ever comes from the full both-diagonal recensus of
  the emitted arrays. Solver "optimal" vs "best_found" is recorded;
  infeasible or limit-without-incumbent emits nothing.

The census here is the planted-size stand-in for engines/selfcross: all
retained quad pairs with quad-origin Chebyshev distance > exclude (=1),
each 2x2 triangle combination under both diagonals, transversality by
the float64 interval predicate `_tri_tri_segment`. Triangle ordering
matches SurfaceGraph.quad_triangles / selfcross.cpp exactly.
"""
from __future__ import annotations

import hashlib
import subprocess
import time
from collections import Counter
from math import fsum
from pathlib import Path

import numpy as np

from .check import PAIR_DTYPE
from .intrinsic import (MAXEDGE_DEFAULT, _tri_tri_segment, oriented_events,
                        retained_quads)

EXCLUDE_DEFAULT = 1                     # selfcross default adjacency exclusion

MISSING = -1.0                          # tifxyz's coordinate missing marker

SHARED_SUPPORT_LABEL = (
    "shared-support topology; displacement and ordinary branch-separation "
    "cutter inapplicable; remesh/junction-excision required")

JUNCTION_EXCISION_LABEL = (
    "junction_excision: the cut passes through shared support (a "
    "self-touching junction). Mathematically excisable and certified the "
    "same way as any other excision, but it is NEVER a branch separation "
    "and the retained side is NOT claimed to be a distinct sheet")

# Round-25 A1: mask-only invalidation is an interoperability defect -- a
# consumer that ignores the sidecar reconstructs the original crossing. The
# artifact contract is therefore HYBRID, and the coordinate guarantee is
# stated over RETAINED cells only.
HYBRID_INVALIDATION = (
    "HYBRID invalidation (round-25 A1): every excised cell is marked BOTH "
    "ways -- validity/mask = 0 AND x = y = z = -1. tifxyz admits two "
    "missing-cell conventions and a consumer implementing only one of them "
    "must still see the cell as absent, otherwise it reconstructs the "
    "crossing this operation removed.")

RETAINED_BIT_IDENTITY = (
    "Every RETAINED coordinate is bit-identical to the input. Coordinates "
    "differ from the input at excised cells ONLY, where they are the -1 "
    "missing marker; no retained point is moved, smoothed or re-quantized.")

STALENESS_WARNING = (
    "Derived artifacts computed from the input mesh (censuses, spectra, "
    "atlases, clearance and geodesic products) are STALE for this output; "
    "recompute them against the output hash before use.")


STAGE2_ALLOWANCE_RULE = (
    "stage-2 area cap = stage-1 proven optimum + allowance; allowance = "
    "1e-6 + 1e-9 * |optimum| (the 1e-6 term absorbs HiGHS's default "
    "absolute MIP gap, for which scipy exposes no documented option; "
    "mip_rel_gap is set to 0.0)")


# ------------------------------------------------------------------ census
def quad_corners(v: int, u: int):
    """The four grid corners of quad (v, u)."""
    return ((v, u), (v + 1, u), (v, u + 1), (v + 1, u + 1))


def quad_triangle_corners(v: int, u: int, diag: int):
    """Grid corners of a quad's two triangles, selfcross's order."""
    if diag == 0:
        return (((v, u), (v, u + 1), (v + 1, u + 1)),
                ((v, u), (v + 1, u + 1), (v + 1, u)))
    return (((v, u), (v, u + 1), (v + 1, u)),
            ((v, u + 1), (v + 1, u + 1), (v + 1, u)))


def census(P: np.ndarray, V: np.ndarray, maxedge: float = MAXEDGE_DEFAULT,
           exclude: int = EXCLUDE_DEFAULT, diagonals=(0, 1)) -> list[dict]:
    """Global transverse crossing census, pure Python, planted sizes.

    Contact key = (diag, v1, u1, v2, u2, k1, k2) with (v1,u1) < (v2,u2);
    identity is stable across masking because coordinates never change.
    """
    P = np.asarray(P, np.float64)
    Q = retained_quads(P, V, maxedge)
    qs = [(int(v), int(u)) for v, u in zip(*np.nonzero(Q))]
    lo, hi = {}, {}
    for v, u in qs:
        C = np.stack([P[v, u], P[v + 1, u], P[v, u + 1], P[v + 1, u + 1]])
        lo[(v, u)], hi[(v, u)] = C.min(axis=0), C.max(axis=0)
    out = []
    for i, (v1, u1) in enumerate(qs):
        for (v2, u2) in qs[i + 1:]:
            if max(abs(v1 - v2), abs(u1 - u2)) <= exclude:
                continue                # census adjacency exclusion
            if (hi[(v1, u1)] < lo[(v2, u2)]).any() \
                    or (hi[(v2, u2)] < lo[(v1, u1)]).any():
                continue                # AABB prefilter
            for diag in diagonals:
                for k1, c1 in enumerate(quad_triangle_corners(v1, u1, diag)):
                    T1 = np.array([P[a] for a in c1])
                    for k2, c2 in enumerate(
                            quad_triangle_corners(v2, u2, diag)):
                        seg = _tri_tri_segment(T1, np.array([P[a] for a in c2]))
                        if seg is not None:
                            out.append({
                                "key": (diag, v1, u1, v2, u2, k1, k2),
                                "q1": (v1, u1), "q2": (v2, u2),
                                "corners1": c1, "corners2": c2,
                                "participants": sorted(set(c1) | set(c2)),
                                # quad-level validity: invalidating ANY
                                # corner of either quad drops that quad
                                # and destroys the crossing triangle
                                "coverage": sorted(
                                    set(quad_corners(v1, u1))
                                    | set(quad_corners(v2, u2))),
                                "seg_len": float(
                                    np.linalg.norm(seg[1] - seg[0])),
                            })
    out.sort(key=lambda c: c["key"])
    return out


def quad_area(P: np.ndarray, v: int, u: int, diag: int) -> float:
    """Area of a quad under diagonal `diag`: sum of its two triangle areas."""
    a = 0.0
    for c in quad_triangle_corners(v, u, diag):
        p0, p1, p2 = (np.asarray(P[i], np.float64) for i in c)
        a += 0.5 * float(np.linalg.norm(np.cross(p1 - p0, p2 - p0)))
    return a


def quad_area_canonical(P: np.ndarray, v: int, u: int) -> float:
    """Declared canonical per-quad area: mean of the d0 and d1 areas."""
    return 0.5 * (quad_area(P, v, u, 0) + quad_area(P, v, u, 1))


# -------------------------------------------------------------------- MILP
def _run_milp(c, A, lb, ub, integrality, lower, upper, options):
    """The single scipy/HiGHS entry point (wrapped so tests can mock
    deterministic solver states without millisecond timing games)."""
    from scipy.optimize import Bounds, LinearConstraint, milp
    return milp(c=c, constraints=LinearConstraint(A, lb, ub),
                integrality=integrality, bounds=Bounds(lower, upper),
                options=options)


def _milp_status(res) -> str:
    if res.status == 0:
        return "optimal"
    if res.status == 1 and getattr(res, "x", None) is not None:
        return "best_found"
    if res.status == 1:
        return "limit_no_incumbent"
    if res.status == 2:
        return "infeasible"
    return "solver_error"


def _solve_milp(constraints: list[dict], P64: np.ndarray, Q_in: np.ndarray,
                protected: set, time_limit: float | None,
                stage2: bool = True) -> dict:
    """One lexicographic (two-stage) vertex-mask solve.

    Stage 1 minimizes excised area ONLY. Iff stage 1 is proven optimal,
    stage 2 adds `area <= optimum + allowance` (STAGE2_ALLOWANCE_RULE,
    recorded per solve) and minimizes the invalidated-vertex count.
    A time-limited stage-1 incumbent skips stage 2: it is an incumbent,
    never "minimum-area", never "lexicographically optimized"; a clean
    claim for it can only come from the full recensus of emitted arrays.
    Tertiary criteria (boundary length, component count) remain a
    documented v2 extension.
    """
    import scipy
    from scipy import sparse

    nv, nu = Q_in.shape[0] + 1, Q_in.shape[1] + 1
    verts = sorted({p for c in constraints for p in c["coverage"]})
    vpos = {v: i for i, v in enumerate(verts)}
    quads = sorted({(v + dv, u + du) for v, u in verts
                    for dv in (-1, 0) for du in (-1, 0)
                    if 0 <= v + dv < nv - 1 and 0 <= u + du < nu - 1
                    and Q_in[v + dv, u + du]})
    qpos = {q: len(verts) + i for i, q in enumerate(quads)}
    n = len(verts) + len(quads)
    areas = [quad_area_canonical(P64, *q) for q in quads]

    rows, cols, data, lb, ub = [], [], [], [], []
    r = 0
    for con in constraints:            # coverage: 8 quad corners, sum >= 1
        for p in con["coverage"]:
            rows.append(r); cols.append(vpos[p]); data.append(1.0)
        lb.append(1.0); ub.append(np.inf)
        r += 1
    for q in quads:                    # coupling: y_q - x_v >= 0
        for cnr in quad_corners(*q):
            if cnr in vpos:
                rows.append(r); cols.append(qpos[q]); data.append(1.0)
                rows.append(r); cols.append(vpos[cnr]); data.append(-1.0)
                lb.append(0.0); ub.append(np.inf)
                r += 1
    lower = np.zeros(n)
    upper = np.ones(n)
    for i, v in enumerate(verts):
        if v in protected:
            upper[i] = 0.0             # x_v forced 0: protected vertex
    integrality = np.ones(n)
    options: dict = {"mip_rel_gap": 0.0}    # documented scipy option
    if time_limit:
        options["time_limit"] = float(time_limit)

    def solve(c_vec, A, lbv, ubv, stage, purpose):
        t0 = time.perf_counter()
        res = _run_milp(np.asarray(c_vec, float), A,
                        np.asarray(lbv, float), np.asarray(ubv, float),
                        integrality, lower, upper, options)
        dt = time.perf_counter() - t0
        fun = getattr(res, "fun", None)
        gap = getattr(res, "mip_gap", None)
        bound = getattr(res, "mip_dual_bound", None)
        rec = {"stage": stage, "purpose": purpose,
               "status": _milp_status(res),
               "raw_status": int(res.status),
               "message": str(getattr(res, "message", "")),
               "objective": (float(fun) if fun is not None else None),
               "dual_bound": (float(bound) if bound is not None else None),
               "mip_gap": (float(gap) if gap is not None else None),
               "solve_time_s": dt, "n_vertex_vars": len(verts),
               "n_quad_vars": len(quads), "n_constraints": int(A.shape[0])}
        return res, rec

    def picks(res):
        return [verts[i] for i in range(len(verts)) if res.x[i] > 0.5]

    # ---- stage 1: minimize excised area ONLY ---------------------------
    A1 = sparse.coo_matrix((data, (rows, cols)), shape=(r, n))
    c1 = np.concatenate([np.zeros(len(verts)), np.asarray(areas, float)])
    res1, rec1 = solve(c1, A1, lb, ub, 1, "min_excised_area")
    records = [rec1]
    chosen: list = []
    lexicographic = False

    if rec1["status"] == "optimal" and not stage2:
        # Round-27 policy: stage 2 is a vertex-count refinement AFTER area is
        # already optimal, and it is not prize-relevant on big components
        # (655 s was measured for one). Skipping it leaves the stage-1
        # optimum, which is area-optimal but not lexicographically optimized.
        rec1["note"] = ("stage 2 skipped by policy (component above the "
                        "stage-2 size threshold): area-optimal, "
                        "vertex count NOT minimized")
        chosen = picks(res1)
    elif rec1["status"] == "optimal":
        # ---- stage 2: minimize invalidated vertices at the optimum ----
        a_star = rec1["objective"]
        allowance = 1e-6 + 1e-9 * abs(a_star)
        A2 = sparse.coo_matrix(
            (data + list(areas),
             (rows + [r] * len(quads), cols + [qpos[q] for q in quads])),
            shape=(r + 1, n))
        c2 = np.concatenate([np.ones(len(verts)), np.zeros(len(quads))])
        res2, rec2 = solve(c2, A2, lb + [-np.inf], ub + [a_star + allowance],
                           2, "min_invalidated_vertices")
        rec2["area_cap"] = a_star + allowance
        rec2["area_allowance"] = allowance
        records.append(rec2)
        if rec2["status"] in ("optimal", "best_found") \
                and getattr(res2, "x", None) is not None:
            chosen = picks(res2)
            lexicographic = rec2["status"] == "optimal"
        else:                          # defensive: stage-1 x is feasible
            rec2["note"] = ("stage 2 returned no usable incumbent; "
                            "falling back to the stage-1 optimum "
                            "(area-optimal, vertex count not minimized)")
            chosen = picks(res1)
    elif rec1["status"] == "best_found":
        rec1["note"] = (
            "time-limited incumbent: NOT minimum-area and NOT "
            "lexicographically optimized; stage 2 skipped; the achieved "
            "area, dual bound and MIP gap are recorded in this record; "
            "any clean claim comes only from the full recensus of the "
            "emitted arrays")
        chosen = picks(res1)

    return {"status": rec1["status"], "chosen": chosen,
            "records": records, "lexicographic": lexicographic,
            "scipy_version": scipy.__version__}


# ------------------------------------------------- constraint reduction (Q2)
REDUCTION_RULE = (
    "Round-26 Q2 reduction, applied to the segment-wide constraint set "
    "before any solve and provably objective-preserving: (a) DEDUPLICATION "
    "-- constraints with identical coverage sets are one constraint; "
    "(b) DOMINANCE -- if coverage(A) is a subset of coverage(B) then any "
    "mask hitting A hits B, so B is redundant and is dropped; "
    "(c) DECOMPOSITION -- constraints are partitioned by connected "
    "components of the shared-variable graph (two constraints are linked "
    "when their incident-quad scopes intersect). Distinct components share "
    "no vertex or quad variable, so the component MILPs are independent and "
    "their optima sum to the global optimum.")


def incident_quads(v: int, u: int, Q_in: np.ndarray):
    """Retained quads having (v, u) as a corner."""
    nv, nu = Q_in.shape
    return [(v + dv, u + du) for dv in (-1, 0) for du in (-1, 0)
            if 0 <= v + dv < nv and 0 <= u + du < nu and Q_in[v + dv, u + du]]


def _dedup(constraints: list[dict]) -> tuple[list[dict], dict]:
    """Collapse constraints with identical coverage sets."""
    groups: dict = {}
    for c in constraints:
        groups.setdefault(frozenset(map(tuple, c["coverage"])), []).append(c)
    kept = []
    for cov in sorted(groups, key=lambda s: (len(s), sorted(s))):
        members = groups[cov]
        rep = dict(members[0])
        rep["multiplicity"] = len(members)
        kept.append(rep)
    return kept, {"n_in": len(constraints), "n_out": len(kept)}


def _dominance(constraints: list[dict]) -> tuple[list[dict], dict]:
    """Drop every coverage set that strictly contains another one.

    Processed in increasing set size, so any subset of the candidate has
    already been decided; if that subset was itself dropped, a still smaller
    kept subset of it also lies inside the candidate. Containment is tested
    by counting, per already-kept set A, how many of A's vertices lie in the
    candidate S: A is a subset of S exactly when that count reaches |A|. Only
    kept sets sharing a vertex with S are ever touched.
    """
    order = sorted(constraints,
                   key=lambda c: (len(c["coverage"]), sorted(c["coverage"])))
    kept: list[dict] = []
    sizes: list[int] = []
    posting: dict = {}
    for c in order:
        S = frozenset(map(tuple, c["coverage"]))
        hits: Counter = Counter()
        for p in S:
            hits.update(posting.get(p, ()))
        if any(n >= sizes[i] for i, n in hits.items()):
            continue                    # dominated: hitting the subset hits S
        idx = len(kept)
        kept.append(c)
        sizes.append(len(S))
        for p in S:
            posting.setdefault(p, []).append(idx)
    return kept, {"n_in": len(constraints), "n_out": len(kept)}


def _components(constraints: list[dict], Q_in: np.ndarray
                ) -> list[list[dict]]:
    """Partition constraints into independent subproblems.

    A constraint's variable scope is its coverage vertices plus every
    retained quad incident to them. Union-find over those quads: two
    constraints land in the same class exactly when their scopes touch, so
    different classes share no decision variable at all.
    """
    parent: dict = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    scopes = []
    for c in constraints:
        scope = sorted({q for p in c["coverage"]
                        for q in incident_quads(*p, Q_in)})
        scopes.append(scope)
        for q in scope[1:]:
            union(scope[0], q)
    groups: dict = {}
    for c, scope in zip(constraints, scopes):
        groups.setdefault(find(scope[0]) if scope else None, []).append(c)
    return [groups[k] for k in sorted(groups, key=lambda k: (k is None, k))]


def reduce_constraints(constraints: list[dict], Q_in: np.ndarray) -> dict:
    """Dedup + dominance + decomposition. Never changes the optimum."""
    n_raw = len(constraints)
    ded, s_ded = _dedup(constraints)
    dom, s_dom = _dominance(ded)
    comps = _components(dom, Q_in) if dom else []
    sizes = [len(g) for g in comps]
    return {
        "components": comps,
        "kept": dom,
        "rule": REDUCTION_RULE,
        "n_raw": n_raw,
        "n_after_dedup": s_ded["n_out"],
        "n_after_dominance": s_dom["n_out"],
        "dedup_removed": n_raw - s_ded["n_out"],
        "dominance_removed": s_ded["n_out"] - s_dom["n_out"],
        "dedup_ratio": (s_ded["n_out"] / n_raw) if n_raw else None,
        "dominance_ratio": ((s_dom["n_out"] / s_ded["n_out"])
                            if s_ded["n_out"] else None),
        "total_ratio": (s_dom["n_out"] / n_raw) if n_raw else None,
        "n_components": len(comps),
        "component_sizes": sorted(sizes, reverse=True),
        "largest_component": (max(sizes) if sizes else 0),
    }


def solve_global(constraints: list[dict], P64: np.ndarray, Q_in: np.ndarray,
                 protected: set, time_limit: float | None,
                 reduce: bool = True) -> dict:
    """ONE segment-wide solve: reduce, then solve each independent component.

    Round-26's core ruling: the excision problem is built from EVERY residual
    transverse row under BOTH diagonals in a single MILP per segment -- no
    event matching, no per-event iteration. Independent components are solved
    separately only because they are provably separable; their objectives sum.

    `time_limit` is a TOTAL budget across all component solves; each solve
    receives what is left of it, and exhausting it before the last component
    is a `limit_no_incumbent` result, never a silent partial mask.
    """
    red = (reduce_constraints(constraints, Q_in) if reduce else
           {"components": ([list(constraints)] if constraints else []),
            "kept": list(constraints), "rule": "reduction disabled",
            "n_raw": len(constraints), "n_after_dedup": len(constraints),
            "n_after_dominance": len(constraints), "dedup_removed": 0,
            "dominance_removed": 0, "dedup_ratio": 1.0,
            "dominance_ratio": 1.0, "total_ratio": 1.0,
            "n_components": (1 if constraints else 0),
            "component_sizes": ([len(constraints)] if constraints else []),
            "largest_component": len(constraints)})

    records: list[dict] = []
    chosen: list = []
    lexicographic = True
    statuses: list[str] = []
    scipy_version = None
    t0 = time.perf_counter()
    for i, group in enumerate(red["components"]):
        left = (None if time_limit is None
                else time_limit - (time.perf_counter() - t0))
        if left is not None and left <= 0.0:
            records.append({"stage": 1, "component": i,
                            "purpose": "min_excised_area",
                            "status": "limit_no_incumbent",
                            "raw_status": 1,
                            "message": ("total MILP time budget exhausted "
                                        "before this component was solved"),
                            "objective": None, "dual_bound": None,
                            "mip_gap": None, "solve_time_s": 0.0,
                            "n_constraints": len(group)})
            statuses.append("limit_no_incumbent")
            break
        sol = _solve_milp(group, P64, Q_in, protected, left)
        for rec in sol["records"]:
            rec["component"] = i
            rec["component_n_constraints"] = len(group)
        records.extend(sol["records"])
        statuses.append(sol["status"])
        lexicographic = lexicographic and sol["lexicographic"]
        scipy_version = sol["scipy_version"]
        if sol["status"] not in ("optimal", "best_found"):
            break
        chosen.extend(sol["chosen"])

    if not red["components"]:
        status = "optimal"
    else:
        bad = [s for s in statuses if s not in ("optimal", "best_found")]
        if bad:
            status = bad[0]
        elif len(statuses) < len(red["components"]):
            status = "limit_no_incumbent"
        else:
            status = ("optimal" if all(s == "optimal" for s in statuses)
                      else "best_found")
    if scipy_version is None:
        import scipy
        scipy_version = scipy.__version__
    return {"status": status, "chosen": sorted(set(chosen)),
            "records": records,
            "lexicographic": (lexicographic and status == "optimal"),
            "scipy_version": scipy_version,
            "reduction": {k: v for k, v in red.items()
                          if k not in ("components", "kept")},
            "total_milp_seconds": time.perf_counter() - t0}


# ==================================== round-27: LP relaxation + 1/k rounding
LP_ROUNDING_RULE = (
    "Round-27 selection strategy, per independent component: solve the "
    "CONTINUOUS relaxation of the same model (x_v, y_q in [0,1]; every "
    "coverage set sums to >= 1; y_q >= x_v for every incident pair; "
    "minimize sum A_q y_q) with scipy.optimize.linprog (HiGHS), then SELECT "
    "every vertex with x_v >= 1/k, where k is the largest coverage-set size "
    "in that component (<= 8). FEASIBLE BY CONSTRUCTION: at most k values "
    "summing to >= 1 have a member >= 1/k, so every coverage constraint "
    "keeps a selected vertex. k-APPROXIMATE ON AREA: every quad charged by "
    "the rounded mask has a corner with x_v >= 1/k, hence its LP y_q >= "
    "1/k, hence A_rounded <= k * A_LP <= k * A_opt. A_LP is a CERTIFIED "
    "LOWER BOUND on that component's minimum excised area, and the "
    "a-posteriori ratio A_achieved / A_LP is reported per component.")

REVERSE_DELETE_RULE = (
    "Reverse-delete post-pass: selected vertices are visited in descending "
    "marginal removed area (ties by grid index) and dropped whenever every "
    "coverage constraint containing them still holds another selected "
    "vertex. Feasibility is invariant by construction -- a vertex is "
    "dropped only when no constraint loses its last cover -- and the "
    "achieved area is non-increasing.")

LOCAL_IMPROVE_RULE = (
    "Bounded local improvement on EXACT marginal removed area: when a "
    "selected vertex v is the last cover of some constraints, any vertex w "
    "lying in all of those constraints' coverage sets is a legal 1-for-1 "
    "swap; it is applied only when the exactly computed freed area (quads "
    "losing their last selected corner) strictly exceeds the exactly "
    "computed added area (quads gaining their first). Runs under a "
    "SEGMENT-WIDE wall-clock budget shared by all components, never a "
    "per-component budget, and is followed by another reverse-delete pass.")

GREEDY_RULE = (
    "Area-aware greedy fallback, used ONLY when the LP itself does not "
    "solve within its cap: repeatedly select the vertex maximizing "
    "(newly covered constraints) / (marginal newly excised area), "
    "zero-marginal-area vertices first. This path is FEASIBLE and nothing "
    "more. The classical ln(n) set-cover guarantee does NOT apply, because "
    "a vertex's cost is the UNION AREA of its incident quads and not an "
    "additive fixed cost, so no approximation claim is attached to it and "
    "components selected this way carry NO lower bound.")

SCHEDULING_NOTE = (
    "Component-size thresholds are a SCHEDULING HEURISTIC -- observed "
    "solve-time behaviour of this MILP family on this data -- and NOT a "
    "hardness claim about the components involved.")

GEOMETRY_STATUS_CLEAN = "transverse_clean_certified"

SELECTION_STATUS_RULE = (
    "geometry_status and selection_status are INDEPENDENT. geometry_status "
    "records whether the emitted arrays, reloaded from disk, censused "
    "transverse-clean under both diagonals. selection_status records how "
    "the cut was CHOSEN: area_optimal (every component proven minimum-area), "
    "lp_rounded (every non-optimal component came from LP rounding, each "
    "with a certified lower bound), heuristic_feasible (greedy fallback "
    "only, no bound), or mixed. The phrase 'minimum-area excision' is "
    "admissible ONLY when selection_status == area_optimal.")

_ROUND_TOL = 1e-9


def quad_area_canonical_grid(P64: np.ndarray) -> np.ndarray:
    """Vectorised canonical per-quad area for the whole grid.

    Same convention as `quad_area_canonical` (mean of the d0 and d1
    two-triangle areas, triangle corner order as in
    `quad_triangle_corners`); pinned equal to the scalar function by the
    test suite. Used so a 400k-quad component is priced in one pass.
    """
    P64 = np.asarray(P64, np.float64)
    p00, p10 = P64[:-1, :-1], P64[1:, :-1]
    p01, p11 = P64[:-1, 1:], P64[1:, 1:]

    def tri(a, b, c):
        return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=-1)

    a0 = tri(p00, p01, p11) + tri(p00, p11, p10)
    a1 = tri(p00, p01, p10) + tri(p01, p11, p10)
    return 0.5 * (a0 + a1)


def _csr(keys: np.ndarray, vals: np.ndarray, n: int):
    """(ptr, idx) grouping `vals` by integer `keys` in [0, n)."""
    order = np.argsort(keys, kind="stable")
    k = keys[order]
    ptr = np.searchsorted(k, np.arange(n + 1), side="left")
    return ptr, vals[order]


class Component:
    """Flat-array view of ONE independent component.

    The shared substrate for the LP relaxation, the 1/k rounding,
    reverse-delete, local improvement and the greedy fallback -- so all of
    them price exactly the same quads with exactly the same areas.
    """

    def __init__(self, constraints: list[dict], P64: np.ndarray,
                 Q_in: np.ndarray, protected: set,
                 area_grid: np.ndarray | None = None):
        nvq, nuq = Q_in.shape
        self.constraints = constraints
        self.n_cons = len(constraints)
        verts = sorted({tuple(map(int, p)) for c in constraints
                        for p in c["coverage"]})
        self.verts = verts
        nvv = len(verts)
        lut = np.full((nvq + 1, nuq + 1), -1, np.int64)
        va = np.asarray(verts, np.int64).reshape(nvv, 2)
        lut[va[:, 0], va[:, 1]] = np.arange(nvv)

        quads = set()
        for v, u in verts:
            for dv in (-1, 0):
                for du in (-1, 0):
                    a, b = v + dv, u + du
                    if 0 <= a < nvq and 0 <= b < nuq and Q_in[a, b]:
                        quads.add((a, b))
        self.quads = sorted(quads)
        nq = len(self.quads)
        qa = (np.asarray(self.quads, np.int64).reshape(nq, 2) if nq
              else np.zeros((0, 2), np.int64))
        if area_grid is None:
            self.areas = np.array([quad_area_canonical(P64, *q)
                                   for q in self.quads], float)
        else:
            self.areas = (np.asarray(area_grid)[qa[:, 0], qa[:, 1]]
                          .astype(float) if nq else np.zeros(0))

        # vertex/quad incidence (only corners that are component vertices)
        if nq:
            cidx = np.stack([lut[qa[:, 0] + dv, qa[:, 1] + du]
                             for dv, du in ((0, 0), (1, 0), (0, 1), (1, 1))],
                            axis=1)
            mask = cidx >= 0
            self.pair_q = np.repeat(np.arange(nq), mask.sum(axis=1))
            self.pair_v = cidx[mask]
        else:
            self.pair_q = np.zeros(0, np.int64)
            self.pair_v = np.zeros(0, np.int64)
        self.vq_ptr, self.vq_idx = _csr(self.pair_v, self.pair_q, nvv)

        # coverage sets, flat in constraint order
        lens = np.array([len(c["coverage"]) for c in constraints], np.int64)
        self.cov_ptr = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
        flat = [tuple(map(int, p)) for c in constraints for p in c["coverage"]]
        fa = (np.asarray(flat, np.int64).reshape(len(flat), 2) if flat
              else np.zeros((0, 2), np.int64))
        self.cov_vidx = (lut[fa[:, 0], fa[:, 1]] if len(flat)
                         else np.zeros(0, np.int64))
        assert (self.cov_vidx >= 0).all(), "coverage vertex outside component"
        self.cov_rows = np.repeat(np.arange(self.n_cons), lens)
        self.vcon_ptr, self.vcon_idx = _csr(self.cov_vidx, self.cov_rows, nvv)
        self.k = int(lens.max()) if self.n_cons else 1

        self.protected = np.zeros(nvv, bool)
        for i, v in enumerate(verts):
            if v in protected:
                self.protected[i] = True
        # a constraint whose whole coverage set is protected can never be hit
        free = ~self.protected[self.cov_vidx]
        self.n_free_per_cons = np.bincount(self.cov_rows[free],
                                           minlength=self.n_cons)
        self.feasible = bool(self.n_cons == 0
                             or (self.n_free_per_cons > 0).all())

    # ---------------------------------------------------------- primitives
    def n_verts(self) -> int:
        return len(self.verts)

    def quad_counts(self, sel: np.ndarray) -> np.ndarray:
        """Per quad: how many of its corners are selected."""
        if not len(self.pair_q):
            return np.zeros(len(self.quads), np.int64)
        return np.bincount(self.pair_q[sel[self.pair_v]],
                           minlength=len(self.quads))

    def cover_counts(self, sel: np.ndarray) -> np.ndarray:
        """Per constraint: how many of its coverage vertices are selected."""
        if not self.n_cons:
            return np.zeros(0, np.int64)
        return np.bincount(self.cov_rows[sel[self.cov_vidx]],
                           minlength=self.n_cons)

    def area(self, sel: np.ndarray) -> float:
        """Excised area of a selection: every quad with a selected corner."""
        qc = self.quad_counts(sel)
        return fsum(self.areas[qc > 0].tolist())

    def is_feasible(self, sel: np.ndarray) -> bool:
        return bool(self.n_cons == 0 or (self.cover_counts(sel) >= 1).all())

    def chosen(self, sel: np.ndarray) -> list:
        return [self.verts[i] for i in np.nonzero(sel)[0]]

    def marginal_areas(self, sel: np.ndarray) -> np.ndarray:
        """Per selected vertex: area freed if it alone were dropped."""
        qc = self.quad_counts(sel)
        out = np.zeros(len(self.verts))
        if not len(self.pair_q):
            return out
        w = np.where(qc[self.pair_q] == 1, self.areas[self.pair_q], 0.0)
        w = w * sel[self.pair_v]
        return np.bincount(self.pair_v, weights=w, minlength=len(self.verts))


# ------------------------------------------------------------------- LP
def lp_relaxation(m: Component, time_limit: float | None = None) -> dict:
    """Continuous relaxation of the component's vertex-mask model.

    Returns the LP status, objective (a CERTIFIED LOWER BOUND on the
    component's minimum excised area when and only when the status is
    `optimal`) and the fractional x vector.
    """
    from scipy import sparse
    from scipy.optimize import linprog

    nvv, nq = len(m.verts), len(m.quads)
    n = nvv + nq
    npair = len(m.pair_v)
    r1, c1 = m.cov_rows, m.cov_vidx
    d1 = -np.ones(len(c1))
    r2 = m.n_cons + np.repeat(np.arange(npair), 2)
    c2 = np.empty(2 * npair, np.int64)
    c2[0::2] = m.pair_v
    c2[1::2] = nvv + m.pair_q
    d2 = np.empty(2 * npair)
    d2[0::2] = 1.0
    d2[1::2] = -1.0
    A = sparse.coo_matrix(
        (np.concatenate([d1, d2]),
         (np.concatenate([r1, r2]), np.concatenate([c1, c2]))),
        shape=(m.n_cons + npair, n)).tocsr()
    b = np.concatenate([-np.ones(m.n_cons), np.zeros(npair)])
    lo = np.zeros(n)
    hi = np.ones(n)
    hi[:nvv][m.protected] = 0.0
    c = np.concatenate([np.zeros(nvv), m.areas])
    options: dict = {"presolve": True}
    if time_limit is not None:
        options["time_limit"] = float(time_limit)
    t0 = time.perf_counter()
    res = linprog(c, A_ub=A, b_ub=b, bounds=np.stack([lo, hi], axis=1),
                  method="highs", options=options)
    dt = time.perf_counter() - t0
    status = {0: "optimal", 1: "limit", 2: "infeasible",
              3: "unbounded"}.get(int(res.status), "numerical_failure")
    x = (np.asarray(res.x[:nvv], float) if getattr(res, "x", None) is not None
         else None)
    return {"status": status, "raw_status": int(res.status),
            "message": str(getattr(res, "message", "")),
            "objective": (float(res.fun) if res.fun is not None else None),
            "x": x, "seconds": dt,
            "n_rows": int(A.shape[0]), "n_cols": n, "n_nnz": int(A.nnz)}


def round_lp(m: Component, x: np.ndarray) -> np.ndarray:
    """1/k rounding: select every vertex with x_v >= 1/k. Feasible by
    construction (see LP_ROUNDING_RULE)."""
    sel = (np.asarray(x, float) >= 1.0 / m.k - _ROUND_TOL) & ~m.protected
    return sel


# -------------------------------------------------------- reverse-delete
def reverse_delete(m: Component, sel: np.ndarray, max_passes: int = 4
                   ) -> np.ndarray:
    """Drop selected vertices whose removal leaves every constraint hit."""
    sel = sel.copy()
    for _ in range(max_passes):
        cnt = m.cover_counts(sel)
        marg = m.marginal_areas(sel)
        idx = np.nonzero(sel)[0]
        order = idx[np.lexsort((idx, -marg[idx]))]
        dropped = 0
        for i in order:
            lo, hi = m.vcon_ptr[i], m.vcon_ptr[i + 1]
            cons = m.vcon_idx[lo:hi]
            if len(cons) and cnt[cons].min() < 2:
                continue
            sel[i] = False
            if len(cons):
                np.subtract.at(cnt, cons, 1)
            dropped += 1
        if not dropped:
            break
    assert m.is_feasible(sel), "reverse-delete broke feasibility"
    return sel


# ----------------------------------------------------- local improvement
def local_improve(m: Component, sel: np.ndarray, deadline: float | None,
                  max_passes: int = 3) -> dict:
    """Bounded 1-for-1 swap search on EXACT marginal removed area."""
    sel = sel.copy()
    swaps = 0
    passes = 0
    stopped = "converged"
    for _ in range(max_passes):
        passes += 1
        cnt = m.cover_counts(sel)
        qc = m.quad_counts(sel)
        crit = np.nonzero(cnt == 1)[0]
        if not len(crit):
            stopped = "no_critical_constraints"
            break
        # the sole selected vertex of each critical constraint
        sole: dict = {}
        for cix in crit:
            lo, hi = m.cov_ptr[cix], m.cov_ptr[cix + 1]
            for vi in m.cov_vidx[lo:hi]:
                if sel[vi]:
                    sole.setdefault(int(vi), []).append(int(cix))
                    break
        applied = 0
        for vi, cons in sorted(sole.items()):
            if deadline is not None and time.perf_counter() > deadline:
                stopped = "budget_exhausted"
                break
            need = len(cons)
            tally: Counter = Counter()
            for cix in cons:
                lo, hi = m.cov_ptr[cix], m.cov_ptr[cix + 1]
                tally.update(int(w) for w in m.cov_vidx[lo:hi])
            cands = [w for w, n in tally.items()
                     if n == need and w != vi and not m.protected[w]
                     and not sel[w]]
            if not cands:
                continue
            qv = set(m.vq_idx[m.vq_ptr[vi]:m.vq_ptr[vi + 1]].tolist())
            freed_all = {q for q in qv if qc[q] == 1}
            best = (0.0, None)
            for w in cands:
                qw = set(m.vq_idx[m.vq_ptr[w]:m.vq_ptr[w + 1]].tolist())
                freed = fsum(m.areas[q] for q in freed_all - qw)
                added = fsum(m.areas[q] for q in qw if qc[q] == 0)
                if freed - added > best[0] + 1e-15:
                    best = (freed - added, w)
            if best[1] is None:
                continue
            w = best[1]
            # revalidate against the CURRENT selection before applying
            lo, hi = m.vcon_ptr[vi], m.vcon_ptr[vi + 1]
            ok = True
            for cix in m.vcon_idx[lo:hi]:
                if cnt[cix] > 1:
                    continue
                a, b = m.cov_ptr[cix], m.cov_ptr[cix + 1]
                if w not in m.cov_vidx[a:b]:
                    ok = False
                    break
            if not ok:
                continue
            sel[vi] = False
            sel[w] = True
            np.subtract.at(cnt, m.vcon_idx[m.vcon_ptr[vi]:m.vcon_ptr[vi + 1]],
                           1)
            np.add.at(cnt, m.vcon_idx[m.vcon_ptr[w]:m.vcon_ptr[w + 1]], 1)
            np.subtract.at(qc, m.vq_idx[m.vq_ptr[vi]:m.vq_ptr[vi + 1]], 1)
            np.add.at(qc, m.vq_idx[m.vq_ptr[w]:m.vq_ptr[w + 1]], 1)
            applied += 1
        swaps += applied
        if stopped == "budget_exhausted" or not applied:
            if stopped != "budget_exhausted":
                stopped = "converged"
            break
    sel = reverse_delete(m, sel)
    assert m.is_feasible(sel), "local improvement broke feasibility"
    return {"sel": sel, "n_swaps": swaps, "passes": passes,
            "stopped": stopped}


# ------------------------------------------------------- greedy fallback
def greedy_cover(m: Component) -> dict:
    """Area-aware greedy. FEASIBLE ONLY -- see GREEDY_RULE."""
    import heapq
    nvv = len(m.verts)
    sel = np.zeros(nvv, bool)
    if not m.n_cons:
        return {"sel": sel, "n_selected": 0, "iterations": 0}
    covered = np.zeros(m.n_cons, bool)
    qcov = np.zeros(len(m.quads), bool)

    def score(i):
        lo, hi = m.vcon_ptr[i], m.vcon_ptr[i + 1]
        cons = m.vcon_idx[lo:hi]
        gain = int((~covered[cons]).sum()) if len(cons) else 0
        if gain == 0:
            return 0.0, 0.0
        qs = m.vq_idx[m.vq_ptr[i]:m.vq_ptr[i + 1]]
        cost = fsum(m.areas[q] for q in qs if not qcov[q])
        return (float("inf") if cost <= 0.0 else gain / cost), cost

    heap = []
    for i in range(nvv):
        if m.protected[i]:
            continue
        s, _ = score(i)
        if s > 0.0:
            heapq.heappush(heap, (-s, i))
    it = 0
    while not covered.all() and heap:
        neg, i = heapq.heappop(heap)
        s, _ = score(i)
        if s <= 0.0:
            continue
        if s < -neg - 1e-12:            # stale key: re-insert and retry
            heapq.heappush(heap, (-s, i))
            continue
        sel[i] = True
        it += 1
        covered[m.vcon_idx[m.vcon_ptr[i]:m.vcon_ptr[i + 1]]] = True
        qcov[m.vq_idx[m.vq_ptr[i]:m.vq_ptr[i + 1]]] = True
    return {"sel": sel, "n_selected": int(sel.sum()), "iterations": it,
            "covered_all": bool(covered.all())}


# ------------------------------------------------------ selection policy
EXACT_MAX_CONSTRAINTS = 700     # scheduling heuristic, see SCHEDULING_NOTE
EXACT_TIME_LIMIT_S = 30.0
EXACT_TOTAL_BUDGET_S = 600.0    # segment-wide, cheapest components first
STAGE2_MAX_CONSTRAINTS = 50
IMPROVE_BUDGET_S = 90.0
LP_TIME_LIMIT_S = 600.0


def _component_select(m: Component, lp_time_limit: float | None) -> dict:
    """Phase A for one component: LP -> 1/k rounding -> reverse-delete,
    with the area-aware greedy as the fallback when the LP does not solve.

    Feasibility NEVER depends on exact effort: this returns a feasible mask
    for every component, whatever the solver did.
    """
    rec: dict = {"n_constraints": m.n_cons, "n_vertex_vars": len(m.verts),
                 "n_quad_vars": len(m.quads), "k": m.k,
                 "lp_seconds": 0.0, "round_seconds": 0.0,
                 "improve_seconds": 0.0, "exact_seconds": 0.0,
                 "greedy_seconds": 0.0}
    if not m.feasible:
        rec.update(method="infeasible", lower_bound=None,
                   lower_bound_kind="none",
                   note=("a coverage set is entirely protected: no mask can "
                         "hit it"))
        return rec
    lp = lp_relaxation(m, lp_time_limit)
    rec["lp_seconds"] = lp["seconds"]
    rec["lp_status"] = lp["status"]
    rec["lp_rows"] = lp["n_rows"]
    rec["lp_cols"] = lp["n_cols"]
    rec["lp_nnz"] = lp["n_nnz"]
    if lp["status"] == "optimal" and lp["x"] is not None:
        t = time.perf_counter()
        sel = round_lp(m, lp["x"])
        raw_area = m.area(sel)
        feasible_raw = m.is_feasible(sel)
        if not feasible_raw:                    # never observed; never silent
            rec["rounding_infeasible"] = True
        else:
            sel = reverse_delete(m, sel)
        rec["round_seconds"] = time.perf_counter() - t
        if feasible_raw:
            rec.update(method="lp_rounded", sel=sel,
                       lower_bound=lp["objective"],
                       lower_bound_kind="lp_objective",
                       area_raw_rounded=raw_area,
                       n_selected_raw_rounded=int(round_lp(m,
                                                           lp["x"]).sum()))
            return rec
    t = time.perf_counter()
    g = greedy_cover(m)
    rec["greedy_seconds"] = time.perf_counter() - t
    sel = reverse_delete(m, g["sel"])
    rec.update(method="greedy_feasible", sel=sel, lower_bound=None,
               lower_bound_kind="none",
               greedy_note=GREEDY_RULE,
               lp_fallback_reason=f"LP status {lp['status']}")
    return rec


def select_global(constraints: list[dict], P64: np.ndarray, Q_in: np.ndarray,
                  protected=(), *,
                  lp_time_limit: float | None = LP_TIME_LIMIT_S,
                  lp_total_budget: float | None = None,
                  exact_max_constraints: int = EXACT_MAX_CONSTRAINTS,
                  exact_time_limit: float = EXACT_TIME_LIMIT_S,
                  exact_total_budget: float | None = EXACT_TOTAL_BUDGET_S,
                  stage2_max_constraints: int = STAGE2_MAX_CONSTRAINTS,
                  improve_budget: float = IMPROVE_BUDGET_S,
                  reduce: bool = True,
                  area_grid: np.ndarray | None = None) -> dict:
    """ONE segment-wide selection under the round-27 policy.

    Phase A -- every component gets an LP-rounded (or, if its LP does not
    solve, a greedy) FEASIBLE mask, followed by reverse-delete. Feasibility
    never depends on exact effort.
    Phase B -- components with at most `exact_max_constraints` constraints
    attempt the stage-1 exact MILP under a SHORT per-component cap; a
    PROVEN-OPTIMAL result replaces the rounded one (a time-limited
    incumbent does not: it would buy nothing and muddy the vocabulary).
    Stage 2 (vertex-count minimization) runs only under
    `stage2_max_constraints`.
    Phase C -- bounded local improvement on the components that were not
    proven optimal, under a SEGMENT-WIDE budget.
    """
    import scipy
    protected = {tuple(map(int, p)) for p in protected}
    t_start = time.perf_counter()
    red = (reduce_constraints(constraints, Q_in) if reduce else
           {"components": ([list(constraints)] if constraints else []),
            "kept": list(constraints), "rule": "reduction disabled",
            "n_raw": len(constraints), "n_after_dedup": len(constraints),
            "n_after_dominance": len(constraints), "dedup_removed": 0,
            "dominance_removed": 0, "dedup_ratio": 1.0,
            "dominance_ratio": 1.0, "total_ratio": 1.0,
            "n_components": (1 if constraints else 0),
            "component_sizes": ([len(constraints)] if constraints else []),
            "largest_component": len(constraints)})
    groups = red["components"]
    if area_grid is None and groups:
        area_grid = quad_area_canonical_grid(P64)

    t = time.perf_counter()
    models = [Component(g, P64, Q_in, protected, area_grid) for g in groups]
    t_model = time.perf_counter() - t

    # ---- phase A: a feasible mask for EVERY component --------------------
    t_a = time.perf_counter()
    recs: list[dict] = []
    lp_spent = 0.0
    for i, m in enumerate(models):
        left = lp_time_limit
        if lp_total_budget is not None:
            rem = lp_total_budget - lp_spent
            left = rem if left is None else min(left, rem)
            if left <= 0.0:
                left = 0.0
        rec = _component_select(m, (None if left is None else max(left, 0.0)))
        lp_spent += rec.get("lp_seconds", 0.0)
        rec["component"] = i
        recs.append(rec)
    t_phase_a = time.perf_counter() - t_a

    if any(r["method"] == "infeasible" for r in recs):
        return {"status": "infeasible", "chosen": [],
                "components": [{k: v for k, v in r.items() if k != "sel"}
                               for r in recs],
                "reduction": {k: v for k, v in red.items()
                              if k not in ("components", "kept")},
                "scipy_version": scipy.__version__,
                "label": ("a coverage set lies entirely inside the protected "
                          "vertex set; no mask can satisfy it")}

    # ---- phase B: exact MILP where scheduling says it is cheap -----------
    # Cheapest components FIRST, under a SEGMENT-WIDE budget: on a segment
    # with 800+ components a per-component cap alone would spend hours
    # proving small optima and never reach the recensus. Skipping is
    # recorded per component; it costs optimality claims, never feasibility.
    t_b = time.perf_counter()
    exact_deadline = (time.perf_counter() + exact_total_budget
                      if exact_total_budget else None)
    milp_records: list[dict] = []
    order_b = sorted(range(len(models)), key=lambda i: (models[i].n_cons, i))
    for i in order_b:
        m, rec = models[i], recs[i]
        if m.n_cons == 0 or m.n_cons > exact_max_constraints:
            continue
        if exact_deadline is not None \
                and time.perf_counter() > exact_deadline:
            rec["exact_status"] = "not_attempted_budget_exhausted"
            continue
        t = time.perf_counter()
        sol = _solve_milp(m.constraints, P64, Q_in, protected,
                          exact_time_limit,
                          stage2=(m.n_cons <= stage2_max_constraints))
        rec["exact_seconds"] = time.perf_counter() - t
        rec["exact_status"] = sol["status"]
        for r in sol["records"]:
            r["component"] = rec["component"]
            r["component_n_constraints"] = m.n_cons
        milp_records.extend(sol["records"])
        stage1 = [r for r in sol["records"] if r["stage"] == 1]
        if stage1 and stage1[0]["dual_bound"] is not None:
            b = float(stage1[0]["dual_bound"])
            if rec["lower_bound"] is None or b > rec["lower_bound"]:
                rec["lower_bound"] = b
                rec["lower_bound_kind"] = (
                    "solver_bound" if sol["status"] != "optimal"
                    else rec["lower_bound_kind"])
        if sol["status"] != "optimal":
            continue
        chosen = {tuple(map(int, v)) for v in sol["chosen"]}
        sel = np.array([v in chosen for v in m.verts], bool)
        if not m.is_feasible(sel):          # never silent
            rec["exact_rejected"] = "MILP solution failed the feasibility re-check"
            continue
        a_exact = m.area(sel)
        a_round = m.area(rec["sel"])
        rec["area_exact"] = a_exact
        rec["area_before_exact"] = a_round
        if a_exact <= a_round + 1e-12:
            rec["sel"] = sel
        rec["method"] = "exact_optimal"
        rec["lower_bound"] = a_exact
        rec["lower_bound_kind"] = "exact_optimum"
        rec["lexicographic"] = bool(sol["lexicographic"])
    t_phase_b = time.perf_counter() - t_b

    # ---- phase C: bounded local improvement, SEGMENT-WIDE budget ---------
    t_c = time.perf_counter()
    deadline = (time.perf_counter() + improve_budget
                if improve_budget and improve_budget > 0 else None)
    order = sorted(
        (i for i, r in enumerate(recs) if r["method"] != "exact_optimal"),
        key=lambda i: -(models[i].area(recs[i]["sel"])
                        - (recs[i]["lower_bound"] or 0.0)))
    for i in order:
        if deadline is not None and time.perf_counter() > deadline:
            recs[i]["improve_status"] = "not_attempted_budget_exhausted"
            continue
        t = time.perf_counter()
        before = models[i].area(recs[i]["sel"])
        imp = local_improve(models[i], recs[i]["sel"], deadline)
        recs[i]["sel"] = imp["sel"]
        recs[i]["improve_seconds"] = time.perf_counter() - t
        recs[i]["improve_status"] = imp["stopped"]
        recs[i]["n_swaps"] = imp["n_swaps"]
        recs[i]["area_before_improve"] = before
    t_phase_c = time.perf_counter() - t_c

    # ---- aggregate --------------------------------------------------------
    chosen: list = []
    bounds = []
    bounded_area = []
    unbounded_area = []
    bound_complete = True
    for m, rec in zip(models, recs):
        sel = rec["sel"]
        assert m.is_feasible(sel), "component mask is not feasible"
        chosen.extend(m.chosen(sel))
        a = m.area(sel)
        rec["achieved_area"] = a
        rec["n_selected"] = int(sel.sum())
        lb = rec["lower_bound"]
        rec["ratio_achieved_over_bound"] = (
            (a / lb) if (lb is not None and lb > 0) else
            (1.0 if (lb is not None and a == 0.0) else None))
        if lb is None:
            bound_complete = False
            unbounded_area.append(a)
        else:
            bounds.append(lb)
            bounded_area.append(a)
    achieved = fsum(r["achieved_area"] for r in recs)
    lower_bound = fsum(bounds) if bounds else 0.0
    a_bounded = fsum(bounded_area) if bounded_area else 0.0
    a_unbounded = fsum(unbounded_area) if unbounded_area else 0.0
    methods = {r["method"] for r in recs}
    if not methods or methods == {"exact_optimal"}:
        selection_status = "area_optimal"
    elif methods <= {"lp_rounded"}:
        selection_status = "lp_rounded"
    elif methods <= {"greedy_feasible"}:
        selection_status = "heuristic_feasible"
    else:
        selection_status = "mixed"
    mix = Counter(r["method"] for r in recs)
    return {
        "status": "ok",
        "chosen": sorted(set(chosen)),
        "selection_status": selection_status,
        "selection_status_rule": SELECTION_STATUS_RULE,
        "method_mix": dict(sorted(mix.items())),
        "achieved_area": achieved,
        "achieved_area_bounded_subset": a_bounded,
        "achieved_area_unbounded_subset": a_unbounded,
        "n_components_without_bound": len(unbounded_area),
        "combined_lower_bound": lower_bound,
        "combined_lower_bound_complete": bound_complete,
        "combined_lower_bound_rule": (
            "sum of the per-component lower bounds (exact optimum, LP "
            "objective, or a valid solver dual bound); valid because "
            "distinct components share no decision variable. Components "
            "selected by the greedy fallback contribute NO bound, and the "
            "sum is then flagged incomplete."),
        # Ratio over the BOUNDED SUBSET ONLY. Dividing the whole achieved
        # area by an incomplete bound would inflate the number by however
        # much area sits in unbounded (greedy) components and is NOT a
        # statement about selection quality; when the bound is incomplete
        # no ratio covers the whole segment and none is offered.
        "ratio_rule": (
            "A_achieved / A_lower_bound over the components that HAVE a "
            "lower bound. When combined_lower_bound_complete is false this "
            "covers only part of the segment: the remaining area "
            "(achieved_area_unbounded_subset) has no bound of any kind and "
            "no ratio may be quoted for it."),
        "ratio_achieved_over_bound": ((a_bounded / lower_bound)
                                      if lower_bound > 0 else None),
        "ratio_covers_area_fraction": ((a_bounded / achieved)
                                       if achieved > 0 else None),
        "k_max": (max(m.k for m in models) if models else None),
        "components": [{k: v for k, v in r.items() if k != "sel"}
                       for r in recs],
        "milp_records": milp_records,
        "reduction": {k: v for k, v in red.items()
                      if k not in ("components", "kept")},
        "rules": {"lp_rounding": LP_ROUNDING_RULE,
                  "reverse_delete": REVERSE_DELETE_RULE,
                  "local_improve": LOCAL_IMPROVE_RULE,
                  "greedy": GREEDY_RULE,
                  "scheduling": SCHEDULING_NOTE},
        "policy": {"lp_time_limit_s": lp_time_limit,
                   "lp_total_budget_s": lp_total_budget,
                   "exact_max_constraints": exact_max_constraints,
                   "exact_time_limit_s": exact_time_limit,
                   "exact_total_budget_s": exact_total_budget,
                   "exact_order": "ascending constraint count, cheapest first",
                   "stage2_max_constraints": stage2_max_constraints,
                   "improve_budget_s": improve_budget},
        "timings": {"model_build_s": t_model, "phase_a_lp_round_s": t_phase_a,
                    "phase_b_exact_s": t_phase_b,
                    "phase_c_improve_s": t_phase_c,
                    "total_s": time.perf_counter() - t_start},
        "scipy_version": scipy.__version__,
    }


# ================== round-28: FROZEN greedy-first scheduling policy
# notes/DECISIONS.md 2026-07-31 round 28 Q2. The prototype (round 27) ran
# the LP FIRST and let it consume 600-1200 s on the giant components, which
# returned neither a mask nor a bound while the greedy produced the same
# component's feasible mask in 0.4-3.7 s. The frozen policy therefore
# builds a feasible incumbent for EVERY component first and treats all
# solver effort as strictly optional improvement.
#
# EVERY NUMBER IN THIS BLOCK IS A SCHEDULING POLICY -- observed solve-time
# behaviour of this MILP family on this data -- and NOT a hardness or
# accuracy claim about the components involved. The block is hashed into
# every certificate so two runs' policies can never be confused.
FROZEN_POLICY_VERSION = "round28-greedy-first-v1"

FROZEN_POLICY: dict = {
    "version": FROZEN_POLICY_VERSION,
    "order": ["greedy_incumbent_for_every_component",
              "reverse_delete_immediately",
              "lp_improvement", "exact_improvement", "local_improvement",
              "reload_and_recensus"],
    "lp_budget_s_per_segment": 60.0,
    "lp_time_limit_s_per_component": 60.0,
    "lp_skip_above_reduced_constraints": 50000,
    "improvement_budget_s_per_segment": 120.0,
    "exact_max_constraints": 700,
    "exact_time_limit_s_per_component": 30.0,
    "stage2_max_constraints": 50,
    "process_limit_s_per_segment": 600.0,
    "replacement_rule": (
        "a candidate replaces the incumbent ONLY when its MEASURED excised "
        "area is strictly lower"),
    "failure_rule": (
        "an optimization failure -- solver error, time limit, infeasible "
        "rounding, exhausted budget -- NEVER removes the feasible "
        "incumbent; it costs an optimality claim, never an artifact"),
    "stage_rule": ("stage 1 (area) only, except components at or under "
                   "stage2_max_constraints"),
    "improvement_budget_split": (
        "the segment-wide improvement budget is spent by the exact MILP "
        "phase first (cheapest components first) and then by bounded local "
        "improvement on whatever remains of it"),
    "scheduling_note": SCHEDULING_NOTE,
}

FROZEN_GREEDY_FIRST_RULE = (
    "Round-28 FROZEN scheduling: (1) an area-aware greedy feasible mask is "
    "built for EVERY component before any solver runs; (2) reverse-delete "
    "is applied immediately; (3) LP and exact improvement are OPTIONAL and "
    "run under segment-wide budgets; (4) a candidate replaces the incumbent "
    "only on strictly lower measured excised area; (5) an optimization "
    "failure never removes the feasible incumbent; (6) the emitted mask is "
    "reloaded from disk and recensused as before. The clean verdict comes "
    "from that recensus and is independent of every number above.")


def frozen_policy_hash(policy: dict | None = None) -> str:
    """sha256 of the canonicalised policy constants block (16 hex chars)."""
    import json as _json
    blob = _json.dumps(policy if policy is not None else FROZEN_POLICY,
                       sort_keys=True, separators=(",", ":"),
                       default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def _greedy_incumbent(m: Component) -> dict:
    """Phase G for one component: area-aware greedy -> reverse-delete.

    This is the ONLY phase that is allowed to fail to produce a mask, and
    it fails only when the component is infeasible by construction (a
    coverage set entirely inside the protected set).
    """
    rec: dict = {"n_constraints": m.n_cons, "n_vertex_vars": len(m.verts),
                 "n_quad_vars": len(m.quads), "k": m.k,
                 "greedy_seconds": 0.0, "lp_seconds": 0.0,
                 "round_seconds": 0.0, "exact_seconds": 0.0,
                 "improve_seconds": 0.0,
                 "lp_attempted": False, "exact_attempted": False,
                 "lower_bound": None, "lower_bound_kind": "none"}
    if not m.feasible:
        rec.update(method="infeasible",
                   note=("a coverage set is entirely protected: no mask can "
                         "hit it"))
        return rec
    t = time.perf_counter()
    g = greedy_cover(m)
    sel = reverse_delete(m, g["sel"])
    rec["greedy_seconds"] = time.perf_counter() - t
    assert m.is_feasible(sel), "greedy incumbent is not feasible"
    rec.update(method="greedy_feasible", sel=sel,
               greedy_incumbent_area=m.area(sel),
               greedy_incumbent_n_selected=int(sel.sum()),
               greedy_construction_seconds=rec["greedy_seconds"],
               greedy_note=GREEDY_RULE)
    return rec


def select_global_frozen(constraints: list[dict], P64: np.ndarray,
                         Q_in: np.ndarray, protected=(), *,
                         policy: dict | None = None,
                         reduce: bool = True,
                         area_grid: np.ndarray | None = None) -> dict:
    """ONE segment-wide selection under the round-28 FROZEN policy.

    Greedy-first: every component holds a feasible incumbent before any
    solver is called, and no later phase may take it away.
    """
    import scipy
    pol = dict(FROZEN_POLICY) if policy is None else dict(policy)
    pol_hash = frozen_policy_hash(pol)
    protected = {tuple(map(int, p)) for p in protected}
    t_start = time.perf_counter()
    red = (reduce_constraints(constraints, Q_in) if reduce else
           {"components": ([list(constraints)] if constraints else []),
            "kept": list(constraints), "rule": "reduction disabled",
            "n_raw": len(constraints), "n_after_dedup": len(constraints),
            "n_after_dominance": len(constraints), "dedup_removed": 0,
            "dominance_removed": 0, "dedup_ratio": 1.0,
            "dominance_ratio": 1.0, "total_ratio": 1.0,
            "n_components": (1 if constraints else 0),
            "component_sizes": ([len(constraints)] if constraints else []),
            "largest_component": len(constraints)})
    groups = red["components"]
    if area_grid is None and groups:
        area_grid = quad_area_canonical_grid(P64)
    t = time.perf_counter()
    models = [Component(g, P64, Q_in, protected, area_grid) for g in groups]
    t_model = time.perf_counter() - t

    # ---- phase G: a feasible incumbent for EVERY component, FIRST --------
    t_g = time.perf_counter()
    recs = [_greedy_incumbent(m) for m in models]
    for i, rec in enumerate(recs):
        rec["component"] = i
    t_phase_g = time.perf_counter() - t_g

    if any(r["method"] == "infeasible" for r in recs):
        return {"status": "infeasible", "chosen": [],
                "components": [{k: v for k, v in r.items() if k != "sel"}
                               for r in recs],
                "reduction": {k: v for k, v in red.items()
                              if k not in ("components", "kept")},
                "policy": pol, "policy_version": pol["version"],
                "policy_hash": pol_hash,
                "scipy_version": scipy.__version__,
                "label": ("a coverage set lies entirely inside the protected "
                          "vertex set; no mask can satisfy it")}

    greedy_total_area = fsum(r["greedy_incumbent_area"] for r in recs)
    greedy_total_seconds = fsum(r["greedy_seconds"] for r in recs)

    # ---- phase LP: OPTIONAL improvement, segment-wide budget -------------
    t_l = time.perf_counter()
    lp_budget = float(pol["lp_budget_s_per_segment"])
    lp_cap = float(pol["lp_time_limit_s_per_component"])
    lp_skip_above = int(pol["lp_skip_above_reduced_constraints"])
    # a zero (or negative) budget is an EXHAUSTED budget, not an absent one:
    # the optional phases must then not run at all
    lp_deadline = time.perf_counter() + lp_budget
    n_lp_improved = 0
    for i in sorted(range(len(models)), key=lambda j: (models[j].n_cons, j)):
        m, rec = models[i], recs[i]
        if m.n_cons == 0:
            rec["lp_skipped_reason"] = "component has no constraints"
            continue
        if m.n_cons > lp_skip_above:
            rec["lp_skipped_reason"] = (
                f"component has {m.n_cons} reduced constraints, above the "
                f"{lp_skip_above} SCHEDULING threshold")
            continue
        left = lp_deadline - time.perf_counter()
        if left <= 0.0:
            rec["lp_skipped_reason"] = "segment-wide LP budget exhausted"
            continue
        rec["lp_attempted"] = True
        t = time.perf_counter()
        try:
            lp = lp_relaxation(m, min(lp_cap, left))
        except Exception as exc:                # never removes the incumbent
            rec["lp_seconds"] = time.perf_counter() - t
            rec["lp_status"] = "error"
            rec["lp_error"] = repr(exc)
            continue
        rec["lp_seconds"] = time.perf_counter() - t
        rec["lp_status"] = lp["status"]
        rec["lp_rows"], rec["lp_cols"] = lp["n_rows"], lp["n_cols"]
        if lp["status"] != "optimal" or lp["x"] is None:
            rec["lp_no_improvement_reason"] = f"LP status {lp['status']}"
            continue
        # the LP objective is a certified lower bound whether or not the
        # rounded mask improves on the incumbent
        if rec["lower_bound"] is None or lp["objective"] > rec["lower_bound"]:
            rec["lower_bound"] = lp["objective"]
            rec["lower_bound_kind"] = "lp_objective"
        t = time.perf_counter()
        cand = round_lp(m, lp["x"])
        rec["area_raw_rounded"] = m.area(cand)
        if not m.is_feasible(cand):             # never observed, never silent
            rec["rounding_infeasible"] = True
            rec["round_seconds"] = time.perf_counter() - t
            continue
        cand = reverse_delete(m, cand)
        rec["round_seconds"] = time.perf_counter() - t
        a_cand, a_inc = m.area(cand), m.area(rec["sel"])
        rec["lp_candidate_area"] = a_cand
        rec["incumbent_area_before_lp"] = a_inc
        if a_cand < a_inc - 1e-12 and m.is_feasible(cand):
            rec["sel"] = cand
            rec["method"] = "lp_improved"
            n_lp_improved += 1
        else:
            rec["lp_no_improvement_reason"] = (
                "LP-rounded candidate did not beat the greedy incumbent")
    t_phase_lp = time.perf_counter() - t_l

    # ---- phase X: OPTIONAL exact improvement, segment-wide budget --------
    t_x = time.perf_counter()
    imp_budget = float(pol["improvement_budget_s_per_segment"])
    imp_deadline = time.perf_counter() + imp_budget      # 0 == exhausted
    exact_max = int(pol["exact_max_constraints"])
    exact_cap = float(pol["exact_time_limit_s_per_component"])
    stage2_max = int(pol["stage2_max_constraints"])
    milp_records: list[dict] = []
    for i in sorted(range(len(models)), key=lambda j: (models[j].n_cons, j)):
        m, rec = models[i], recs[i]
        if m.n_cons == 0:
            continue
        if m.n_cons > exact_max:
            rec["exact_skipped_reason"] = (
                f"component has {m.n_cons} reduced constraints, above the "
                f"{exact_max} SCHEDULING threshold")
            continue
        if time.perf_counter() > imp_deadline:
            rec["exact_skipped_reason"] = (
                "segment-wide improvement budget exhausted")
            continue
        rec["exact_attempted"] = True
        t = time.perf_counter()
        try:
            sol = _solve_milp(m.constraints, P64, Q_in, protected, exact_cap,
                              stage2=(m.n_cons <= stage2_max))
        except Exception as exc:                # never removes the incumbent
            rec["exact_seconds"] = time.perf_counter() - t
            rec["exact_status"] = "error"
            rec["exact_error"] = repr(exc)
            continue
        rec["exact_seconds"] = time.perf_counter() - t
        rec["exact_status"] = sol["status"]
        for r in sol["records"]:
            r["component"] = rec["component"]
            r["component_n_constraints"] = m.n_cons
        milp_records.extend(sol["records"])
        stage1 = [r for r in sol["records"] if r["stage"] == 1]
        if stage1 and stage1[0]["dual_bound"] is not None:
            b = float(stage1[0]["dual_bound"])
            if rec["lower_bound"] is None or b > rec["lower_bound"]:
                rec["lower_bound"] = b
                rec["lower_bound_kind"] = (
                    "solver_bound" if sol["status"] != "optimal"
                    else rec["lower_bound_kind"])
        if sol["status"] != "optimal":
            continue
        chosen_v = {tuple(map(int, v)) for v in sol["chosen"]}
        cand = np.array([v in chosen_v for v in m.verts], bool)
        if not m.is_feasible(cand):             # never silent
            rec["exact_rejected"] = (
                "MILP solution failed the feasibility re-check")
            continue
        a_cand, a_inc = m.area(cand), m.area(rec["sel"])
        rec["exact_candidate_area"] = a_cand
        rec["incumbent_area_before_exact"] = a_inc
        if a_cand < a_inc - 1e-12:
            rec["sel"] = cand                   # strictly lower MEASURED area
        # the component is now PROVEN optimal whichever mask is held: the
        # incumbent either was replaced by the optimum or already matched it
        rec["method"] = "exact_optimal"
        rec["lower_bound"] = a_cand
        rec["lower_bound_kind"] = "exact_optimum"
        rec["lexicographic"] = bool(sol["lexicographic"])
    t_phase_x = time.perf_counter() - t_x

    # ---- phase C: bounded local improvement on what is left of the budget
    t_c = time.perf_counter()
    for i in sorted((j for j, r in enumerate(recs)
                     if r["method"] != "exact_optimal"),
                    key=lambda j: -(models[j].area(recs[j]["sel"])
                                    - (recs[j]["lower_bound"] or 0.0))):
        if time.perf_counter() > imp_deadline:
            recs[i]["improve_status"] = "not_attempted_budget_exhausted"
            continue
        t = time.perf_counter()
        before = models[i].area(recs[i]["sel"])
        try:
            imp = local_improve(models[i], recs[i]["sel"], imp_deadline)
        except Exception as exc:                # never removes the incumbent
            recs[i]["improve_status"] = "error"
            recs[i]["improve_error"] = repr(exc)
            continue
        recs[i]["improve_seconds"] = time.perf_counter() - t
        recs[i]["improve_status"] = imp["stopped"]
        recs[i]["n_swaps"] = imp["n_swaps"]
        recs[i]["area_before_improve"] = before
        a_cand = models[i].area(imp["sel"])
        if a_cand < before - 1e-12 and models[i].is_feasible(imp["sel"]):
            recs[i]["sel"] = imp["sel"]
        else:
            recs[i]["improve_no_improvement"] = True
    t_phase_c = time.perf_counter() - t_c

    # ---- aggregate --------------------------------------------------------
    chosen: list = []
    bounds, bounded_area, unbounded_area = [], [], []
    bound_complete = True
    for m, rec in zip(models, recs):
        sel = rec["sel"]
        assert m.is_feasible(sel), "component mask is not feasible"
        chosen.extend(m.chosen(sel))
        a = m.area(sel)
        rec["achieved_area"] = a
        rec["n_selected"] = int(sel.sum())
        rec["improvement_over_greedy"] = rec["greedy_incumbent_area"] - a
        lb = rec["lower_bound"]
        rec["ratio_achieved_over_bound"] = (
            (a / lb) if (lb is not None and lb > 0) else
            (1.0 if (lb is not None and a == 0.0) else None))
        if lb is None:
            bound_complete = False
            unbounded_area.append(a)
        else:
            bounds.append(lb)
            bounded_area.append(a)
    achieved = fsum(r["achieved_area"] for r in recs)
    lower_bound = fsum(bounds) if bounds else 0.0
    a_bounded = fsum(bounded_area) if bounded_area else 0.0
    a_unbounded = fsum(unbounded_area) if unbounded_area else 0.0
    methods = {r["method"] for r in recs}
    all_optimal = bool(methods) and methods == {"exact_optimal"}
    selection_status = "area_optimal" if all_optimal or not methods else "mixed"
    mix = Counter(r["method"] for r in recs)
    return {
        "status": "ok",
        "chosen": sorted(set(chosen)),
        "policy_version": pol["version"],
        "policy_hash": pol_hash,
        "policy": pol,
        "scheduling_rule": FROZEN_GREEDY_FIRST_RULE,
        "selection_status": selection_status,
        "selection_status_rule": SELECTION_STATUS_RULE,
        "minimum_area_claim_admissible": all_optimal,
        "minimum_area_claim_rule": (
            "the phrase 'minimum-area excision' is admissible ONLY when "
            "EVERY component is proven optimal (selection_status == "
            "area_optimal); a greedy or LP-improved component anywhere in "
            "the segment forbids it"),
        "method_mix": dict(sorted(mix.items())),
        "greedy_incumbent_area": greedy_total_area,
        "greedy_construction_seconds": greedy_total_seconds,
        "n_components_lp_attempted": sum(1 for r in recs if r["lp_attempted"]),
        "n_components_lp_skipped": sum(1 for r in recs
                                       if "lp_skipped_reason" in r),
        "n_components_lp_improved": n_lp_improved,
        "n_components_exact_attempted": sum(1 for r in recs
                                            if r["exact_attempted"]),
        "n_components_exact_skipped": sum(1 for r in recs
                                          if "exact_skipped_reason" in r),
        "achieved_area": achieved,
        "improvement_over_greedy": greedy_total_area - achieved,
        "achieved_area_bounded_subset": a_bounded,
        "achieved_area_unbounded_subset": a_unbounded,
        "n_components_without_bound": len(unbounded_area),
        "combined_lower_bound": lower_bound,
        "combined_lower_bound_complete": bound_complete,
        "combined_lower_bound_rule": (
            "sum of the per-component lower bounds (exact optimum, LP "
            "objective, or a valid solver dual bound); valid because "
            "distinct components share no decision variable. Components "
            "whose LP was skipped or did not solve contribute NO bound, and "
            "the sum is then flagged incomplete."),
        "ratio_rule": (
            "A_achieved / A_lower_bound over the components that HAVE a "
            "lower bound. When combined_lower_bound_complete is false this "
            "covers only part of the segment: the remaining area "
            "(achieved_area_unbounded_subset) has no bound of any kind and "
            "no ratio may be quoted for it."),
        "ratio_achieved_over_bound": ((a_bounded / lower_bound)
                                      if lower_bound > 0 else None),
        "ratio_covers_area_fraction": ((a_bounded / achieved)
                                       if achieved > 0 else None),
        "k_max": (max(m.k for m in models) if models else None),
        "components": [{k: v for k, v in r.items() if k != "sel"}
                       for r in recs],
        "milp_records": milp_records,
        "reduction": {k: v for k, v in red.items()
                      if k not in ("components", "kept")},
        "rules": {"greedy_first": FROZEN_GREEDY_FIRST_RULE,
                  "greedy": GREEDY_RULE,
                  "lp_rounding": LP_ROUNDING_RULE,
                  "reverse_delete": REVERSE_DELETE_RULE,
                  "local_improve": LOCAL_IMPROVE_RULE,
                  "scheduling": SCHEDULING_NOTE},
        "timings": {"model_build_s": t_model,
                    "phase_g_greedy_s": t_phase_g,
                    "phase_lp_s": t_phase_lp,
                    "phase_exact_s": t_phase_x,
                    "phase_improve_s": t_phase_c,
                    "total_s": time.perf_counter() - t_start},
        "scipy_version": scipy.__version__,
    }


# ------------------------------------------------------------- certificate
def _mesh_hash(P: np.ndarray, V: np.ndarray) -> str:
    h = hashlib.sha256()
    h.update(repr(P.shape).encode())
    h.update(str(P.dtype).encode())
    h.update(np.ascontiguousarray(P).tobytes())
    h.update(np.ascontiguousarray(V.astype(np.uint8)).tobytes())
    return h.hexdigest()


def _code_provenance() -> dict:
    """Provenance a public reader can verify: code version, frozen policy
    and a content digest of the published source tree. Never a commit sha
    -- the reader of a release has no repository to resolve one against."""
    from .provenance import release_provenance
    return release_provenance()


def _quad_components(Q: np.ndarray) -> list[list[tuple[int, int]]]:
    """Connected components of the retained-quad complex (shared corners),
    returned as lists of quads so distributions can be reported."""
    qs = [(int(v), int(u)) for v, u in zip(*np.nonzero(Q))]
    if not qs:
        return []
    pos = {q: i for i, q in enumerate(qs)}
    parent = list(range(len(qs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    corner_map: dict = {}
    for v, u in qs:
        for cnr in quad_corners(v, u):
            corner_map.setdefault(cnr, []).append(pos[(v, u)])
    for members in corner_map.values():
        r0 = find(members[0])
        for m in members[1:]:
            rm = find(m)
            if rm != r0:
                parent[rm] = r0
    comps: dict = {}
    for i, q in enumerate(qs):
        comps.setdefault(find(i), []).append(q)
    return list(comps.values())


def _quad_side_edges(v: int, u: int):
    return (((v, u), (v, u + 1)), ((v + 1, u), (v + 1, u + 1)),
            ((v, u), (v + 1, u)), ((v, u + 1), (v + 1, u + 1)))


def _certificate(P_in, P64, V_in, V_out, maxedge, exclude, input_contacts,
                 output_contacts, invalidated, removed, kept, solves,
                 iterations, protected, time_limit, identity,
                 lexicographic, scipy_version, reduction) -> dict:
    inval = set(invalidated)

    # --- removed contacts with destroyed-triangle witnesses --------------
    # Witness scopes (round-24): "triangle_participant" -- the invalidated
    # vertex is a corner of the destroyed crossing triangle itself;
    # "quad_retention" -- it is a corner of the participating quad but not
    # of the crossing triangle (quad-level validity: dropping the quad
    # destroys the triangle all the same).
    out_keys = Counter(c["key"] for c in output_contacts)
    in_keys = Counter(c["key"] for c in input_contacts)
    removed_contacts, all_witnessed = [], True
    for c in input_contacts:
        if out_keys[c["key"]] >= in_keys[c["key"]]:
            continue                    # contact survived (only if not clean)
        diag, v1, u1, v2, u2, k1, k2 = c["key"]
        wit = None
        for (qv, qu), k, tri in (((v1, u1), k1, c["corners1"]),
                                 ((v2, u2), k2, c["corners2"])):
            tri_set = set(tri)
            for cnr in tri:             # prefer a triangle participant
                if cnr in inval:
                    wit = {"invalidated_vertex": list(cnr),
                           "destroyed_triangle": [diag, qv, qu, k],
                           "witness_scope": "triangle_participant"}
                    break
            if wit is None:
                for cnr in quad_corners(qv, qu):
                    if cnr in inval and cnr not in tri_set:
                        wit = {"invalidated_vertex": list(cnr),
                               "destroyed_triangle": [diag, qv, qu, k],
                               "witness_scope": "quad_retention"}
                        break
            if wit:
                break
        if wit is None:
            all_witnessed = False
        removed_contacts.append({"key": list(c["key"]), "witness": wit})

    # --- triangle-identity multisets -------------------------------------
    tri_before = sorted((d, v, u, k) for v, u in
                        [(int(a), int(b)) for a, b in
                         zip(*np.nonzero(retained_quads(P64, V_in, maxedge)))]
                        for d in (0, 1) for k in (0, 1))
    tri_after = sorted((d, v, u, k) for v, u in kept
                       for d in (0, 1) for k in (0, 1))
    after_subset = not (Counter(tri_after) - Counter(tri_before))

    # --- area accounting --------------------------------------------------
    retained_in = sorted(kept) + sorted(removed)
    unresolved = sorted({q for c in output_contacts for q in (c["q1"],
                                                              c["q2"])})
    kept_only = [q for q in kept if q not in set(unresolved)]

    def block(area_fn):
        a_in = fsum(area_fn(*q) for q in retained_in)
        a_clean = fsum(area_fn(*q) for q in kept_only)
        a_exc = fsum(area_fn(*q) for q in removed)
        a_unres = fsum(area_fn(*q) for q in unresolved)
        return {"A_input": a_in, "A_clean": a_clean, "A_excised": a_exc,
                "A_unresolved": a_unres,
                "clean_recovery_fraction": (a_clean / a_in) if a_in else None}

    area = {"convention": ("quad area under diagonal d = sum of its two "
                           "triangle areas, float64, input coordinates; "
                           "canonical = mean(d0, d1)"),
            "d0": block(lambda v, u: quad_area(P64, v, u, 0)),
            "d1": block(lambda v, u: quad_area(P64, v, u, 1)),
            "canonical": block(lambda v, u: quad_area_canonical(P64, v, u))}

    # --- cut boundary -----------------------------------------------------
    kept_edges = {tuple(sorted(e)) for q in kept for e in _quad_side_edges(*q)}
    boundary = sorted({tuple(sorted(e)) for q in removed
                       for e in _quad_side_edges(*q)}
                      & kept_edges)
    boundary_length = fsum(
        float(np.linalg.norm(P64[a] - P64[b])) for a, b in boundary)

    # --- component distributions ------------------------------------------
    def comp_summary(Q):
        comps = _quad_components(Q)
        dist = sorted(
            ({"n_quads": len(comp),
              "area_canonical": fsum(quad_area_canonical(P64, v, u)
                                     for v, u in comp)} for comp in comps),
            key=lambda d: (-d["area_canonical"], -d["n_quads"]))
        return len(comps), dist

    n_before, dist_before = comp_summary(retained_quads(P64, V_in, maxedge))
    n_after, dist_after = comp_summary(retained_quads(P64, V_out, maxedge))

    d0_in = sum(1 for c in input_contacts if c["key"][0] == 0)
    d1_in = sum(1 for c in input_contacts if c["key"][0] == 1)
    d0_out = sum(1 for c in output_contacts if c["key"][0] == 0)
    d1_out = sum(1 for c in output_contacts if c["key"][0] == 1)

    stage1 = [s for s in solves if s["stage"] == 1]
    solver_status = ("not_required" if not solves else
                     ("optimal" if all(s["status"] == "optimal"
                                       for s in stage1) else "best_found"))
    solver_block = {
        "backend": "scipy.optimize.milp (HiGHS)",
        "scipy_version": scipy_version,
        "status": solver_status,
        "lexicographic": (None if not solves else lexicographic),
        "objective_policy": ("two-stage lexicographic: stage 1 minimizes "
                             "excised area only; stage 2 minimizes "
                             "invalidated vertices subject to the recorded "
                             "area cap; tertiary (boundary length, "
                             "component count) is a documented v2"),
        "stage2_allowance_rule": STAGE2_ALLOWANCE_RULE,
        "constraint_scope": (
            "SEGMENT-WIDE (round 26): every censused transverse pair under "
            "BOTH diagonals is a constraint in ONE MILP; there is no event "
            "matching and no per-event iteration"),
        "constraint_reduction": reduction,
        "time_limit_s": time_limit,
        "iterations": iterations,
        "solves": solves,               # per-solve objective/bound/gap
        "final_stage1_objective": (stage1[-1]["objective"] if stage1
                                   else None),
        "final_stage1_dual_bound": (stage1[-1]["dual_bound"] if stage1
                                    else None),
        "optimality_gap": stage1[-1]["mip_gap"] if stage1 else None,
        "protected_vertices": sorted(map(list, protected)),
    }
    if solver_status == "best_found":
        solver_block["note"] = (
            "incumbent-based: the recorded stage-1 objective is the "
            "ACHIEVED area, not a proven minimum, and the cut is not "
            "lexicographically optimized; the clean verdict comes only "
            "from the full recensus of the emitted arrays")

    return {
        "operation": "certified excision",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        # grid shape and dtype are input==output by contract; recorded once
        "grid_shape": list(V_in.shape),
        "dtype": str(P_in.dtype),
        "input": {"sha256": _mesh_hash(P_in, V_in),
                  "n_valid": int(V_in.sum())},
        "output": {"sha256": _mesh_hash(P_in, V_out),
                   "n_valid": int(V_out.sum())},
        "provenance": _code_provenance(),
        "census_params": {"maxedge": maxedge, "exclude": exclude,
                          "diagonals": [0, 1], "arithmetic": "float64",
                          "predicate": ("pure-Python planted-size census: "
                                        "_tri_tri_segment interval method "
                                        "over all retained quad pairs")},
        "solver": solver_block,
        "excision": {"invalidated_vertices": sorted(map(list, invalidated)),
                     "removed_quads": sorted(map(list, removed)),
                     "cut_boundary_edges": [[list(a), list(b)]
                                            for a, b in boundary],
                     "cut_boundary_length": boundary_length},
        "removed_contacts": removed_contacts,
        "input_census": {"d0_contacts": d0_in, "d1_contacts": d1_in,
                         "n_contacts": len(input_contacts)},
        "output_census": {"d0_contacts": d0_out, "d1_contacts": d1_out,
                          "n_contacts": len(output_contacts),
                          "clean_both_diagonals": not output_contacts},
        "triangle_multisets": {
            "identity": "(diag, v, u, triangle_index)",
            "before": [list(t) for t in tri_before],
            "after": [list(t) for t in tri_after],
            "after_subset_of_before": after_subset},
        "contacts": {"n_input_keys": len(input_contacts),
                     "n_output_keys": len(output_contacts),
                     "input_keys": sorted(list(c["key"])
                                          for c in input_contacts),
                     "output_keys": sorted(list(c["key"])
                                           for c in output_contacts),
                     "output_submultiset_of_input":
                         not (out_keys - in_keys),
                     "missing_all_witnessed": all_witnessed},
        "area": area,
        "topology": {
            "components_before": n_before,
            "components_after": n_after,
            "component_distribution_before": dist_before,
            "component_distribution_after": dist_after},
        "validity": {"changes_only_valid_to_invalid":
                         not bool((V_out & ~V_in).any()),
                     "n_invalidated": int((V_in & ~V_out).sum())},
        # COMPUTED (never hardcoded): asserted in excise() over the actual
        # output array's bytes, dtype and shape versus the input
        "invalidation": {
            "carrier": HYBRID_INVALIDATION,
            "mask_cleared_at_excised": bool(identity["mask_cleared"]),
            "coordinates_stamped_missing_at_excised":
                bool(identity["stamped"]),
            "missing_marker": MISSING},
        "retained_coordinate_bit_identity":
            bool(identity["retained_bit_identical"]),
        "coordinates_changed_only_at_excised_cells":
            bool(identity["changed_only_at_excised"]),
        "coordinate_guarantee": RETAINED_BIT_IDENTITY,
        "coordinate_bit_identity_method": (
            "computed: the output points array's dtype and shape are "
            "compared to the input's, the retained cells' bytes are compared "
            "cell-for-cell, and the set of cells whose coordinates differ is "
            "compared to the excised set; all asserted before emission"),
        "staleness_warning": STALENESS_WARNING,
    }


# -------------------------------------------------------------------- main
def excise(P: np.ndarray, V: np.ndarray, *,
           maxedge: float = MAXEDGE_DEFAULT, exclude: int = EXCLUDE_DEFAULT,
           protected=(), time_limit: float | None = None,
           max_iter: int = 8, refuse_shared_support: bool = True,
           reduce: bool = True, strategy: str = "exact",
           selection_options: dict | None = None) -> dict:
    """Certified excision of a planted mesh (see notes/CUTTER-SPEC.md).

    Returns a dict. On success: status "clean" with `points` (the input
    coordinates with the excised cells stamped to the -1 missing marker and
    every retained coordinate bit-identical), `valid` (input mask with
    excised vertices invalidated), and `certificate` -- the round-25 A1
    HYBRID invalidation. Shared-support events refuse with status
    "refused_shared_support" and SHARED_SUPPORT_LABEL unless
    `refuse_shared_support=False`, in which case they are cut and labelled
    JUNCTION_EXCISION_LABEL (round-26 Q1). Solver failure paths
    ("infeasible", "limit_no_incumbent", "solver_error", "iteration_limit",
    "stalled") report evidence and emit NO output arrays and NO clean claim.
    """
    P_in = np.asarray(P)
    P64 = np.asarray(P, np.float64)
    V_in = np.asarray(V, bool)
    Q_in = retained_quads(P64, V_in, maxedge)
    protected = {tuple(map(int, p)) for p in protected}

    input_contacts = census(P64, V_in, maxedge, exclude)

    # third class: shared-support events are refused, no output ----------
    junction = []
    if input_contacts:
        pairs = sorted({(c["q1"], c["q2"]) for c in input_contacts})
        rec = np.array([(v1, u1, v2, u2, 0.0, 0.0)
                        for (v1, u1), (v2, u2) in pairs], dtype=PAIR_DTYPE)
        touching = [e for e in oriented_events(rec) if e["self_touching"]]
        if touching and refuse_shared_support:
            return {"status": "refused_shared_support",
                    "label": SHARED_SUPPORT_LABEL,
                    "events": [{"region_a": sorted(e["region_a"]),
                                "region_b": sorted(e["region_b"])}
                               for e in touching]}
        junction = [{"region_a": sorted(e["region_a"]),
                     "region_b": sorted(e["region_b"])} for e in touching]

    constraints = {c["key"]: c for c in input_contacts}
    solves: list[dict] = []
    lexicographic = True
    scipy_version = None
    reduction: dict = {"n_raw": 0, "n_after_dedup": 0, "n_after_dominance": 0,
                       "n_components": 0, "component_sizes": [],
                       "rule": REDUCTION_RULE}
    invalidated: list = []
    selection: dict | None = None
    V_out = V_in.copy()
    remaining = input_contacts
    iterations = 0
    while remaining:
        if iterations >= max_iter:
            return {"status": "iteration_limit", "solver": {"solves": solves},
                    "label": "no clean claim: iteration budget exhausted",
                    "remaining_contacts": [list(c["key"]) for c in remaining]}
        iterations += 1
        if strategy in ("lp_round", "greedy_first"):
            if strategy == "greedy_first":
                sel = select_global_frozen(list(constraints.values()), P64,
                                           Q_in, protected, reduce=reduce,
                                           **(selection_options or {}))
            else:
                sel = select_global(list(constraints.values()), P64, Q_in,
                                    protected, reduce=reduce,
                                    **(selection_options or {}))
            selection = sel
            sol = {"status": ("optimal" if sel["status"] == "ok"
                              else sel["status"]),
                   "chosen": sel["chosen"],
                   "records": sel.get("milp_records", []),
                   "lexicographic": False,
                   "scipy_version": sel["scipy_version"],
                   "reduction": sel["reduction"]}
        else:
            sol = solve_global(list(constraints.values()), P64, Q_in,
                               protected, time_limit, reduce=reduce)
        for rec in sol["records"]:
            rec["iteration"] = iterations
        solves.extend(sol["records"])
        lexicographic = lexicographic and sol["lexicographic"]
        scipy_version = sol["scipy_version"]
        reduction = sol["reduction"]
        if sol["status"] not in ("optimal", "best_found"):
            return {"status": sol["status"],
                    "solver": {"solves": solves},
                    "label": ("no clean claim: solver returned no usable "
                              "incumbent")}
        invalidated = sorted(sol["chosen"])
        V_out = V_in.copy()
        for v, u in invalidated:
            V_out[v, u] = False
        remaining = census(P64, V_out, maxedge, exclude)
        new = [c for c in remaining if c["key"] not in constraints]
        if remaining and not new:
            # every remaining contact was already constrained: the mask the
            # solver returned cannot satisfy its own constraints -- internal
            # error, never a silent pass
            return {"status": "stalled", "solver": {"solves": solves},
                    "label": ("no clean claim: recensus dirty with no new "
                              "constraints (internal error)"),
                    "remaining_contacts": [list(c["key"]) for c in remaining]}
        for c in new:
            constraints[c["key"]] = c

    inval = set(invalidated)
    removed = []
    kept = []
    for v, u in ((int(a), int(b)) for a, b in zip(*np.nonzero(Q_in))):
        (removed if set(quad_corners(v, u)) & inval else kept).append((v, u))

    # HYBRID invalidation (round-25 A1): clear the mask AND stamp the
    # coordinates. The guarantee is stated over RETAINED cells, and both
    # halves of it are COMPUTED over the actual output array, never
    # hardcoded -- CUTTER-SPEC round-24 amendment 9.3, amended by section 10.
    P_out = P_in.copy()
    stamp = np.zeros(V_in.shape, bool)
    for v, u in invalidated:
        P_out[v, u] = MISSING
        stamp[v, u] = True
    changed = np.any(P_out != P_in, axis=-1)
    excised_cells = stamp
    identity = {
        "retained_bit_identical": bool(
            P_out.dtype == P_in.dtype and P_out.shape == P_in.shape
            and np.ascontiguousarray(P_out[V_out]).tobytes()
            == np.ascontiguousarray(P_in[V_out]).tobytes()),
        "changed_only_at_excised": bool(np.array_equal(changed & ~excised_cells,
                                                       np.zeros_like(changed))),
        "stamped": bool(np.all(P_out[excised_cells] == MISSING)),
        "mask_cleared": bool(not V_out[excised_cells].any()),
    }
    assert identity["retained_bit_identical"], \
        "retained-coordinate bit-identity violated"
    assert identity["changed_only_at_excised"], \
        "coordinates changed outside the excised set"
    assert identity["stamped"] and identity["mask_cleared"], \
        "hybrid invalidation incomplete"

    cert = _certificate(P_in, P64, V_in, V_out, maxedge, exclude,
                        input_contacts, remaining, invalidated, removed,
                        kept, solves, iterations, protected, time_limit,
                        identity, lexicographic, scipy_version, reduction)
    cert["geometry_status"] = (GEOMETRY_STATUS_CLEAN if not remaining
                               else "residual_transverse")
    if selection is not None:
        cert["selection"] = {
            "strategy": strategy,
            "policy_version": selection.get("policy_version"),
            "policy_hash": selection.get("policy_hash"),
            "greedy_incumbent_area": selection.get("greedy_incumbent_area"),
            "greedy_construction_seconds":
                selection.get("greedy_construction_seconds"),
            "improvement_over_greedy":
                selection.get("improvement_over_greedy"),
            "minimum_area_claim_admissible":
                selection.get("minimum_area_claim_admissible"),
            "n_components_lp_attempted":
                selection.get("n_components_lp_attempted"),
            "n_components_lp_skipped":
                selection.get("n_components_lp_skipped"),
            "n_components_exact_attempted":
                selection.get("n_components_exact_attempted"),
            "n_components_exact_skipped":
                selection.get("n_components_exact_skipped"),
            "selection_status": selection["selection_status"],
            "selection_status_rule": SELECTION_STATUS_RULE,
            "method_mix": selection["method_mix"],
            "achieved_area_canonical": selection["achieved_area"],
            "achieved_area_bounded_subset":
                selection["achieved_area_bounded_subset"],
            "achieved_area_unbounded_subset":
                selection["achieved_area_unbounded_subset"],
            "combined_lower_bound": selection["combined_lower_bound"],
            "combined_lower_bound_complete":
                selection["combined_lower_bound_complete"],
            "combined_lower_bound_rule":
                selection["combined_lower_bound_rule"],
            "ratio_achieved_over_bound":
                selection["ratio_achieved_over_bound"],
            "ratio_rule": selection["ratio_rule"],
            "k_max": selection["k_max"],
            "per_component": selection["components"],
            "rules": selection["rules"], "policy": selection["policy"],
            "timings": selection["timings"]}
        cert["selection_status"] = selection["selection_status"]
        cert["solver"]["selection_strategy"] = strategy
        cert["solver"]["status"] = ("area_optimal"
                                    if selection["selection_status"]
                                    == "area_optimal" else "not_proven_optimal")
        cert["solver"]["lexicographic"] = False
        cert["solver"]["note"] = (
            ("The cut was CHOSEN under the round-28 FROZEN greedy-first "
             "policy: an area-aware greedy feasible incumbent for EVERY "
             "component, reverse-deleted immediately, then OPTIONAL LP and "
             "exact improvement under segment-wide budgets which may "
             "replace the incumbent only on strictly lower MEASURED excised "
             "area. `solves` lists only the exact attempts. The clean "
             "verdict comes from the recensus of the emitted arrays and is "
             "independent of selection optimality.")
            if strategy == "greedy_first" else
            ("The cut was CHOSEN by the round-27 deterministic constructor "
             "(LP relaxation + 1/k rounding + reverse-delete + bounded local "
             "improvement), with the exact MILP used only where the "
             "scheduling policy allowed it. `solves` lists only those exact "
             "attempts. The clean verdict comes from the recensus of the "
             "emitted arrays and is independent of selection optimality."))
    else:
        cert["selection_status"] = (
            "area_optimal" if lexicographic or (solves and all(
                s["status"] == "optimal" for s in solves if s["stage"] == 1))
            else "mixed")
    if junction:
        cert["junction_excision"] = {"label": JUNCTION_EXCISION_LABEL,
                                     "events": junction}
    return {"status": "clean",
            "points": P_out,
            "valid": V_out,
            "invalidated_vertices": invalidated,
            "removed_quads": removed,
            "retained_quads": kept,
            "junction_excision": bool(junction),
            "certificate": cert}
