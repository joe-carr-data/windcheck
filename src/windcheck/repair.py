"""From a crossing event to a certified rigid-clearance candidate.

Bridges the spectrum's events (grid regions) to the clearance primitive
(triangles, directions, neighbourhood): extract each side's retained
triangles, propose candidate directions, collect the static neighbourhood
the moving patch could sweep through, and return the best-found rigid
translation per movable side.

The result is a CANDIDATE: acceptance downstream is only ever the census
predicate reporting zero crossings under both triangulations after the
displacement is applied, plus the fidelity/quality gates (round-13 rule).
Candidate directions are principally perpendicular to the pair intersection
lines -- translating parallel to the line where two triangle planes meet
cannot separate the planes.
"""
from __future__ import annotations

from collections import Counter
from time import perf_counter as _pc

import numpy as np

from .clearance import event_rigid_clearance
from .intrinsic import SurfaceGraph


# ---------------------------------------------------------------------------
# Kernel profiler (reviewer spec A4). Pure observation: counters and wall
# clocks around existing calls; NO numerical code path is altered. Near-zero
# overhead (Counter increments + perf_counter pairs around operations that
# each cost >= tens of microseconds). Callers that want per-scope numbers
# (e.g. one transaction) call kernel_profile_reset() at scope entry and
# kernel_profile_snapshot() at scope exit.
#
# The wall clocks are INCLUSIVE: instrumented kernels nest (search_repair
# contains apply_field, local_field_contacts, event_pairs_local and
# quantize_f32; apply_field contains harmonic_field), so a child's seconds
# are counted again inside every ancestor's. The per-key values MUST NOT be
# summed to a total (round-25 hardening) -- the snapshot therefore exposes
# them as "inclusive_seconds" with an explicit "clocks_are_inclusive"
# marker. The profiler is PROCESS-LOCAL module state and is NOT thread-safe:
# concurrent threads would interleave reset/accumulate non-atomically.
# Parallel drivers must use process workers (the W3 corpus driver does).
KERNEL_PROFILE: dict[str, Counter] = {
    "calls": Counter(),      # invocation counts per instrumented kernel
    "seconds": Counter(),    # INCLUSIVE wall seconds per instrumented kernel
    "counts": Counter(),     # auxiliary quantities (cache hits, candidates,
}                            # triangle tests, ...)


def kernel_profile_reset() -> None:
    """Zero every counter (start of a profiling scope)."""
    for c in KERNEL_PROFILE.values():
        c.clear()


def kernel_profile_snapshot() -> dict:
    """Plain-dict copy of the current counters (end of a profiling scope).

    "inclusive_seconds" values nest (clocks_are_inclusive=True): never sum
    them across keys. Process-local, not thread-safe -- see the
    KERNEL_PROFILE note above."""
    return {"calls": dict(KERNEL_PROFILE["calls"]),
            "inclusive_seconds": {
                k: round(float(v), 6)
                for k, v in KERNEL_PROFILE["seconds"].items()},
            "counts": dict(KERNEL_PROFILE["counts"]),
            "clocks_are_inclusive": True}


def region_triangles(g: SurfaceGraph, X: np.ndarray, region) -> list[np.ndarray]:
    """Triangles of RETAINED quads only, bounds-checked.

    Round 14: without the bounds + g.Q check, negative indices wrap to the
    opposite grid edge and maxedge-dropped quads are resurrected from their
    corners -- manufactured triangles in the obstacle set.
    """
    nv, nu = g.Q.shape
    out = []
    for v, u in region:
        if not (0 <= v < nv and 0 <= u < nu and g.Q[v, u]):
            continue
        for t in g.quad_triangles(v, u):
            if (t >= 0).all():
                out.append(X[t])
    return out


def patch_normal(tris: list[np.ndarray]) -> np.ndarray:
    n = np.zeros(3)
    for T in tris:
        n += np.cross(T[1] - T[0], T[2] - T[0])
    L = np.linalg.norm(n)
    return n / L if L > 1e-12 else np.array([0.0, 0.0, 1.0])


def candidate_directions(tris_a, tris_b, sphere: int = 48) -> list[np.ndarray]:
    """Normals of both patches, their sums/differences, the perpendiculars to
    the mean plane-intersection line, plus a DETERMINISTIC sphere sample.

    Slice 1 measured why the sphere matters: two real events had cheap
    clearances (0.16 / 0.48 vx) that the normal-derived set missed entirely
    (round 14 predicted exactly this failure of aggregate normals on curved
    regions). Fixed seed so results are reproducible.
    """
    na, nb = patch_normal(tris_a), patch_normal(tris_b)
    line = np.cross(na, nb)
    cands = [na, -na, nb, -nb, na + nb, na - nb, nb - na, -na - nb]
    if np.linalg.norm(line) > 1e-9:
        line = line / np.linalg.norm(line)
        for base in (na, nb):
            perp = np.cross(line, base)
            if np.linalg.norm(perp) > 1e-9:
                cands += [perp, -perp]
    if sphere:
        rng = np.random.default_rng(7)
        S = rng.normal(size=(sphere, 3))
        cands += list(S / np.linalg.norm(S, axis=1, keepdims=True))
    return [c for c in cands if np.linalg.norm(c) > 1e-9]


def swept_neighbourhood(g: SurfaceGraph, X: np.ndarray, moving_tris,
                        moving, exclude, t_max: float) -> list[np.ndarray]:
    """Static triangles whose AABBs overlap the volume swept by the moving
    patch over ANY direction of length <= t_max (round-14 hole 2).

    Grid vicinity is not a valid collision neighbourhood -- the project's
    core finding is that grid-distant patches share 3D space. The sweep is
    bounded by the fidelity budget t_max: inflate the moving patch's AABB by
    t_max and take every retained quad whose AABB intersects it, excluding
    the moving region and crossing partner (handled as crossing pairs) --
    by grid identity, so shared edges are not self-collisions.
    """
    lo = np.min([T.min(axis=0) for T in moving_tris], axis=0) - t_max
    hi = np.max([T.max(axis=0) for T in moving_tris], axis=0) + t_max
    nv, nu = g.Q.shape
    qv, qu = np.nonzero(g.Q)
    # quad AABBs from the four corners, vectorised
    C = np.stack([g.P[qv, qu], g.P[qv + 1, qu],
                  g.P[qv, qu + 1], g.P[qv + 1, qu + 1]])
    qlo, qhi = C.min(axis=0), C.max(axis=0)
    hit = np.all(qhi >= lo, axis=1) & np.all(qlo <= hi, axis=1)
    skip = set(moving) | set(exclude)
    near = [(int(v), int(u)) for v, u in zip(qv[hit], qu[hit])
            if (int(v), int(u)) not in skip]
    return region_triangles(g, X, near)


def event_clearance(g: SurfaceGraph, X: np.ndarray, ev: dict,
                    voxel_um: float, t_max_vx: float = 8.0) -> dict | None:
    """Best rigid translation over both movable sides, within budget t_max.

    The obstacle set is the swept-AABB spatial neighbourhood, not grid
    vicinity. Returns {side, t_vx, t_mm, direction} for the cheaper side, or
    None when no sampled direction clears either side within t_max_vx (the
    fidelity budget in voxels). Upper bound per sampled directions.
    """
    ta = region_triangles(g, X, ev["region_a"])
    tb = region_triangles(g, X, ev["region_b"])
    if not ta or not tb:
        return None
    dirs = candidate_directions(ta, tb)
    best = None
    for side, mov, stat, mreg, sreg in (("a", ta, tb, ev["region_a"], ev["region_b"]),
                                        ("b", tb, ta, ev["region_b"], ev["region_a"])):
        nb = swept_neighbourhood(g, X, mov, mreg, sreg, t_max_vx)
        r = event_rigid_clearance(mov, stat, dirs,
                                  neighbourhood=[(m, s) for m in mov for s in nb])
        if r and r["t"] <= t_max_vx and (best is None or r["t"] < best["t_vx"]):
            best = {"side": side, "t_vx": r["t"],
                    "t_mm": r["t"] * voxel_um / 1000.0,
                    "direction": r["direction"]}
    return best



def harmonic_field(g: SurfaceGraph, core: set[int], partner: set[int],
                   w_core: float, w_partner: float,
                   support_vx: float = 40.0) -> dict[int, float]:
    """Constrained scalar displacement field over the retained-surface graph.

    Round 15: the radial ring dragged the partner branch along and cancelled
    the repair; and Chebyshev dilation crossed invalid gaps. The field is
    harmonic on the LOCAL SURFACE SUBGRAPH (vertices within surface distance
    support_vx of either branch), with Dirichlet conditions: core at w_core,
    partner at w_partner, the support boundary and beyond at 0. An island
    with no surface path to the core gets no weight, however close its grid
    indices sit.
    """
    from scipy import sparse
    from scipy.sparse.csgraph import dijkstra
    from scipy.sparse.linalg import spsolve

    seeds = sorted(core | partner)
    dfield = dijkstra(g.g, directed=False, indices=seeds, min_only=True,
                      limit=support_vx)
    inside = np.nonzero(np.isfinite(dfield))[0]
    fixed = {**{n: w_core for n in core}, **{n: w_partner for n in partner}}
    interior = [int(n) for n in inside if int(n) not in fixed]
    if not interior:
        return dict(fixed)
    pos = {n: i for i, n in enumerate(interior)}
    indptr, indices = g.g.indptr, g.g.indices
    rows, cols, data = [], [], []
    b = np.zeros(len(interior))
    for i, n in enumerate(interior):
        nbrs = indices[indptr[n]:indptr[n + 1]].tolist()
        rows.append(i); cols.append(i); data.append(float(len(nbrs)))
        for m in nbrs:
            if m in pos:
                rows.append(i); cols.append(pos[m]); data.append(-1.0)
            elif m in fixed:
                b[i] += fixed[m]
            # neighbours outside the support: Dirichlet 0, no contribution
    L = sparse.csr_matrix((data, (rows, cols)),
                          shape=(len(interior), len(interior)))
    w = spsolve(L.tocsc(), b)
    out = dict(fixed)
    out.update({n: float(w[i]) for i, n in enumerate(interior)})
    return out


def merge_events(ev0: dict, ev1: dict) -> dict:
    """Consistently oriented union of the same event under both diagonals.

    Round 16: w094's d1 event spans one more row than its d0 counterpart,
    and a d0-only repair left two d1 crossings. Orientation is chosen by
    region overlap/adjacency: aligned if A0~A1 and B0~B1, swapped if
    crossed; ambiguous overlap raises rather than guesses.
    """
    def near(R, S):
        return sum(1 for a in R for dv in (-1, 0, 1) for du in (-1, 0, 1)
                   if (a[0] + dv, a[1] + du) in S)

    straight = near(ev0["region_a"], ev1["region_a"]) \
        + near(ev0["region_b"], ev1["region_b"])
    swapped = near(ev0["region_a"], ev1["region_b"]) \
        + near(ev0["region_b"], ev1["region_a"])
    if straight == swapped:
        raise ValueError("cannot orient the d0/d1 event union")
    a1, b1 = (ev1["region_a"], ev1["region_b"]) if straight > swapped \
        else (ev1["region_b"], ev1["region_a"])
    return {"region_a": set(ev0["region_a"]) | set(a1),
            "region_b": set(ev0["region_b"]) | set(b1)}


def split_weights(mode) -> tuple[float, float]:
    """Round-20 asymmetric split: wA = lam, wB = -(1-lam). "pinned" = 1.0,
    "symmetric" = 0.5; any float in [0,1] is a valid split."""
    if mode == "pinned":
        lam = 1.0
    elif mode == "symmetric":
        lam = 0.5
    else:
        lam = float(mode)
        if not 0.0 <= lam <= 1.0:
            raise ValueError(f"split {lam} outside [0,1]")
    return lam, -(1.0 - lam)


def split_hi_weight(mode) -> float:
    """max(|w_core|, |w_partner|) for a split -- bitwise the historical
    ``max(lam, 1.0 - lam)`` scalar cap denominator."""
    wa, wb = split_weights(mode)
    return max(abs(wa), abs(wb))


def field_weights(gf: SurfaceGraph, ev: dict, mode, support_vx: float = 40.0,
                  field_cache: dict | None = None) -> dict[int, float]:
    """Harmonic displacement weights for (event, split) on the neutral
    graph `gf`. The field depends on (graph, regions, split, support) but
    NOT on direction or t; both apply_field and the remaining-budget
    interval draw from here so a shared `field_cache` gives them the one
    identical solve. Raises ValueError on a self-touching event exactly
    like apply_field always has."""
    core = {c for q in ev["region_a"] for c in gf.quad_corners(*q)}
    partner = {c for q in ev["region_b"] for c in gf.quad_corners(*q)}
    if core & partner:
        raise ValueError("core and partner share vertices; event is "
                         "self-touching -- not a displacement repair")
    wa, wb = split_weights(mode)
    ckey = (wa, wb, float(support_vx))
    if field_cache is not None and ckey in field_cache:
        KERNEL_PROFILE["counts"]["field_cache_hits"] += 1
        return field_cache[ckey]
    _ts = _pc()
    w = harmonic_field(gf, core, partner, wa, wb, support_vx)
    KERNEL_PROFILE["calls"]["harmonic_field"] += 1
    KERNEL_PROFILE["seconds"]["harmonic_field"] += _pc() - _ts
    if field_cache is not None:
        field_cache[ckey] = w
    return w


# Round-25 Part C policy (recorded): the analytic upper root of the
# remaining-budget quadratic is NOT shrunk by any ULP allowance. A candidate
# whose emitted float32 coordinates land an ULP outside the budget is
# rejected honestly by the authoritative emitted-float32 cumulative gate in
# the executor; the analytic interval is a float64 admissibility bound only.
BUDGET_ULP_ALLOWANCE_VX = 0.0


def field_residual_arrays(w: dict[int, float], gf: SurfaceGraph, P_cur,
                          P_orig) -> tuple[np.ndarray, np.ndarray | None]:
    """(weights, residuals) over the field's NONZERO-weight vertices.

    residuals r_i = P_cur_i - P_orig_i in float64. Zero-weight vertices
    never move, so they impose no budget constraint and are excluded
    entirely (however large their residual). P_orig=None means r_i = 0
    everywhere (returned as residuals=None)."""
    items = [(n, wn) for n, wn in w.items() if wn != 0.0]
    warr = np.array([wn for _, wn in items], np.float64)
    if P_orig is None or P_orig is P_cur or not items:
        return warr, None
    gv, gu = np.nonzero(gf.idx >= 0)
    at = {int(gf.idx[v, u]): (int(v), int(u)) for v, u in zip(gv, gu)}
    A = np.asarray(P_cur, np.float64)
    O = np.asarray(P_orig, np.float64)
    R = np.array([A[at[n]] - O[at[n]] for n, _ in items], np.float64)
    return warr, R


def budget_interval_from_arrays(warr: np.ndarray, R: np.ndarray | None,
                                d: np.ndarray, budget_vx: float,
                                w_hi: float) -> tuple[float, float] | None:
    """Exact remaining per-vertex budget interval for the candidate scalar
    t (round-24 Part C, approved spec).

    Per vertex i with nonzero weight w_i, residual r_i and unit direction
    d, the cumulative displacement after a further relative move t is
    |r_i + w_i t d|; staying within the per-vertex budget B means

        w_i^2 t^2 + 2 w_i (r_i . d) t + (|r_i|^2 - B^2) <= 0.

    The admissible interval is the intersection of every vertex's root
    interval with t >= 0. Returns (lo, hi) or None when empty. Vertices
    with residual EXACTLY zero are bounded analytically by B / w_hi
    (discrete maximum principle: |w_i| <= w_hi = max(|w_core|,
    |w_partner|)), which makes the no-residual case bitwise identical to
    the historical scalar cap budget / max(lam, 1-lam). The upper root is
    NOT shrunk (BUDGET_ULP_ALLOWANCE_VX = 0.0): the emitted-float32
    cumulative gate stays the authoritative rejector at the boundary."""
    B = float(budget_vx)
    if len(warr) == 0:
        return (0.0, B / float(w_hi))
    if R is None:
        return (0.0, B / float(w_hi))
    d = np.asarray(d, np.float64)
    rr = np.einsum("ij,ij->i", R, R)
    zero = rr == 0.0
    lo, hi = 0.0, np.inf
    if zero.any():
        hi = B / float(w_hi)
    m = ~zero
    if m.any():
        wv = warr[m]
        a = wv * wv
        b = 2.0 * wv * (R[m] @ d)
        c = rr[m] - B * B
        disc = b * b - 4.0 * a * c
        if (disc < 0.0).any():
            return None                # some vertex can never re-enter budget
        sq = np.sqrt(disc)
        lo = max(lo, float(((-b - sq) / (2.0 * a)).max()))
        hi = min(hi, float(((-b + sq) / (2.0 * a)).min()))
    if not np.isfinite(hi) or hi < lo:
        return None
    return (lo, hi)


def remaining_budget_interval(w: dict[int, float], gf: SurfaceGraph,
                              P_cur, P_orig, direction,
                              budget_vx: float, mode) -> tuple | None:
    """Convenience wrapper: weights+residuals -> admissible t interval for
    one (direction, split). See budget_interval_from_arrays."""
    d = np.asarray(direction, np.float64)
    d = d / np.linalg.norm(d)
    warr, R = field_residual_arrays(w, gf, P_cur, P_orig)
    return budget_interval_from_arrays(warr, R, d, budget_vx,
                                       split_hi_weight(mode))


def _workspace_base_token(P) -> tuple:
    """Identity token for the apply_field workspace's base mesh: object id
    plus a cheap strided content fingerprint (round-25 hardening). id()
    alone is unsafe (ids recycle after gc); a few strided samples make an
    accidental same-shape collision vanishingly unlikely without touching
    the whole array."""
    A = np.asarray(P)
    n = int(A.size)
    idx = np.linspace(0, n - 1, num=min(16, n)).astype(np.int64)
    samples = tuple(float(A.flat[int(i)]) for i in idx)
    return (id(P), A.shape, str(A.dtype), samples)


def apply_field(g: SurfaceGraph, P: np.ndarray, V: np.ndarray, ev: dict,
                direction, t_vx: float, mode="symmetric",
                support_vx: float = 40.0, gf: SurfaceGraph | None = None,
                field_cache: dict | None = None,
                workspace: dict | None = None) -> tuple[np.ndarray, dict]:
    """Displace an event by t*d RELATIVE motion under a harmonic field.

    mode="pinned": core at +1, partner pinned at 0. mode="symmetric": core
    +1/2, partner -1/2 (round 16's recommendation for w094: same relative
    separation, half the maximum displacement). The field is solved on the
    TRIANGULATION-NEUTRAL perimeter graph (diagonal=-1) so it cannot depend
    on the d0/d1 choice. All arithmetic is float64 regardless of the input
    dtype (real tifxyz is float32 -- round 16's second blocker); the caller
    quantizes exactly once at the certification boundary.

    `workspace` (optional, caller-owned, valid for ONE base mesh P like
    gf/field_cache): patch-not-copy buffer for probe loops. First use
    stores a float64 base copy of P plus a working array; every later call
    resets ONLY the vertices the previous call moved back to their base
    values and applies the new displacements in place, recording the new
    moved-vertex list. Unmoved vertices are bitwise the base values and
    moved vertices see the identical `+=` arithmetic on identical inputs,
    so the result is byte-identical to the fresh-copy path (pinned by
    test_apply_field_workspace). The returned array ALIASES the workspace
    buffer: callers must quantize/copy before the next probe, exactly as
    the existing probe loops already do.
    """
    _t0 = _pc()
    d = np.asarray(direction, float)
    n = np.linalg.norm(d)
    if n < 1e-12:
        raise ValueError("zero-length direction")
    d = d / n
    if gf is None:
        # the deformation graph depends only on (P, V, maxedge): callers
        # probing many (direction, lambda, t) candidates against one base
        # mesh should build it once and pass it in
        gf = SurfaceGraph(np.asarray(P, np.float64), V, diagonal=-1,
                          maxedge=g.maxedge)
    partner = {c for q in ev["region_b"] for c in gf.quad_corners(*q)}
    # the field depends on (graph, regions, split, support) but NOT on
    # direction or t -- a caller-owned cache (valid for one base mesh +
    # one event) collapses the per-probe solve to a lookup; field_weights
    # also validates core/partner disjointness (raises like before)
    w = field_weights(gf, ev, mode, support_vx, field_cache)

    nv, nu = gf.idx.shape
    gv, gu = np.nonzero(gf.idx >= 0)
    at = {int(gf.idx[v, u]): (int(v), int(u)) for v, u in zip(gv, gu)}
    if workspace is None:
        P2 = np.asarray(P, np.float64).copy()
    else:
        # patch-not-copy: reset only the previously moved vertices to the
        # cached base values, then apply the new displacements below with
        # the IDENTICAL in-place arithmetic the fresh-copy path uses.
        # Round-25 hardening: the workspace is valid for ONE base mesh; a
        # base identity token (id + strided content fingerprint) rejects
        # reuse against a different same-shaped base.
        token = _workspace_base_token(P)
        if "base" not in workspace:
            workspace["base"] = np.asarray(P, np.float64).copy()
            workspace["P2"] = workspace["base"].copy()
            workspace["touched"] = []
            workspace["base_token"] = token
        base, P2 = workspace["base"], workspace["P2"]
        assert workspace["base_token"] == token, \
            "apply_field workspace reused against a different base mesh"
        assert base.shape == np.shape(P), "apply_field workspace mesh mismatch"
        for v, u in workspace["touched"]:
            P2[v, u] = base[v, u]
    moved = []
    for node, wn in w.items():
        if wn != 0.0:
            v, u = at[node]
            P2[v, u] += t_vx * wn * d
            moved.append((v, u, wn))
    if workspace is not None:
        workspace["touched"] = [(v, u) for v, u, _ in moved]
    report = {"mode": str(mode), "split_lambda": split_weights(mode)[0],
              "field_nonzero_vertices": len(moved),
              "field_partner_nonzero": sum(1 for c in partner
                                           if w.get(c, 0.0) != 0.0)}
    KERNEL_PROFILE["calls"]["apply_field"] += 1
    KERNEL_PROFILE["seconds"]["apply_field"] += _pc() - _t0
    return P2, report


def quantize32(P2: np.ndarray) -> np.ndarray:
    """float32 quantization at the certification boundary, timed. The
    operation is EXACTLY the historical inline idiom
    ``P2.astype(np.float32).astype(np.float64)``; this helper only adds the
    kernel-profile clock around it."""
    _t0 = _pc()
    out = P2.astype(np.float32).astype(np.float64)
    KERNEL_PROFILE["calls"]["quantize_f32"] += 1
    KERNEL_PROFILE["seconds"]["quantize_f32"] += _pc() - _t0
    return out


def displacement_stats(P_in: np.ndarray, P_out: np.ndarray, V: np.ndarray,
                       maxedge: float) -> dict:
    """Quality gates measured on the ACTUAL final coordinates (round 16:
    never on the requested scalar field), under both triangulations."""
    from .intrinsic import retained_quads
    A = np.asarray(P_in, np.float64)
    B = np.asarray(P_out, np.float64)
    delta = np.linalg.norm(B - A, axis=-1)
    moved = delta > 0
    Qb, Qa = retained_quads(A, V, maxedge), retained_quads(B, V, maxedge)
    out = {
        "quantized_moved_vertices": int(moved.sum()),
        "max_disp_vx": float(delta.max()),
        "rms_disp_vx": float(np.sqrt((delta[moved] ** 2).mean()))
            if moved.any() else 0.0,
        "quads_newly_dropped": int((Qb & ~Qa).sum()),
        "quads_newly_retained": int((~Qb & Qa).sum()),
    }
    nv, nu = V.shape
    mv, mu = np.nonzero(moved)
    quads = {(v + dv, u + du) for v, u in zip(mv, mu)
             for dv in (-1, 0) for du in (-1, 0)}
    quads = {(v, u) for v, u in quads
             if 0 <= v < nv - 1 and 0 <= u < nu - 1 and Qb[v, u] and Qa[v, u]}
    for diag in (0, 1):
        inv = 0
        for v, u in quads:
            combos = (((v, u), (v, u + 1), (v + 1, u + 1)),
                      ((v, u), (v + 1, u + 1), (v + 1, u))) if diag == 0 else \
                     (((v, u), (v, u + 1), (v + 1, u)),
                      ((v, u + 1), (v + 1, u + 1), (v + 1, u)))
            for c in combos:
                n0 = np.cross(A[c[1]] - A[c[0]], A[c[2]] - A[c[0]])
                n1 = np.cross(B[c[1]] - B[c[0]], B[c[2]] - B[c[0]])
                if float(n0 @ n1) < 0:
                    inv += 1
        out[f"triangle_inversions_d{diag}"] = inv
    return out


def event_crossing_pairs(g: SurfaceGraph, X: np.ndarray, ev: dict) -> int:
    """Local predicate: intersecting triangle combos between the regions.

    A region quad missing from g.Q raises -- a repair that drops one of the
    event's own quads is a failed quality gate, never a silent pass.
    """
    from .intrinsic import _tri_tri_segment
    n = 0
    for qa in ev["region_a"]:
        for qb in ev["region_b"]:
            if not (g.Q[qa] and g.Q[qb]):
                raise ValueError(f"event quad dropped from complex: {qa} {qb}")
            for t1 in g.quad_triangles(*qa):
                for t2 in g.quad_triangles(*qb):
                    if (t1 >= 0).all() and (t2 >= 0).all() and \
                       _tri_tri_segment(X[t1], X[t2]) is not None:
                        n += 1
    return n


def _quad_retained(B: np.ndarray, V: np.ndarray, v: int, u: int,
                   maxedge: float) -> bool:
    """Per-quad retention on coordinates alone: four valid corners, six
    edges within maxedge -- the same predicate as retained_quads/selfcross,
    evaluated for one quad without touching the rest of the grid."""
    if not (V[v, u] and V[v + 1, u] and V[v, u + 1] and V[v + 1, u + 1]):
        return False
    if maxedge and maxedge > 0:
        c00, c10 = B[v, u], B[v + 1, u]
        c01, c11 = B[v, u + 1], B[v + 1, u + 1]
        for a, b in ((c00, c01), (c01, c11), (c11, c10), (c00, c10),
                     (c00, c11), (c10, c01)):
            if np.linalg.norm(a - b) > maxedge:
                return False
    return True


def _quad_tris_at(B: np.ndarray, v: int, u: int, diag: int):
    """Triangle corner coordinates of a quad, selfcross's order (identical
    combos to SurfaceGraph.quad_triangles)."""
    if diag == 0:
        cs = (((v, u), (v, u + 1), (v + 1, u + 1)),
              ((v, u), (v + 1, u + 1), (v + 1, u)))
    else:
        cs = (((v, u), (v, u + 1), (v + 1, u)),
              ((v, u + 1), (v + 1, u + 1), (v + 1, u)))
    return [np.array([B[c] for c in combo]) for combo in cs]


def event_pairs_local(B: np.ndarray, V: np.ndarray, ev: dict,
                      maxedge: float, diag: int) -> int:
    """event_crossing_pairs without the SurfaceGraph: retention is a
    per-quad property and the triangle corners come straight from the
    coordinates, so the probe cost is O(|region_a| x |region_b|) instead of
    O(grid). Raises on a dropped region quad exactly like the graph
    version (a repair that drops an event's own quad is a failed gate)."""
    from .intrinsic import _tri_tri_segment
    _t0 = _pc()
    try:
        B = np.asarray(B, np.float64)
        n = 0
        for qa in ev["region_a"]:
            for qb in ev["region_b"]:
                if not (_quad_retained(B, V, *qa, maxedge)
                        and _quad_retained(B, V, *qb, maxedge)):
                    raise ValueError(
                        f"event quad dropped from complex: {qa} {qb}")
                for T1 in _quad_tris_at(B, *qa, diag):
                    for T2 in _quad_tris_at(B, *qb, diag):
                        if _tri_tri_segment(T1, T2) is not None:
                            n += 1
        return n
    finally:
        KERNEL_PROFILE["calls"]["event_pairs_local"] += 1
        KERNEL_PROFILE["seconds"]["event_pairs_local"] += _pc() - _t0


def _both_diag_clean(P32, V, ev, maxedge):
    for diag in (0, 1):
        if event_pairs_local(P32, V, ev, maxedge, diag) != 0:
            return False
    return True


def certified_repair(g: SurfaceGraph, P: np.ndarray, V: np.ndarray,
                     ev: dict, direction, t_vx: float,
                     mode: str = "symmetric", budget_vx: float = 8.0,
                     support_vx: float = 40.0, P_orig=None,
                     budget_point_vx: float | None = None) -> tuple | None:
    """Grow the displacement until the FLOAT32-QUANTIZED mesh is locally
    clean under BOTH triangulations, then bisect the dirty/clean bracket to
    the numerical resolution (a multiple of the local float32 ULP).

    `ev` should be the d0/d1 union from merge_events. The extra displacement
    beyond the LP exit is reported as
    `additional_displacement_to_local_clearance_vx` -- round 16: it is
    geometric (the harmonic repair differs from the rigid LP model), not
    rounding -- with both bracket endpoints recorded.

    Round-25 Part C: with `budget_point_vx` set (the executor passes
    STRICT_POINT_VX=1.0 together with the ORIGINAL mesh `P_orig`), every
    probed t is constrained to the exact remaining per-vertex budget
    interval (see budget_interval_from_arrays), intersected with the
    relative cap `budget_vx` and the rigid lower bound `t_vx`; an empty
    intersection returns None before any probe. With budget_point_vx=None
    (the default) the interval machinery is bypassed entirely and the
    behaviour is bitwise the historical one; omitting P_orig means
    r_i = 0 (the base mesh IS the original).
    """
    corner_grid = {(a, b)
                   for reg in (ev["region_a"], ev["region_b"])
                   for v, u in reg
                   for a, b in ((v, u), (v + 1, u), (v, u + 1),
                                (v + 1, u + 1))
                   if 0 <= a < P.shape[0] and 0 <= b < P.shape[1]}
    pts = np.asarray([P[a, b] for a, b in corner_grid], np.float64)
    ulp = float(np.max(np.spacing(np.abs(pts).astype(np.float32)))) \
        if len(pts) else 1e-3
    res = max(4 * ulp, 1e-3)

    gf = SurfaceGraph(np.asarray(P, np.float64), V, diagonal=-1,
                      maxedge=g.maxedge)
    fcache: dict = {}
    ws: dict = {}

    iv = None
    if budget_point_vx is not None:
        w = field_weights(gf, ev, mode, support_vx, fcache)
        iv = remaining_budget_interval(w, gf, P, P_orig, direction,
                                       budget_point_vx, mode)
        if iv is None or iv[1] < t_vx:
            return None   # empty intersection with the rigid lower bound:
                          # candidate inadmissible before any probe

    def clean_at(t):
        P2, _ = apply_field(g, P, V, ev, direction, t, mode, support_vx,
                            gf=gf, field_cache=fcache, workspace=ws)
        P32 = quantize32(P2)
        return P32, _both_diag_clean(P32, V, ev, g.maxedge)

    t_dirty, t = None, t_vx if iv is None else max(t_vx, iv[0])
    found = None
    for _ in range(60):                       # coarse search
        if t > budget_vx:
            return None
        if iv is not None and t > iv[1]:
            return None                       # remaining budget exhausted
        P32, ok = clean_at(t)
        if ok:
            found = t
            break
        t_dirty, t = t, t + max(res, 0.05 * t_vx)
    if found is None:
        return None
    lo = t_dirty if t_dirty is not None else 0.0
    hi = found
    while hi - lo > res:                      # bracket refinement
        mid = 0.5 * (lo + hi)
        _, ok = clean_at(mid)
        if ok:
            hi = mid
        else:
            lo = mid
    P32, _ = clean_at(hi)
    P2, rep = apply_field(g, P, V, ev, direction, hi, mode, support_vx,
                          gf=gf, field_cache=fcache, workspace=ws)
    P32 = quantize32(P2)
    stats = displacement_stats(P, P32, V, g.maxedge)
    stats.update(rep)
    partner_grid = {(a, b) for v, u in ev["region_b"]
                    for a, b in ((v, u), (v + 1, u), (v, u + 1),
                                 (v + 1, u + 1))
                    if 0 <= a < P.shape[0] and 0 <= b < P.shape[1]}
    quant_partner = sum(
        1 for a, b in partner_grid
        if not np.array_equal(P32[a, b],
                              np.asarray(P, np.float64)[a, b]))
    stats.update({
        "quantized_partner_moved": quant_partner,
        "candidate_lp_exit_vx": t_vx,
        "candidate_diagonal": g.diagonal,
        "applied_relative_vx": hi,
        "additional_displacement_to_local_clearance_vx": hi - t_vx,
        "bracket_dirty_vx": lo,
        "bracket_clean_vx": hi,
        "numerical_resolution_vx": res,
    })
    if iv is not None:
        stats["remaining_budget_interval_vx"] = [iv[0], iv[1]]
        stats["budget_interval_ulp_allowance_vx"] = BUDGET_ULP_ALLOWANCE_VX
    return P32, stats


class BroadphaseCache:
    """Base-mesh broad-phase state for local_field_contacts (W4 item).

    Retention and quad AABBs are strictly per-quad properties of the four
    corner coordinates plus the validity mask. A deformation that moves a
    set of vertices can therefore only change the retention or AABB of the
    quads incident to a moved vertex; every other quad's values under the
    deformed coordinates B are bitwise those of the base mesh A. Building
    the full retained-quad mask and corner-AABB stack once per transaction
    and patching only the touched quads per probe is exactly equivalent to
    the per-call full rebuild.

    A spatial bucket index (CELL = 40.0, the census cell) maps every cell
    to the base-retained quads whose BASE AABB overlaps it, so a probe
    tests only the quads bucketed in the cells its swept AABB overlaps
    (plus the moved quads themselves, whose patched AABBs may have left
    their base cells) instead of scanning all N quad AABBs. The buckets
    are a conservative prefilter only: the exact AABB overlap test is
    still applied to every candidate, so the contact set is unchanged.
    """

    CELL = 40.0

    def __init__(self, P_in: np.ndarray, V: np.ndarray, maxedge: float):
        from .intrinsic import retained_quads
        A = np.asarray(P_in, np.float64)
        self.shape = V.shape
        self.maxedge = float(maxedge) if maxedge else maxedge
        self.Q0 = retained_quads(A, V, maxedge)
        self.qv0, self.qu0 = np.nonzero(self.Q0)
        C = np.stack([A[self.qv0, self.qu0], A[self.qv0 + 1, self.qu0],
                      A[self.qv0, self.qu0 + 1],
                      A[self.qv0 + 1, self.qu0 + 1]])
        self.qlo0, self.qhi0 = C.min(axis=0), C.max(axis=0)
        # lexicographic flat keys, ascending by construction of nonzero
        self.flat0 = self.qv0.astype(np.int64) * (V.shape[1] - 1) \
            + self.qu0.astype(np.int64)
        # ---- spatial bucket index over the BASE AABBs -------------------
        # every base-retained quad i is registered in EVERY cell its base
        # AABB [qlo0[i], qhi0[i]] overlaps; buckets: cell key -> ascending
        # int64 array of indices into flat0/qlo0/qhi0 order
        self.buckets: dict = {}
        n0 = len(self.flat0)
        if n0:
            lo_c = np.floor(self.qlo0 / self.CELL).astype(np.int64)
            hi_c = np.floor(self.qhi0 / self.CELL).astype(np.int64)
            span = hi_c - lo_c                # cells covered - 1, per axis
            keys, idxs = [], []
            for dx in range(int(span[:, 0].max()) + 1):
                for dy in range(int(span[:, 1].max()) + 1):
                    for dz in range(int(span[:, 2].max()) + 1):
                        m = ((span[:, 0] >= dx) & (span[:, 1] >= dy)
                             & (span[:, 2] >= dz))
                        i = np.nonzero(m)[0]
                        if not len(i):
                            continue
                        keys.append(self._cell_key(lo_c[i, 0] + dx,
                                                   lo_c[i, 1] + dy,
                                                   lo_c[i, 2] + dz))
                        idxs.append(i)
            key = np.concatenate(keys)
            idx = np.concatenate(idxs).astype(np.int64)
            order = np.lexsort((idx, key))    # by cell, ascending idx within
            key, idx = key[order], idx[order]
            uk, starts = np.unique(key, return_index=True)
            self.buckets = dict(zip(uk.tolist(), np.split(idx, starts[1:])))

    @staticmethod
    def _cell_key(cx, cy, cz):
        """Pack integer cell coords into one int64 (21 bits per axis)."""
        return (((cx + (1 << 20)) << 42) | ((cy + (1 << 20)) << 21)
                | (cz + (1 << 20)))

    def query(self, blo: np.ndarray, bhi: np.ndarray) -> np.ndarray:
        """Ascending base-quad indices bucketed in cells overlapping
        [blo, bhi]. Complete for every quad whose CURRENT AABB equals its
        base AABB (i.e. all unmoved quads): overlapping AABBs share at
        least one cell. Moved quads must be added by the caller."""
        lo = np.floor(np.asarray(blo, np.float64) / self.CELL).astype(np.int64)
        hi = np.floor(np.asarray(bhi, np.float64) / self.CELL).astype(np.int64)
        hits = []
        for cx in range(int(lo[0]), int(hi[0]) + 1):
            for cy in range(int(lo[1]), int(hi[1]) + 1):
                for cz in range(int(lo[2]), int(hi[2]) + 1):
                    b = self.buckets.get(int(self._cell_key(cx, cy, cz)))
                    if b is not None:
                        hits.append(b)
        if not hits:
            return np.empty(0, dtype=np.int64)
        if len(hits) == 1:
            return hits[0]                    # ascending unique already
        return np.unique(np.concatenate(hits))


def local_field_contacts(P_in: np.ndarray, P32: np.ndarray, V: np.ndarray,
                         ev: dict, maxedge: float,
                         broadphase: BroadphaseCache | None = None) -> list:
    """Round-20 Fix A: test the ACTUAL quantized deformation locally.

    Moved quads = retained quads incident to any vertex whose quantized
    position changed. Obstacles = retained quads whose current AABB overlaps
    a moved quad's before-union-after AABB, excluding grid adjacency
    (Chebyshev <= 1 -- the census's own exclude=1 semantics). Every triangle
    combo is tested under BOTH diagonals on the quantized coordinates, plus
    the target-event pairs themselves. The rigid LP is never consulted here:
    a sheared transition triangle is evaluated as it actually deforms.
    Returns [(diag, moved_quad, obstacle_quad), ...] -- empty means locally
    clean.

    `broadphase` (optional) is a BroadphaseCache built once from the BASE
    mesh P_in: retention and AABBs are then recomputed for only the quads
    incident to a moved vertex and patched over the cached base arrays --
    the patched mask/AABBs equal the full recompute from B exactly (both
    are per-quad functions of the corner coordinates, and unmoved corners
    are bitwise identical between A and B). With broadphase=None the
    behaviour is the original full per-call rebuild.
    """
    from .intrinsic import retained_quads, _tri_tri_segment
    _t0 = _pc()
    _n_bp_cand = 0
    _n_tri_tests = 0
    A = np.asarray(P_in, np.float64)
    B = np.asarray(P32, np.float64)
    nv, nu = V.shape
    delta = np.linalg.norm(B - A, axis=-1)
    mv, mu = np.nonzero(delta > 0)
    moved_q0 = {(v + dv, u + du) for v, u in zip(mv.tolist(), mu.tolist())
                for dv in (-1, 0) for du in (-1, 0)}
    if broadphase is None:
        Q = retained_quads(B, V, maxedge)
    else:
        bp = broadphase
        assert bp.shape == V.shape, "broadphase cache grid mismatch"
        Q = bp.Q0.copy()
        for v, u in moved_q0:                 # only touched quads can differ
            if 0 <= v < nv - 1 and 0 <= u < nu - 1:
                Q[v, u] = _quad_retained(B, V, v, u, maxedge)
    moved_q = {(v, u) for v, u in moved_q0
               if 0 <= v < nv - 1 and 0 <= u < nu - 1 and Q[v, u]}
    qv, qu = np.nonzero(Q)
    if broadphase is None:
        C = np.stack([B[qv, qu], B[qv + 1, qu], B[qv, qu + 1],
                      B[qv + 1, qu + 1]])
        qlo, qhi = C.min(axis=0), C.max(axis=0)
    else:
        # gather base AABBs for base-retained quads, then overwrite every
        # retained touched quad with its AABB under B (this superset covers
        # all newly-retained quads: retention only changes at touched quads)
        flat = qv.astype(np.int64) * (nu - 1) + qu.astype(np.int64)
        n = len(flat)
        qlo, qhi = np.empty((n, 3)), np.empty((n, 3))
        if len(bp.flat0):
            pos = np.minimum(np.searchsorted(bp.flat0, flat),
                             len(bp.flat0) - 1)
            in_base = bp.flat0[pos] == flat
            qlo[in_base] = bp.qlo0[pos[in_base]]
            qhi[in_base] = bp.qhi0[pos[in_base]]
        else:
            in_base = np.zeros(n, dtype=bool)
        patch = ~in_base
        mpos = np.empty(0, dtype=np.int64)    # moved-quad positions in flat
        if moved_q:
            mf = np.asarray(sorted(v * (nu - 1) + u for v, u in moved_q),
                            dtype=np.int64)
            mpos = np.searchsorted(flat, mf)
            patch[mpos] = True
        idx = np.nonzero(patch)[0]
        pv, pu = qv[idx], qu[idx]
        Cp = np.stack([B[pv, pu], B[pv + 1, pu], B[pv, pu + 1],
                       B[pv + 1, pu + 1]])
        qlo[idx], qhi[idx] = Cp.min(axis=0), Cp.max(axis=0)

    def tris(P_arr, v, u, diag):
        if diag == 0:
            cs = (((v, u), (v, u + 1), (v + 1, u + 1)),
                  ((v, u), (v + 1, u + 1), (v + 1, u)))
        else:
            cs = (((v, u), (v, u + 1), (v + 1, u)),
                  ((v, u + 1), (v + 1, u + 1), (v + 1, u)))
        return [np.array([P_arr[a] for a in c]) for c in cs]

    contacts = []
    for v, u in moved_q:
        c0 = np.array([A[v, u], A[v + 1, u], A[v, u + 1], A[v + 1, u + 1]])
        c1 = np.array([B[v, u], B[v + 1, u], B[v, u + 1], B[v + 1, u + 1]])
        blo = np.minimum(c0.min(0), c1.min(0)) - 1e-6
        bhi = np.maximum(c0.max(0), c1.max(0)) + 1e-6
        if broadphase is None:
            m = np.all(qhi >= blo, axis=1) & np.all(qlo <= bhi, axis=1)
            sv, su = qv[m], qu[m]
            _n_bp_cand += len(qv)          # full scan: every quad is tested
        else:
            # bucket prefilter: base-retained quads whose base cells the
            # swept AABB overlaps, mapped to current positions (quads that
            # lost retention under B drop out), PLUS every moved quad --
            # a moved obstacle's patched AABB may have left its base cells.
            # Unmoved quads keep their base AABB bitwise, so the union is
            # a superset of the full-scan hits; the identical AABB test
            # below then reproduces the full scan exactly.
            cand = bp.query(blo, bhi)
            _n_bp_cand += len(cand)
            if len(cand):
                keys = bp.flat0[cand]
                pos = np.searchsorted(flat, keys)
                ok = pos < len(flat)
                pos, keys = pos[ok], keys[ok]
                pos = pos[flat[pos] == keys]
            else:
                pos = np.empty(0, dtype=np.int64)
            cp = np.union1d(pos, mpos)        # ascending, matches scan order
            m = (np.all(qhi[cp] >= blo, axis=1)
                 & np.all(qlo[cp] <= bhi, axis=1))
            sel = cp[m]
            sv, su = qv[sel], qu[sel]
        for a, b_ in zip(sv.tolist(), su.tolist()):
            if max(abs(a - v), abs(b_ - u)) <= 1:
                continue                       # census adjacency exclusion
            if (a, b_) in moved_q and (a, b_) < (v, u):
                continue                       # each moved-moved pair once
            for diag in (0, 1):
                for T1 in tris(B, v, u, diag):
                    for T2 in tris(B, a, b_, diag):
                        _n_tri_tests += 1
                        if _tri_tri_segment(T1, T2) is not None:
                            contacts.append((diag, (v, u), (a, b_)))
    KERNEL_PROFILE["calls"]["local_field_contacts"] += 1
    KERNEL_PROFILE["seconds"]["local_field_contacts"] += _pc() - _t0
    KERNEL_PROFILE["counts"]["broadphase_candidates"] += _n_bp_cand
    KERNEL_PROFILE["counts"]["tri_tri_tests"] += _n_tri_tests
    return contacts


def search_repair(g: SurfaceGraph, P: np.ndarray, V: np.ndarray, ev: dict,
                  budget_point_vx: float = 1.0, support_vx: float = 40.0,
                  lams=(0.5, 0.4, 0.6, 0.3, 0.7, 0.15, 0.85),
                  max_candidates: int = 6, t_steps: int = 5,
                  extra_candidates: int = 0, P_orig=None,
                  rel_cap_vx: float | None = None,
                  skipped_log: list | None = None) -> list[dict]:
    """Round-20 ranked joint direction x split search, locally pre-checked.

    Proposals: per direction, the rigid required displacement over the
    CROSSING pairs only (LP as ranking, per round 20 -- obstacles are
    handled by the actual-field pre-check, not rigid intervals). For each of
    the cheapest directions, each admissible split (|lam|*t and |1-lam|*t
    within the point budget), and a bounded set of displacements from the
    proposal up to the admissible maximum, the ACTUAL harmonic deformation
    is applied, quantized, and locally recensused; locally-clean candidates
    return best-first for the engine oracle (which remains the only
    acceptance authority).

    Round-25 Part C: per (direction, lambda) -- the weights change with
    the split -- the probed t range is the exact remaining per-vertex
    budget interval (budget_interval_from_arrays) against `P_orig`,
    intersected with t>=0, the per-transaction relative cap `rel_cap_vx`
    and the candidate's rigid lower bound `need`. When the caller omits
    P_orig the residuals are zero and the interval upper bound is bitwise
    the historical scalar cap budget_point_vx / max(lam, 1-lam), so the
    default behaviour is unchanged. A combination whose interval is empty
    (or entirely below the rigid lower bound) where the historical scalar
    cap would NOT have skipped it is recorded in `skipped_log` -- the
    candidate is skipped and labelled, never probed.

    `extra_candidates` (round-24 search reuse): collect up to that many
    FURTHER candidates in the same pass, so an oracle-driven caller never
    has to re-run the search. ORDER COMPATIBILITY GUARANTEE: the returned
    list is head + tail, where head is byte-for-byte the list an
    extra_candidates=0 call returns (same candidates, same order) -- a
    candidate goes to the head only while the plain scan would still have
    been running (direction rank < max_candidates*3 AND fewer than
    max_candidates found), and head/tail are sorted separately, so the
    larger budget can never reorder or replace the first
    `max_candidates` results.
    """
    from .clearance import pair_intervals_batched

    KERNEL_PROFILE["calls"]["search_repair"] += 1
    _t0 = _pc()
    ta = region_triangles(g, None if False else _coords(g), ev["region_a"])
    tb = region_triangles(g, _coords(g), ev["region_b"])
    if not ta or not tb:
        return []
    dirs = candidate_directions(ta, tb)
    As = np.array([TA for TA in ta for _ in tb])
    Bs = np.array([TB for _ in ta for TB in tb])
    props = []
    for d in dirs:
        d = d / np.linalg.norm(d)
        _ts = _pc()
        lo, hi, ok = pair_intervals_batched(As, Bs, d)
        KERNEL_PROFILE["calls"]["pair_intervals_batched_crossing"] += 1
        KERNEL_PROFILE["seconds"]["pair_intervals_batched_crossing"] += \
            _pc() - _ts
        inter = ok & (lo <= 0.0) & (hi >= 0.0)
        need = float(hi[inter].max()) if inter.any() else 0.0
        if need > 0:
            props.append((need, d))
    props.sort(key=lambda x: x[0])

    head: list[dict] = []          # what the extra_candidates=0 scan returns
    tail: list[dict] = []          # additional same-pass candidates
    budget = max_candidates + extra_candidates
    gf = SurfaceGraph(np.asarray(P, np.float64), V, diagonal=-1,
                      maxedge=g.maxedge)
    fcache: dict = {}
    ws: dict = {}
    # base broad-phase state (retained mask + quad AABBs) is invariant
    # across all probes of this search: build once, patch per probe
    bpc = BroadphaseCache(P, V, g.maxedge)
    # remaining-budget prep is per LAMBDA (weights change with the split);
    # residual projections then vary per direction only
    iv_prep: dict = {}

    def prep(lam):
        if lam not in iv_prep:
            w = field_weights(gf, ev, lam, support_vx, fcache)
            iv_prep[lam] = field_residual_arrays(w, gf, P, P_orig)
        return iv_prep[lam]

    for pi, (need, d) in enumerate(props[:budget * 3]):
        for lam in lams:
            hi_w = max(lam, 1.0 - lam)
            warr, R = prep(lam)
            iv = budget_interval_from_arrays(warr, R, d, budget_point_vx,
                                             hi_w)
            if iv is not None and rel_cap_vx is not None \
                    and iv[1] > rel_cap_vx:
                iv = (iv[0], float(rel_cap_vx))
                if iv[0] > iv[1]:
                    iv = None
            if iv is None or need > iv[1]:
                # historical scalar-cap skips stay silent (they are not a
                # residual effect); a skip the scalar cap would NOT have
                # made is a remaining-budget exclusion: label it
                if not need > budget_point_vx / hi_w:
                    KERNEL_PROFILE["counts"]["budget_interval_skips"] += 1
                    if skipped_log is not None:
                        skipped_log.append(
                            {"direction": d.tolist(), "lam": lam,
                             "rigid_lower_bound_vx": float(need),
                             "remaining_budget_interval_vx":
                                 None if iv is None else [iv[0], iv[1]],
                             "reason": "empty_remaining_budget_interval"})
                continue
            for t in np.linspace(max(need * 1.001, iv[0]), iv[1], t_steps):
                P2, rep = apply_field(g, P, V, ev, d, float(t), lam,
                                      support_vx, gf=gf, field_cache=fcache,
                                      workspace=ws)
                P32 = quantize32(P2)
                # target event must be clear on both diagonals first
                clear = True
                for diag in (0, 1):
                    try:
                        if event_pairs_local(P32, V, ev, g.maxedge,
                                             diag) != 0:
                            clear = False
                    except ValueError:
                        clear = False          # dropped quad = failed gate
                    if not clear:
                        break
                if not clear:
                    continue
                cts = local_field_contacts(P, P32, V, ev, g.maxedge,
                                           broadphase=bpc)
                if not cts:
                    cand = {"direction": d.tolist(), "lam": lam,
                            "t_rel": float(t),
                            "max_point_vx": float(t) * hi_w,
                            "P32": P32, "field_report": rep}
                    # head only where the plain max_candidates scan would
                    # still have been running; everything else is tail,
                    # marked "extra" so callers can tell the two apart
                    if pi < max_candidates * 3 \
                            and len(head) < max_candidates:
                        head.append(cand)
                    else:
                        cand["extra"] = True
                        tail.append(cand)
                    break                       # cheapest t for this lam
            if len(head) + len(tail) >= budget:
                break
        if len(head) + len(tail) >= budget:
            break
    head.sort(key=lambda c: c["max_point_vx"])
    tail.sort(key=lambda c: c["max_point_vx"])
    KERNEL_PROFILE["seconds"]["search_repair"] += _pc() - _t0
    return head + tail


def _coords(g: SurfaceGraph) -> np.ndarray:
    gv, gu = np.nonzero(g.idx >= 0)
    X = np.empty((g.n, 3))
    X[g.idx[gv, gu]] = np.asarray(g.P, np.float64)[gv, gu]
    return X


def match_events(evs0: list[dict], evs1: list[dict]) -> dict:
    """General cross-diagonal event matcher (round 21: merge_events only
    ORIENTS an already-matched pair; it is not a matcher).

    Score = region adjacency overlap in the better of the two orientations.
    Greedy assignment by descending score; a tie between competing partners
    within 20% relative score is AMBIGUOUS (reported, never guessed).
    Returns {"matched": [(i0, i1, score, swapped)], "unmatched_d0": [...],
    "unmatched_d1": [...], "ambiguous": [(i0, [i1 candidates])]} with
    every score recorded.
    """
    def near(R, S):
        return sum(1 for a in R for dv in (-1, 0, 1) for du in (-1, 0, 1)
                   if (a[0] + dv, a[1] + du) in S)

    # Spatial prefilter: a pair can only score > 0 when some quad of one
    # event lies within Chebyshev 1 of some quad of the other, so bucket
    # evs1 quads and score only 1-ring neighbours. Zero-score pairs never
    # entered `scores`, so the output is identical to the full n0 x n1
    # scan -- this is a pure cost cut for event-rich segments.
    bucket: dict = {}
    for j, e1 in enumerate(evs1):
        for q in set(e1["region_a"]) | set(e1["region_b"]):
            bucket.setdefault(q, set()).add(j)
    scores = {}
    for i, e0 in enumerate(evs0):
        cand = set()
        for q in set(e0["region_a"]) | set(e0["region_b"]):
            for dv in (-1, 0, 1):
                for du in (-1, 0, 1):
                    cand |= bucket.get((q[0] + dv, q[1] + du), set())
        for j in sorted(cand):        # index order, as the full scan had
            e1 = evs1[j]
            straight = near(e0["region_a"], e1["region_a"]) \
                + near(e0["region_b"], e1["region_b"])
            swapped = near(e0["region_a"], e1["region_b"]) \
                + near(e0["region_b"], e1["region_a"])
            s = max(straight, swapped)
            if s > 0:
                scores[(i, j)] = (s, swapped > straight)
    # Round-23 blocker 2: RECIPROCAL matching with total ownership. The
    # greedy one-sided rule could hand one d1 event both to a matched
    # unit and to an ambiguous group. An edge is accepted only when each
    # endpoint is the other's unique best partner outside the 20%
    # ambiguity margin; every remaining connected conflict component is
    # ONE ambiguous group; every event index lands in exactly one bucket.
    by0: dict = {}
    by1: dict = {}
    for (i, j), (s, _) in scores.items():
        by0.setdefault(i, []).append((s, j))
        by1.setdefault(j, []).append((s, i))

    def unique_best(cands):
        cands = sorted(cands, reverse=True)
        if len(cands) > 1 and cands[1][0] >= 0.8 * cands[0][0]:
            return None
        return cands[0][1]

    matched, used0, used1 = [], set(), set()
    for i, cands in by0.items():
        j = unique_best(cands)
        if j is not None and unique_best(by1[j]) == i:
            matched.append((i, j, scores[(i, j)][0], scores[(i, j)][1]))
            used0.add(i)
            used1.add(j)
    # conflict components over score edges among unassigned events
    parent: dict = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (i, j) in scores:
        if i not in used0 and j not in used1:
            a, b = find(("0", i)), find(("1", j))
            if a != b:
                parent[a] = b
    groups: dict = {}
    for (i, j) in scores:
        if i not in used0 and j not in used1:
            groups.setdefault(find(("0", i)), [set(), set()])
            g0, g1 = groups[find(("0", i))]
            g0.add(i)
            g1.add(j)
    ambiguous = [(sorted(g0), sorted(g1)) for g0, g1 in groups.values()]
    in_amb0 = {i for g0, _ in ambiguous for i in g0}
    in_amb1 = {j for _, g1 in ambiguous for j in g1}
    out = {"matched": sorted(matched),
           "unmatched_d0": [i for i in range(len(evs0))
                            if i not in used0 and i not in in_amb0],
           "unmatched_d1": [j for j in range(len(evs1))
                            if j not in used1 and j not in in_amb1],
           "ambiguous": ambiguous}
    n0 = len(out["matched"]) + len(out["unmatched_d0"]) + len(in_amb0)
    n1 = len(out["matched"]) + len(out["unmatched_d1"]) + len(in_amb1)
    assert n0 == len(evs0) and n1 == len(evs1), \
        "matcher ownership invariant violated"
    return out
