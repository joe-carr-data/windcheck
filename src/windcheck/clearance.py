"""Rigid directional clearance: the exact repair-displacement primitive.

For two triangles A (movable) and B (static), the set of translations t*d
that keep them intersecting is {t : t*d in (B - A)}, the Minkowski
difference of convex sets -- an interval in t. Its endpoints come from a
small linear program (round 13's recommended formulation):

    max/min t   s.t.   t*d = sum_j beta_j b_j - sum_i alpha_i a_i,
                       alpha, beta >= 0, each summing to 1.

`t_hi` is the exit distance along +d: displace A by anything beyond it and
the pair no longer intersects. This is EXACT for the declared pair and
direction; optimising over sampled directions is a best-found upper bound on
unrestricted displacement, and must be reported as such ("rigid clearance",
never "exact persistence").

Event-level: a rigid translation of the whole A-side patch by
max-over-pairs t_hi clears every original pair; the swept patch must then be
checked against its static NEIGHBOURHOOD (not only the currently crossing
pairs), avoiding the union of their intersection intervals -- a translation
can clear one crossing and land in another. Acceptance downstream is only
ever the census predicate reporting zero crossings under BOTH triangulations.
"""
from __future__ import annotations

import numpy as np


def pair_interval(A: np.ndarray, B: np.ndarray, d: np.ndarray,
                  ) -> tuple[float, float] | None:
    """[t_lo, t_hi] such that (A + t*d) intersects B, or None if never.

    Closed-form separating-axis computation (W4 item 1): two convex sets
    intersect iff their projections overlap on every axis of the complete
    SAT set -- both face normals, the nine edge-pair cross products, and
    the six in-plane edge-normal axes (needed for coplanar pairs).
    Translation by t*d shifts A's projection linearly, so each axis
    yields a t-interval and the answer is their intersection. Exact for
    non-degenerate triangles; ~10^3x cheaper than the LP, which is kept
    below as `pair_interval_lp`, the reference oracle.
    """
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    d = np.asarray(d, float)
    eA = (A[1] - A[0], A[2] - A[1], A[0] - A[2])
    eB = (B[1] - B[0], B[2] - B[1], B[0] - B[2])
    nA = np.cross(eA[0], eA[1])
    nB = np.cross(eB[0], eB[1])
    axes = [nA, nB]
    axes += [np.cross(a, b) for a in eA for b in eB]
    axes += [np.cross(nA, a) for a in eA]
    axes += [np.cross(nB, b) for b in eB]
    scale = max(float(np.abs(np.concatenate((eA, eB))).max()), 1.0)
    lo, hi = -np.inf, np.inf
    for ax in axes:
        L = float(np.linalg.norm(ax))
        if L < 1e-12 * scale * scale:      # degenerate axis, no information
            continue
        ax = ax / L
        pa = A @ ax
        pb = B @ ax
        s = float(d @ ax)
        c1 = float(pb.min() - pa.max())    # overlap iff c1 <= t*s <= c2
        c2 = float(pb.max() - pa.min())
        if abs(s) < 1e-12:
            if c1 > 0.0 or c2 < 0.0:       # disjoint on this axis for all t
                return None
            continue
        t1, t2 = c1 / s, c2 / s
        if t1 > t2:
            t1, t2 = t2, t1
        lo, hi = max(lo, t1), min(hi, t2)
        if lo > hi:
            return None
    if not (np.isfinite(lo) and np.isfinite(hi)):
        raise ArithmeticError(             # mirror of the LP's unbounded
            "pair_interval unbounded along d")  # status (round 14: never
    return lo, hi                          # silently 'no collision')


def pair_intervals_batched(As: np.ndarray, Bs: np.ndarray, d: np.ndarray):
    """Vectorised pair_interval over N pairs for one direction.

    Returns (lo, hi, feasible): arrays of shape (N,). Pairs with
    feasible=False never intersect along d (the scalar version's None).
    Same SAT axis set; per-axis intervals are scale-invariant so axes are
    not normalised (the division by s cancels the axis length). Raises on
    an unbounded feasible interval, mirroring the scalar/LP behaviour.
    """
    As = np.asarray(As, float).reshape(-1, 3, 3)
    Bs = np.asarray(Bs, float).reshape(-1, 3, 3)
    d = np.asarray(d, float)
    N = len(As)
    if N == 0:
        z = np.zeros(0)
        return z, z, np.zeros(0, bool)
    CHUNK = 100_000
    if N > CHUNK:
        # bounded memory on huge obstacle neighbourhoods (a 1.66M-pair
        # segment allocated multi-GB (N,17) intermediates and OOMed under
        # the per-worker ulimit); per-pair results are independent, so
        # chunking is bit-identical
        parts = [pair_intervals_batched(As[i:i + CHUNK], Bs[i:i + CHUNK], d)
                 for i in range(0, N, CHUNK)]
        return (np.concatenate([p[0] for p in parts]),
                np.concatenate([p[1] for p in parts]),
                np.concatenate([p[2] for p in parts]))
    eA = As[:, [1, 2, 0]] - As                      # (N,3,3)
    eB = Bs[:, [1, 2, 0]] - Bs
    nA = np.cross(eA[:, 0], eA[:, 1])[:, None]      # (N,1,3)
    nB = np.cross(eB[:, 0], eB[:, 1])[:, None]
    ee = np.cross(eA[:, :, None, :], eB[:, None, :, :]).reshape(N, 9, 3)
    axes = np.concatenate([nA, nB, ee,
                           np.cross(nA, eA), np.cross(nB, eB)], axis=1)
    L2 = (axes * axes).sum(2)                       # (N,17) squared length
    scale = np.maximum(np.abs(np.concatenate([eA, eB], axis=1)
                              ).max(axis=(1, 2)), 1.0)
    valid = L2 >= (1e-12 * scale * scale)[:, None] ** 2
    pa = np.einsum("nax,ncx->nac", axes, As)        # (N,17,3)
    pb = np.einsum("nax,ncx->nac", axes, Bs)
    c1 = pb.min(2) - pa.max(2)                      # overlap: c1<=t*s<=c2
    c2 = pb.max(2) - pa.min(2)
    s = axes @ d                                    # (N,17)
    ax_scale = np.sqrt(np.maximum(L2, 1e-300))
    static = np.abs(s) < 1e-12 * ax_scale
    dead = valid & static & ((c1 > 0.0) | (c2 < 0.0))
    infeasible = dead.any(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        t1 = np.where(valid & ~static, c1 / s, -np.inf)
        t2 = np.where(valid & ~static, c2 / s, np.inf)
    tlo = np.minimum(t1, t2)
    thi = np.maximum(t1, t2)
    tlo[~(valid & ~static)] = -np.inf
    thi[~(valid & ~static)] = np.inf
    lo = tlo.max(axis=1)
    hi = thi.min(axis=1)
    feasible = ~infeasible & (lo <= hi)
    if np.any(feasible & (~np.isfinite(lo) | ~np.isfinite(hi))):
        raise ArithmeticError("pair_intervals_batched: unbounded along d")
    return lo, hi, feasible


def pair_interval_lp(A: np.ndarray, B: np.ndarray, d: np.ndarray,
                     ) -> tuple[float, float] | None:
    """Reference oracle for pair_interval: two LPs over the Minkowski
    difference. Retained for regression only (the analytic SAT version
    above is the production path)."""
    from scipy.optimize import linprog

    d = np.asarray(d, float)
    # variables: t, alpha(3), beta(3)
    A_eq = np.zeros((5, 7))
    A_eq[:3, 0] = -d
    for i in range(3):
        A_eq[:3, 1 + i] = -A[i]
        A_eq[:3, 4 + i] = B[i]
    A_eq[3, 1:4] = 1.0
    A_eq[4, 4:7] = 1.0
    b_eq = np.array([0.0, 0.0, 0.0, 1.0, 1.0])
    bounds = [(None, None)] + [(0.0, 1.0)] * 6
    out = []
    for sign in (1.0, -1.0):
        c = np.zeros(7)
        c[0] = sign
        r = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
        if r.status == 2:            # infeasible: the pair never intersects
            return None              # along this line -- a real answer
        if not r.success:            # numerical/unbounded failure is NOT
            raise ArithmeticError(   # "no collision" (round 14)
                f"pair_interval LP failed: status {r.status} {r.message}")
        out.append(float(r.x[0]))
    lo, hi = sorted(out)
    return lo, hi


def clear_time(required: float, blocked: list[tuple[float, float]],
               margin: float = 1e-9) -> float:
    """Smallest t >= required outside every blocking interval."""
    t = required
    changed = True
    while changed:
        changed = False
        for lo, hi in blocked:
            if lo - margin <= t < hi + margin:   # half-open: the advanced
                t = hi + margin                  # t itself must not re-match
                changed = True
    return t


def event_rigid_clearance(tris_a: list[np.ndarray], tris_b: list[np.ndarray],
                          directions: list[np.ndarray],
                          neighbourhood: list[tuple[np.ndarray, np.ndarray]]
                          = (), margin: float = 1e-6) -> dict | None:
    """Best-found rigid translation of the A patch clearing all pairs.

    For each candidate direction: required displacement = max over crossing
    pairs of that pair's exit distance; then advanced past the union of
    intersection intervals with every (moving, static) neighbourhood pair.
    Returns the minimising direction, or None when no direction admits a
    finite clearance. The value is an upper bound on the unrestricted
    minimal displacement -- exact only per declared direction.
    """
    # batched interval evaluation: one vectorised call per direction over
    # all crossing pairs and all neighbourhood pairs (W4 item 1 -- the
    # scalar loop was thousands of tiny calls per direction)
    from time import perf_counter as _pc
    from .repair import KERNEL_PROFILE as _KP    # observation only
    cross_A = np.array([TA for TA in tris_a for _ in tris_b]) \
        if tris_a and tris_b else np.zeros((0, 3, 3))
    cross_B = np.array([TB for _ in tris_a for TB in tris_b]) \
        if tris_a and tris_b else np.zeros((0, 3, 3))
    nb_M = np.array([TM for TM, _ in neighbourhood]) \
        if len(neighbourhood) else np.zeros((0, 3, 3))
    nb_S = np.array([TS for _, TS in neighbourhood]) \
        if len(neighbourhood) else np.zeros((0, 3, 3))
    best = None
    for d in directions:
        d = np.asarray(d, float)
        n = np.linalg.norm(d)
        if n < 1e-12:
            continue
        d = d / n
        _t0 = _pc()
        lo, hi, ok = pair_intervals_batched(cross_A, cross_B, d)
        _KP["calls"]["pair_intervals_batched_crossing"] += 1
        _KP["seconds"]["pair_intervals_batched_crossing"] += _pc() - _t0
        inter = ok & (lo <= 0.0) & (hi >= 0.0)   # currently intersecting
        need = float(hi[inter].max()) if inter.any() else 0.0
        ahead = ok & (lo > 0.0)                  # partner ahead along d
        blocked = list(zip(lo[ahead].tolist(), hi[ahead].tolist()))
        _t0 = _pc()
        nlo, nhi, nok = pair_intervals_batched(nb_M, nb_S, d)
        _KP["calls"]["pair_intervals_batched_neighbourhood"] += 1
        _KP["seconds"]["pair_intervals_batched_neighbourhood"] += _pc() - _t0
        blocked += list(zip(nlo[nok].tolist(), nhi[nok].tolist()))
        t = clear_time(need + margin, blocked, margin)
        if np.isfinite(t) and (best is None or t < best["t"]):
            best = {"t": float(t), "direction": d.tolist()}
    return best


def pair_min_exit(A: np.ndarray, B: np.ndarray,
                  witness: bool = False):
    """Minimum-norm relative translation that separates two crossing
    triangles: the distance from the origin to the boundary of the convex
    Minkowski difference hull(B - A) (round 19/20's lower-bound primitive).

    Round-21 hardening: facet offsets are divided by their own normal
    norms explicitly (never trusting upstream normalization), and the
    minimizing facet is returned as a checkable witness on request.
    Returns 0.0 when the origin is not strictly inside (already separated
    or degenerate) -- conservative: never overstates.
    """
    from scipy.spatial import ConvexHull, QhullError

    pts = np.array([b - a for a in A for b in B])
    try:
        h = ConvexHull(pts)
    except QhullError:
        return (0.0, None) if witness else 0.0
    n = h.equations[:, :3]
    off = h.equations[:, 3]
    norms = np.linalg.norm(n, axis=1)
    good = norms > 1e-300
    dist = np.where(good, -off / np.maximum(norms, 1e-300), 0.0)
    if np.any(dist <= 1e-12):
        return (0.0, None) if witness else 0.0  # on/outside/degenerate
    i = int(np.argmin(dist))
    L = float(dist[i])
    if not witness:
        return L
    return L, {"facet_normal": (n[i] / norms[i]).tolist(),
               "facet_plane_distance": L}


def min_exit_from_pairs(tri_pairs, allowance_vx: float):
    """Certificate-grade event lower bound from the ACTUAL crossing
    triangle pairs (round 21: never the Cartesian product of regions).

    Returns (L_raw, L_safe, witness): L_raw = max over pairs of the
    per-pair minimum exit; L_safe = L_raw - allowance (the number compared
    to the budget); witness identifies the maximizing pair and facet.
    """
    L, wit = 0.0, None
    for i, (TA, TB) in enumerate(tri_pairs):
        li, wi = pair_min_exit(TA, TB, witness=True)
        if li > L:
            L, wit = li, {"pair_index": i, **(wi or {})}
    return L, max(0.0, L - allowance_vx), wit


def event_min_exit(tris_a, tris_b) -> float:
    """Lower bound on ANY rigid relative translation clearing the event:
    every crossing pair must be exited, so the norm must be at least the
    largest per-pair minimum exit. If this exceeds the admissible relative
    budget, NO rigid translation works -- a certificate of infeasibility,
    which direction sampling can never provide (round 20)."""
    L = 0.0
    for TA in tris_a:
        for TB in tris_b:
            iv = pair_interval(TA, TB, np.array([0.0, 0.0, 1.0]))
            if iv and iv[0] <= 0.0 <= iv[1]:    # currently intersecting
                L = max(L, pair_min_exit(TA, TB))
    return L
