"""How far apart, along the papyrus, are the two patches that meet?

A transverse self-intersection has two preimages in the mesh grid that occupy
one place in space. The distance between them measured ALONG the represented
surface is intrinsic: no revolution estimate, no axis, no centroid, no
neighbourhood radius. It is defined for a half-turn fragment and an
eighteen-turn trace alike.

Two review rounds shaped this file; their faults are pinned so they are not
reintroduced (notes/DECISIONS.md, rounds 11 and 12):

- **stride 2 was not topology-preserving** (round 11) -- it bridged holes and
  erased narrow connectors. Stride 1 only.
- **the graph must be the censused triangle complex, not the valid-vertex
  grid** (round 12). `selfcross` builds triangles only from quads whose four
  corners are all valid AND whose six pairwise corner distances are all within
  `maxedge`; a vertex chain that supports no triangle is not surface. Every
  edge here comes from a retained quad: its four sides and the chosen
  diagonal. The result is the exact shortest path on the triangulation's edge
  graph -- an upper approximation to the continuous polyhedral geodesic, since
  paths cannot cross triangle interiors.
- **event grouping must preserve branch identity** (round 12). Pooling every
  participating quad into one cluster produced self-paired and
  mixed-orientation events whose two "regions" contained each other -- ~26% of
  spectrum events, ~61% of the sub-0.1 mm mode. Grouping is now in product
  space with an orientation parity: two crossing pairs join one event only
  when their branches correspond side-to-side (allowing a consistent swap),
  and parity conflicts are flagged ambiguous, not measured.
- **endpoints are the intersection points themselves** (round 12). Corner-set
  endpoints carry an absolute error of up to two quad diameters, which is the
  entire size of the small-separation mode. Each pair's triangle-triangle
  intersection segment is recomputed and its endpoints enter the graph as
  virtual nodes tied barycentrically to their OWN triangle's corners, so the
  measured quantity is the walk along the sheet from the crossing locus back
  to the crossing locus via the other branch. Corner endpoints remain only as
  a flagged fallback for numerically degenerate pairs.

The per-pair scalar the census stores (`pen`) is the Euclidean length of the
triangle intersection segment, NOT a penetration depth; it is reported as
`intersection_length_vx` and must never be read as a displacement bound.
"""
from __future__ import annotations

import numpy as np

MAXEDGE_DEFAULT = 60.0          # selfcross default; census.json params confirm


def centreline_arclength(P: np.ndarray, V: np.ndarray) -> np.ndarray:
    """Cumulative distance along the surface's middle, per grid column.

    Kept as a cheap per-column coordinate (figure axes, candidate ranking).
    Not a substitute for the graph distance.
    """
    cols, have = [], []
    for u in range(P.shape[1]):
        m = V[:, u]
        if m.sum() >= 3:
            cols.append(np.median(P[m, u], axis=0))
            have.append(u)
    if len(have) < 2:
        return np.zeros(P.shape[1])
    C = np.asarray(cols)
    step = np.linalg.norm(np.diff(C, axis=0), axis=1)
    s_have = np.concatenate([[0.0], np.cumsum(step)])
    return np.interp(np.arange(P.shape[1]), np.asarray(have), s_have)


def retained_quads(P: np.ndarray, V: np.ndarray,
                   maxedge: float = MAXEDGE_DEFAULT) -> np.ndarray:
    """The quads selfcross keeps: four valid corners, six edges within maxedge.

    The six distances are the quad's four sides and both diagonals, exactly as
    engines/selfcross.cpp computes them before building triangles.
    """
    Q = V[:-1, :-1] & V[1:, :-1] & V[:-1, 1:] & V[1:, 1:]
    if maxedge and maxedge > 0:
        p00, p10 = P[:-1, :-1], P[1:, :-1]
        p01, p11 = P[:-1, 1:], P[1:, 1:]
        e = np.zeros(Q.shape)
        for a, b in ((p00, p01), (p01, p11), (p11, p10), (p00, p10),
                     (p00, p11), (p10, p01)):
            e = np.maximum(e, np.linalg.norm(a - b, axis=-1))
        Q &= e <= maxedge
    return Q


class SurfaceGraph:
    """The censused triangulation's edge graph, stride 1, full resolution.

    Vertices: grid points belonging to at least one retained quad. Edges: the
    four sides and the chosen diagonal of every retained quad -- `diagonal=0`
    joins (v,u)-(v+1,u+1) matching selfcross's d0 tessellation, `diagonal=1`
    joins (v+1,u)-(v,u+1). A quad selfcross dropped contributes nothing, so a
    path can only go where the censused surface exists.
    """

    def __init__(self, P: np.ndarray, V: np.ndarray, diagonal: int = 0,
                 maxedge: float = MAXEDGE_DEFAULT):
        from scipy import sparse
        from scipy.sparse.csgraph import connected_components

        nv, nu = V.shape
        Q = retained_quads(P, V, maxedge)
        self.Q = Q

        used = np.zeros((nv, nu), dtype=bool)
        used[:-1, :-1] |= Q
        used[1:, :-1] |= Q
        used[:-1, 1:] |= Q
        used[1:, 1:] |= Q
        idx = -np.ones((nv, nu), np.int64)
        vv, uu = np.nonzero(used)
        idx[vv, uu] = np.arange(len(vv))
        self.idx, self.n = idx, len(vv)
        self.P = P

        # An edge exists iff one of the (at most two) quads containing it is
        # retained; each edge is emitted exactly once (duplicates would sum in
        # the sparse build).
        rows, cols, data = [], [], []

        def add(a, b, pa, pb, mask):
            m = mask
            if m.any():
                rows.append(a[m]); cols.append(b[m])
                data.append(np.linalg.norm(pa[m] - pb[m], axis=1))

        H = np.zeros((nv, nu - 1), dtype=bool)      # (v,u)-(v,u+1)
        H[:-1, :] |= Q
        H[1:, :] |= Q
        add(idx[:, :-1], idx[:, 1:], P[:, :-1], P[:, 1:], H)
        W = np.zeros((nv - 1, nu), dtype=bool)      # (v,u)-(v+1,u)
        W[:, :-1] |= Q
        W[:, 1:] |= Q
        add(idx[:-1, :], idx[1:, :], P[:-1, :], P[1:, :], W)
        if diagonal == 0:
            add(idx[:-1, :-1], idx[1:, 1:], P[:-1, :-1], P[1:, 1:], Q)
        elif diagonal == 1:
            add(idx[1:, :-1], idx[:-1, 1:], P[1:, :-1], P[:-1, 1:], Q)
        # diagonal == -1: perimeter edges only -- the triangulation-neutral
        # graph used for deformation fields (round 16: the d0-vs-d1 field
        # weights differed by up to 0.23 on w094)

        if rows:
            r = np.concatenate(rows); c = np.concatenate(cols)
            w = np.concatenate(data)
        else:
            r = c = np.zeros(0, np.int64); w = np.zeros(0)
        self.g = sparse.coo_matrix(
            (np.concatenate([w, w]),
             (np.concatenate([r, c]), np.concatenate([c, r]))),
            shape=(self.n, self.n)).tocsr()
        self.ncomp, self.comp = connected_components(self.g, directed=False)
        self.diagonal = diagonal
        self.maxedge = maxedge

    # -- lookups ----------------------------------------------------------
    def vertex(self, v: int, u: int) -> int:
        return int(self.idx[v, u])

    def quad_corners(self, v: int, u: int) -> list[int]:
        nv, nu = self.idx.shape
        out = []
        for a, b in ((v, u), (v + 1, u), (v, u + 1), (v + 1, u + 1)):
            if 0 <= a < nv and 0 <= b < nu and self.idx[a, b] >= 0:
                out.append(int(self.idx[a, b]))
        return out

    def quad_triangles(self, v: int, u: int) -> list[np.ndarray]:
        """Vertex indices of the quad's two triangles, selfcross's order."""
        if self.diagonal not in (0, 1):
            raise ValueError("no triangulation on a perimeter-only graph")
        i = self.idx
        if self.diagonal == 0:
            t = [(i[v, u], i[v, u + 1], i[v + 1, u + 1]),
                 (i[v, u], i[v + 1, u + 1], i[v + 1, u])]
        else:
            t = [(i[v, u], i[v, u + 1], i[v + 1, u]),
                 (i[v, u + 1], i[v + 1, u + 1], i[v + 1, u])]
        return [np.array(x, np.int64) for x in t]

    def distance_field(self, sources: list[int]) -> np.ndarray:
        from scipy.sparse.csgraph import dijkstra
        return dijkstra(self.g, directed=False, indices=sources, min_only=True)

    def vertex_distance(self, a: tuple, b: tuple) -> float:
        """Exact graph distance between two grid vertices (test hook)."""
        d = self.distance_field([self.vertex(*a)])
        return float(d[self.vertex(*b)])

    def _super_source_field(self, offsets: dict[int, float]) -> np.ndarray:
        """Exact seeded field: distances from a virtual node tied to each
        seed vertex at its own offset weight. d[t] = min_a (w_a + d(a, t)).

        Zero offsets are clamped to a tiny epsilon so the sparse build cannot
        drop them as implicit zeros.
        """
        from scipy import sparse
        from scipy.sparse.csgraph import dijkstra

        n = self.n
        nodes = np.fromiter(offsets.keys(), np.int64, len(offsets))
        w = np.maximum(np.fromiter(offsets.values(), float, len(offsets)),
                       1e-12)
        col = sparse.coo_matrix((w, (nodes, np.zeros(len(nodes), np.int64))),
                                shape=(n, 1))
        aug = sparse.bmat([[self.g, col], [col.T, None]], format="csr")
        return dijkstra(aug, directed=False, indices=[n], min_only=True)[:n]

    def seeded_distance(self, seeds_a: list[tuple[int, float]],
                        seeds_b: list[tuple[int, float]],
                        max_pops: int | None = 2_000_000
                        ) -> tuple[float, bool]:
        """Min over virtual endpoints of (in-A weight + path + in-B weight).

        `seeds_a` / `seeds_b` are (vertex, offset) pairs: virtual points tied
        to graph vertices at a known extra distance (barycentric attachment).
        Exact early exit: Dijkstra runs from the A seeds; every settled vertex
        that carries a B offset updates the best answer, and the search stops
        once the heap minimum can no longer beat it.

        Returns (distance, method). If the pop budget runs out -- only possible
        on separations whose search ball exceeds ~max_pops vertices -- the
        computation restarts as a scipy full field from a virtual super-source
        tied to every A vertex at its own offset weight, which is exact for
        the seeded problem (round 13 rejected the earlier offset-mixing
        approximation). exact stays True on both paths; the flag remains for
        any future inexact strategy.
        """
        import heapq

        b_off: dict[int, float] = {}
        for n, w in seeds_b:
            b_off[n] = min(b_off.get(n, np.inf), w)
        comp_a = {int(self.comp[n]) for n, _ in seeds_a}
        if not (comp_a & {int(self.comp[n]) for n in b_off}):
            return float("inf"), "early_exit"

        indptr, indices, weights = self.g.indptr, self.g.indices, self.g.data
        dist: dict[int, float] = {}
        heap: list[tuple[float, int]] = []
        a_off: dict[int, float] = {}
        for n, w in seeds_a:
            a_off[n] = min(a_off.get(n, np.inf), w)
            if w < dist.get(n, np.inf):
                dist[n] = w
                heapq.heappush(heap, (w, n))
        best = np.inf
        pops = 0
        while heap:
            d, n = heapq.heappop(heap)
            if d >= best:
                break
            if d > dist.get(n, np.inf):
                continue
            if n in b_off:
                best = min(best, d + b_off[n])
            pops += 1
            if max_pops is not None and pops > max_pops:
                f = self._super_source_field(a_off)
                cand = min(float(f[t]) + wb for t, wb in b_off.items())
                return float(min(cand, best)), "full_super_source"
            lo, hi = indptr[n], indptr[n + 1]
            for m, w in zip(indices[lo:hi].tolist(), weights[lo:hi].tolist()):
                nd = d + w
                if nd < dist.get(m, np.inf):
                    dist[m] = nd
                    heapq.heappush(heap, (nd, m))
        return float(best), "early_exit"


# ---------------------------------------------------------------- grouping
def _adj(a: tuple, b: tuple) -> bool:
    return abs(a[0] - b[0]) <= 1 and abs(a[1] - b[1]) <= 1


def oriented_events(rec: np.ndarray) -> list[dict]:
    """Crossing pairs grouped into events WITH branch identity preserved.

    Two pairs belong to one event only when their branches correspond: side 1
    adjacent to side 1 and side 2 adjacent to side 2 (parity 0), or both
    swapped (parity 1). Parity is propagated through a union-find; a component
    whose constraints conflict is returned with `ambiguous=True` and must be
    reported, not measured. Region A/B are the parity-normalised quad sets.
    """
    n = len(rec)
    if n == 0:
        return []
    q1 = [(int(r["v1"]), int(r["u1"])) for r in rec]
    q2 = [(int(r["v2"]), int(r["u2"])) for r in rec]

    parent = list(range(n))
    parity = [0] * n            # orientation of row relative to its root
    conflict = [False] * n

    def find(x):
        chain = []
        while parent[x] != x:
            chain.append(x)
            x = parent[x]
        p = 0
        for c in reversed(chain):          # accumulate parity root-down
            p ^= parity[c]
            parent[c] = x
            parity[c] = p
        return x, (parity[chain[0]] if chain else 0)

    def union(a, b, rel):
        ra, pa = find(a)
        rb, pb = find(b)
        if ra == rb:
            if pa ^ pb != rel:
                conflict[ra] = True
            return
        parent[rb] = ra
        parity[rb] = pa ^ pb ^ rel
        conflict[ra] |= conflict[rb]

    cell: dict = {}
    for i in range(n):
        for q in (q1[i], q2[i]):
            cell.setdefault(q, []).append(i)
    for i in range(n):
        seen = set()
        for q in (q1[i], q2[i]):
            for dv in (-1, 0, 1):
                for du in (-1, 0, 1):
                    for j in cell.get((q[0] + dv, q[1] + du), ()):
                        if j <= i or j in seen:
                            continue
                        seen.add(j)
                        straight = _adj(q1[i], q1[j]) and _adj(q2[i], q2[j])
                        swapped = _adj(q1[i], q2[j]) and _adj(q2[i], q1[j])
                        if straight and swapped:
                            union(i, j, 0)
                            union(i, j, 1)   # forces the conflict flag
                        elif straight:
                            union(i, j, 0)
                        elif swapped:
                            union(i, j, 1)

    groups: dict = {}
    for i in range(n):
        root, p = find(i)
        groups.setdefault(root, []).append((i, p))
    out = []
    for root, members in groups.items():
        A, B, idxs = set(), set(), []
        for i, p in members:
            idxs.append(i)
            (A if p == 0 else B).add(q1[i])
            (B if p == 0 else A).add(q2[i])
        touching = bool(A & B) or any(
            (a[0] + dv, a[1] + du) in B
            for a in A for dv in (-1, 0, 1) for du in (-1, 0, 1))
        out.append({"rows": idxs,
                    "flip": [p for _, p in members],
                    "region_a": A, "region_b": B,
                    "ambiguous": bool(conflict[root]),
                    "self_touching": touching})
    return out


# ------------------------------------------------- barycentric endpoints
def _tri_tri_segment(A: np.ndarray, B: np.ndarray,
                     eps: float = 1e-12) -> tuple | None:
    """Endpoints of the intersection segment of two triangles, or None.

    Interval method: each triangle is clipped by the other's plane; both
    clipped segments lie on the planes' intersection line; the overlap of
    their parameter intervals is the shared segment.
    """
    def clip(T, n, d0):
        s = T @ n + d0
        pts = []
        for i in range(3):
            j = (i + 1) % 3
            si, sj = s[i], s[j]
            if si == 0.0:
                pts.append(T[i])
            if (si > 0) != (sj > 0) and si != 0.0 and sj != 0.0:
                t = si / (si - sj)
                pts.append(T[i] + t * (T[j] - T[i]))
        return pts

    n1 = np.cross(A[1] - A[0], A[2] - A[0])
    n2 = np.cross(B[1] - B[0], B[2] - B[0])
    d = np.cross(n1, n2)
    L = np.linalg.norm(d)
    if L < eps * max(np.linalg.norm(n1), 1.0) * max(np.linalg.norm(n2), 1.0):
        return None                              # near-coplanar: fallback
    d = d / L
    segA = clip(A, n2, -float(n2 @ B[0]))
    segB = clip(B, n1, -float(n1 @ A[0]))
    if len(segA) < 2 or len(segB) < 2:
        return None
    ta = sorted((float(p @ d), tuple(p)) for p in segA[:2])
    tb = sorted((float(p @ d), tuple(p)) for p in segB[:2])
    lo, hi = max(ta[0][0], tb[0][0]), min(ta[1][0], tb[1][0])
    if hi <= lo:
        return None
    a0, a1 = np.array(ta[0][1]), np.array(ta[1][1])
    span = ta[1][0] - ta[0][0]
    if span < eps:
        return None
    p_lo = a0 + (lo - ta[0][0]) / span * (a1 - a0)
    p_hi = a0 + (hi - ta[0][0]) / span * (a1 - a0)
    return p_lo, p_hi


def _pair_seeds(g: SurfaceGraph, X: np.ndarray, v1, u1, v2, u2):
    """Barycentric seeds for one crossing pair: ((A seeds, B seeds), exact).

    Recomputes which triangle combination of the two quads intersects and
    attaches each intersection-segment endpoint to its own triangle's three
    corners at straight-line in-plane distance (a triangle is flat, so that is
    the exact within-triangle geodesic). Falls back to zero-offset corner
    seeds when no combination yields a clean transverse segment.
    """
    combos = []
    for t1 in g.quad_triangles(v1, u1):
        for t2 in g.quad_triangles(v2, u2):
            if (t1 < 0).any() or (t2 < 0).any():
                continue
            seg = _tri_tri_segment(X[t1], X[t2])
            if seg is not None:
                combos.append((t1, t2, seg))
    if not combos:
        a = [(c, 0.0) for c in g.quad_corners(v1, u1)]
        b = [(c, 0.0) for c in g.quad_corners(v2, u2)]
        return (a, b), False
    sa, sb = [], []
    for t1, t2, (p, q) in combos:
        for x in (p, q):
            sa += [(int(c), float(np.linalg.norm(X[c] - x))) for c in t1]
            sb += [(int(c), float(np.linalg.norm(X[c] - x))) for c in t2]
    return (sa, sb), True


def event_separation(g: SurfaceGraph, rec: np.ndarray, ev: dict,
                     X: np.ndarray) -> dict:
    """Surface distance from the event's crossing locus back to itself.

    Sources: the intersection points attached to their region-A triangles;
    targets: the same points attached to their region-B triangles. The
    measured number is the shortest walk along the censused surface from the
    locus via one branch to the locus via the other -- "how far the trace
    travelled before landing back here", with no corner-set slack.
    """
    seeds_a, seeds_b = [], []
    exact = True
    for i, p in zip(ev["rows"], ev["flip"]):
        r = rec[i]
        (sa, sb), ok = _pair_seeds(g, X, int(r["v1"]), int(r["u1"]),
                                   int(r["v2"]), int(r["u2"]))
        exact &= ok
        if p == 1:
            sa, sb = sb, sa
        seeds_a += sa
        seeds_b += sb
    if not seeds_a or not seeds_b:
        return {"separation_vx": None, "same_component": None,
                "endpoint_exact": False}
    same = len({int(g.comp[n]) for n, _ in seeds_a}
               | {int(g.comp[n]) for n, _ in seeds_b}) == 1
    sep, method = g.seeded_distance(seeds_a, seeds_b)
    return {"separation_vx": (sep if np.isfinite(sep) else None),
            "same_component": bool(same and np.isfinite(sep)),
            "endpoint_exact": bool(exact),      # tri-tri reconstruction held
            "distance_exact": True,              # both methods are exact
            # operational provenance, distinct from exactness (round 14)
            "distance_method": method}


def segment_spectrum(P: np.ndarray, V: np.ndarray, rec: np.ndarray,
                     voxel_um: float, diagonal: int = 0,
                     max_events: int | None = None,
                     maxedge: float = MAXEDGE_DEFAULT) -> list[dict]:
    """Per-event intrinsic severities for one segment: the spectrum.

    Ambiguous events (orientation conflict) are reported with flags and no
    separation. Events are processed largest-first (by pair count) so a cap
    keeps the ones that carry the most surface; the cap is the caller's to
    report, never silent.
    """
    if len(rec) == 0:
        return []
    g = SurfaceGraph(P, V, diagonal, maxedge)
    vv, uu = np.nonzero(g.idx >= 0)
    X = np.empty((g.n, 3))
    X[g.idx[vv, uu]] = P[vv, uu]

    evs = sorted(oriented_events(rec), key=lambda e: -len(e["rows"]))
    if max_events is not None and len(evs) > max_events:
        evs = evs[:max_events]
    mm = voxel_um / 1000.0
    out = []
    for ev in evs:
        arr = rec[ev["rows"]]
        du = np.abs(arr["u1"].astype(int) - arr["u2"].astype(int))
        row = {
            "n_pairs": int(len(arr)),
            "ambiguous": ev["ambiguous"],
            "self_touching": ev["self_touching"],
            "du_max": int(du.max()),
            "median_intersection_length_vx": float(np.nanmedian(arr["pen"])),
        }
        if ev["ambiguous"]:
            row.update({"separation_mm": None, "same_component": None,
                        "endpoint_exact": None, "distance_exact": None})
        else:
            r = event_separation(g, rec, ev, X)
            row.update({
                "separation_mm": (round(r["separation_vx"] * mm, 4)
                                  if r["separation_vx"] is not None else None),
                "same_component": r["same_component"],
                "endpoint_exact": r["endpoint_exact"],
                "distance_exact": r["distance_exact"],
            })
        out.append(row)
    return out
