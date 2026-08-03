"""Do two surfaces pass through each other?

A per-surface verdict says nothing about assembly. Two surfaces that are
each free of self-intersection can still interpenetrate, and any merge of
such a pair either self-intersects or must discard part of one. Mergers
that build a graph over patches and join along its edges need to know
which proposed edges are geometrically impossible BEFORE assembling.

There are two ways to put a pair in front of the census, and both are
here. The engine can load several surfaces from one atlas, scope its
adjacency exclusion to within a surface, and tag every contact with the
two surface ids it came from; that is `classify_many`, and it answers a
whole edge list in two engine passes per batch. Or both surfaces can be
written into a single grid separated by a band of invalid rows far wider
than the adjacency exclusion; that is `classify`, one pair at a time.

The second exists because it came first and because it shares nothing
with the engine's surface bookkeeping, so it is a genuine second opinion
rather than a rephrasing of the first. `test_check_pairs_batch_matches_
stitch` holds them against each other.

Either way each surface's own contacts are identified as its own and
never blamed on the pair: a patch that self-intersects is a separate
finding, reported in `self_a` / `self_b`.

The clean verdict is `no_transverse_conflict`, deliberately not
"compatible". This test is silent about parametrisation, seams, scale
agreement and every failure mode that is not interpenetration.

It is also silent about coordinate frames. Two surfaces are compared as
sets of numbers, and nothing here establishes that those numbers refer to
the same volume at the same voxel scale. Surfaces traced from different
scans could be numerically coincident and reported as a conflict they do
not physically have. A candidate edge list from one workflow satisfies
this by construction; one assembled by hand does not.
"""
from __future__ import annotations

import shutil
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import tifffile

from . import tifxyz

# Comfortably wider than the exclusion of 1, so the two surfaces can never
# be adjacent in the grid however they are shaped.
GAP_ROWS = 8

NO_CONFLICT = "no_transverse_conflict"
CONFLICT = "transverse_conflict"
NOT_TESTABLE = "not_testable"

VERDICT_NOTE = {
    NO_CONFLICT: ("the two surfaces do not pass through each other; this is "
                  "not a statement that they are compatible to merge"),
    CONFLICT: ("the two surfaces interpenetrate, so any merge that keeps both "
               "in full contains a self-intersection"),
    NOT_TESTABLE: ("too little of either surface lies within the other's "
                   "extent to decide; reported rather than silently dropped"),
}

# Batch TARGETS, not hard caps: both endpoints of an edge always share an
# atlas, so a single edge may exceed either. Triangles rather than grid
# cells, because a grid cell need not produce one -- it is the triangle
# count that sets the engine's memory, at roughly 88 bytes each plus its
# broad-phase bucket membership, so 12M triangles is on the order of 2 GB.
# A hundred patches are small; three full traces are not.
MAX_BATCH_SURFACES = 64
MAX_BATCH_TRIANGLES = 12_000_000


@dataclass
class PairResult:
    a: str
    b: str
    verdict: str
    transverse_d0: int = 0
    transverse_d1: int = 0
    transverse_both: int = 0
    self_a: int = 0
    self_b: int = 0
    n_valid_a: int = 0
    n_valid_b: int = 0
    #: Contact counts are triangle-pair rows summed over BOTH triangulations,
    #: not distinct crossing events: one crossing region yields many rows, and
    #: a quad contributes under each diagonal. They rank severity; they do not
    #: count places.
    overlap_points: int = 0
    angles_deg: dict = field(default_factory=dict)
    max_penetration_vx: float = 0.0
    reason: str = ""

    def as_dict(self) -> dict:
        d = dict(self.__dict__)
        d["verdict_means"] = VERDICT_NOTE[self.verdict]
        return d


# --------------------------------------------------------------- contacts

WANT = "transverse"


def _read_contacts(csv: Path, want: str = WANT
                   ) -> Iterator[tuple[int, int, int, int, float, float]]:
    """Yield `(v1, v2, surf1, surf2, penetration, angle)` for one verdict.

    The engine writes `surf1`/`surf2` only when the atlas held more than one
    surface, so the header decides whether they are read or defaulted to
    zero. Anything else is a schema change and must fail loudly rather than
    be papered over with `.get`.

    A census CSV is mostly grazing contacts -- a real edge list runs to
    hundreds of thousands of rows, of which a handful matter -- so rows are
    rejected by a substring test before anything is split or converted. No
    verdict name contains another, so the test cannot pass a row it should
    have dropped, and the columns are located by name rather than position.
    """
    with Path(csv).open() as fh:
        header = fh.readline().strip().split(",")
        col = {name: i for i, name in enumerate(header)}
        if not {"v1", "v2", "verdict", "penetration", "angle_deg"} <= col.keys():
            raise ValueError(f"unexpected census schema in {csv}: {header}")
        iv1, iv2 = col["v1"], col["v2"]
        ipen, iang, ivd = col["penetration"], col["angle_deg"], col["verdict"]
        is1 = col.get("surf1", -1)
        is2 = col.get("surf2", -1)
        if (is1 < 0) != (is2 < 0):
            raise ValueError(f"census schema in {csv} has one surface column "
                             f"but not the other: {header}")
        for line in fh:
            if want not in line:
                continue
            p = line.rstrip("\n").split(",")
            if p[ivd] != want:
                continue
            yield (int(p[iv1]), int(p[iv2]),
                   int(p[is1]) if is1 >= 0 else 0,
                   int(p[is2]) if is2 >= 0 else 0,
                   float(p[ipen]), float(p[iang]))


def _accumulate(csv: Path, diag: int, key_of: Callable, stats: dict) -> None:
    """Bin one diagonal's transverse contacts by which surfaces they join."""
    for v1, v2, s1, s2, pen, ang in _read_contacts(csv):
        i, j = key_of(v1, v2, s1, s2)
        e = stats.setdefault((i, j) if i <= j else (j, i),
                             {"d0": 0, "d1": 0, "total": 0,
                              "_ang": [], "_pen": []})
        e[f"d{diag}"] += 1
        e["total"] += 1
        e["_ang"].append(ang)
        e["_pen"].append(pen)


def _finalise(stats: dict) -> dict:
    """Turn the accumulated angle and penetration lists into summaries."""
    for e in stats.values():
        angs = np.asarray(e.pop("_ang"), float)
        pens = np.asarray(e.pop("_pen"), float)
        if angs.size:
            e["angles"] = {
                "median": float(np.median(angs)),
                "p10": float(np.percentile(angs, 10)),
                "p90": float(np.percentile(angs, 90)),
                "fraction_below_10_deg": float(np.mean(angs < 10))}
            e["max_penetration"] = float(pens.max())
    return stats


# ------------------------------------------------------- the batched path

def census_batch(paths: Sequence[Path], work: Path, tag: str, *,
                 surfaces: Sequence | None = None,
                 params: dict | None = None) -> dict:
    """Census a whole set of surfaces at once. Two engine passes, not 2N.

    Returns a dict keyed by ordered surface-index pairs: `(i, i)` is what
    surface `i` does to itself, `(i, j)` with `i < j` is what surfaces `i`
    and `j` do to each other. Only TRANSVERSE contacts are counted.

    Which surfaces share the atlas does not change any transverse count.
    The engine's broad phase owns each triangle pair at the cell containing
    the minimum corner of their overlapping bounding boxes; two triangles
    that genuinely cross have overlapping boxes, that corner lies inside
    both, and cell ranges are inclusive, so both are in that bucket and the
    pair is tested there and rejected everywhere else. Adding surfaces moves
    the grid origin, which renumbers the owner cell without removing it.
    That is what makes it legitimate to read a surface's own census out of a
    run it shares with others.

    The guarantee is for transverse contacts only, and deliberately not for
    the tolerance-classified ones. Two nearly coincident sheets can be
    separated in space yet within the grazing tolerance; whether they land
    in a shared bucket then depends on where the cell boundaries fall, and
    the broad phase does not inflate boxes by that tolerance. Counting
    coplanar or grazing contacts here would need it to.
    """
    from . import pipeline
    n = len(paths)
    stats: dict = {}
    for d in (0, 1):
        csvp, _ = pipeline.run_engine_multi(list(paths), f"{tag}_b", work, d,
                                            params, surfaces=surfaces)
        _accumulate(csvp, d, lambda v1, v2, s1, s2: (s1, s2), stats)
    # An id outside the atlas would silently become a pair nobody asked
    # about, or worse be read as a different edge's verdict.
    bad = [k for k in stats if not (0 <= k[0] < n and 0 <= k[1] < n)]
    if bad:
        raise ValueError(f"census reported surface ids outside the atlas of "
                         f"{n}: {sorted(bad)[:5]}")
    return _finalise(stats)


def _batches(edges: Sequence[tuple[Path, Path]], max_surfaces: int,
             max_triangles: int, maxedge: float | None) -> Iterator[tuple]:
    """Group edges so both endpoints of every edge share one atlas.

    Greedy and order-preserving: an edge joins the open batch unless its new
    surfaces would push it past either bound, in which case the batch is
    closed first. A surface appearing in several batches is simply censused
    in each of them, which costs a repeat of its self-census and nothing
    else. Each batch yields `(paths, index, surfaces, edges)`, where the
    position of a path in `paths` is the surface id the engine will give it.

    Surfaces are identified by their resolved path, so one surface reached
    by two spellings -- a symlink, a relative and an absolute form -- is one
    atlas entry rather than two. Loaded twice it would be its own perfect
    duplicate, and every contact between the copies would be coplanar rather
    than transverse, so the verdict would come back clean for a pair that was
    never really tested. Edges keep the caller's spelling for reporting.
    """
    order: list[Path] = []
    index: dict[Path, int] = {}
    cache: dict[Path, object] = {}
    pending: list[tuple[int, Path, Path, Path, Path]] = []
    tris = 0
    sized: dict[Path, int] = {}

    def close():
        return (list(order), dict(index), dict(cache), list(pending),
                {p: sized[p] for p in order})

    for ei, (a, b) in enumerate(edges):
        ka, kb = Path(a).resolve(), Path(b).resolve()
        # Read before deciding to flush: the surfaces this edge needs are
        # wanted either way, and the flush is what drops the old cache.
        need = list(dict.fromkeys((ka, kb)))
        fresh = {p: tifxyz.read(p) for p in need if p not in index}
        for p, s in fresh.items():
            if p not in sized:
                sized[p] = n_triangles(s, maxedge)
        add = sum(sized[p] for p in fresh)
        if order and (len(order) + len(fresh) > max_surfaces
                      or tris + add > max_triangles):
            yield close()
            order, index, cache, pending, tris = [], {}, {}, [], 0
            # The flush took the cache with it, so an endpoint that was
            # already in the closed batch is not carried over -- it has to
            # be read again for the new one. Both endpoints of an edge share
            # an atlas unconditionally: an edge is atomic, and a batch that
            # holds only one may exceed either bound.
            for p in need:
                if p not in fresh:
                    fresh[p] = tifxyz.read(p)
        for p, s in fresh.items():
            index[p] = len(order)
            order.append(p)
            cache[p] = s
            tris += sized[p]
        pending.append((ei, a, b, ka, kb))

    if pending:
        yield close()


def classify_many(edges: Sequence[tuple[Path, Path]], work: Path, *,
                  tag: str = "batch", min_overlap_points: int = 8,
                  max_batch_surfaces: int = MAX_BATCH_SURFACES,
                  max_batch_triangles: int = MAX_BATCH_TRIANGLES,
                  params: dict | None = None,
                  on_result: Callable | None = None) -> list[PairResult]:
    """Classify a whole candidate edge list. Results follow the input order.

    Two engine passes per batch replace six per edge, and a surface that
    appears in many edges is read and censused once for its whole batch
    rather than once per edge it touches.
    """
    results: list[PairResult | None] = [None] * len(edges)

    maxedge = (params or {}).get("maxedge")
    for nb, (paths, index, cache, pending, sized) in enumerate(
            _batches(edges, max_batch_surfaces, max_batch_triangles,
                     maxedge), 1):
        stats = census_batch(paths, work, f"{tag}{nb:04d}",
                             surfaces=[cache[p] for p in paths], params=params)
        for ei, a, b, ka, kb in pending:
            r = _pair_result(a, b, cache[ka], cache[kb], stats,
                             index[ka], index[kb], min_overlap_points,
                             sized[ka], sized[kb])
            results[ei] = r
            if on_result is not None:
                on_result(ei, r)

    missing = [i for i, r in enumerate(results) if r is None]
    if missing:
        raise RuntimeError(f"{len(missing)} edges left unclassified: {missing[:5]}")
    return results  # type: ignore[return-value]


# -------------------------------------------------------- the stitched path

def stitch(a: Path, b: Path, dst: Path, surfaces: Sequence | None = None
           ) -> tuple[int, int]:
    """Write both surfaces into one grid, separated by invalid rows."""
    sa, sb = surfaces if surfaces is not None else (tifxyz.read(a),
                                                    tifxyz.read(b))
    va, ua = sa.shape
    vb, ub = sb.shape
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    V, U = va + GAP_ROWS + vb, max(ua, ub)
    keep = np.zeros((V, U), bool)
    keep[:va, :ua] = sa.valid
    keep[va + GAP_ROWS:, :ub] = sb.valid
    for i, ax in enumerate("xyz"):
        out = np.full((V, U), -1.0, np.float32)
        pa = np.asarray(sa.points[..., i], np.float32).copy()
        pa[~sa.valid] = -1.0
        out[:va, :ua] = pa
        pb = np.asarray(sb.points[..., i], np.float32).copy()
        pb[~sb.valid] = -1.0
        out[va + GAP_ROWS:, :ub] = pb
        tifffile.imwrite(dst / f"{ax}.tif", out)
    tifxyz.write_meta(a, dst, np.stack(
        [np.asarray(tifffile.imread(dst / f"{ax}.tif")) for ax in "xyz"], -1),
        keep)
    return sa.n_valid, sb.n_valid


def classify(a: Path, b: Path, work: Path, tag: str, *,
             min_overlap_points: int = 8,
             params: dict | None = None) -> PairResult:
    """Classify one candidate edge on its own, via the stitched grid.

    `classify_many` is the fast path and the one the CLI uses. This one
    reaches the same verdict without the engine ever being told that more
    than one surface is present: the two are laid out in a single grid and
    a contact belongs to the pair exactly when its two ends fall on
    opposite sides of the invalid band. Keeping it exercised is what makes
    the equivalence test worth anything.
    """
    sa, sb = tifxyz.read(a), tifxyz.read(b)
    res = PairResult(a=a.name, b=b.name, verdict=NOT_TESTABLE,
                     n_valid_a=sa.n_valid, n_valid_b=sb.n_valid)
    maxedge = (params or {}).get("maxedge")
    if _untestable(res, sa, sb, n_triangles(sa, maxedge),
                   n_triangles(sb, maxedge)):
        return res

    dst = work / f"pair_{tag}.tifxyz"
    stitch(a, b, dst, surfaces=(sa, sb))
    split = sa.shape[0]
    try:
        from . import pipeline
        stats: dict = {}
        for d in (0, 1):
            csvp, _ = pipeline.run_engine(dst, f"{tag}_ab", work, d, params)
            _accumulate(csvp, d,
                        lambda v1, v2, s1, s2: (int(v1 >= split),
                                                int(v2 >= split)), stats)
        _finalise(stats)
    finally:
        shutil.rmtree(dst, ignore_errors=True)

    _fill(res, stats, 0, 1)
    if _sparse_overlap(res, min_overlap_points):
        res.verdict = NOT_TESTABLE
    return res


# ------------------------------------------------------------------ shared

def extents_meet(a: Path, b: Path) -> tuple[bool, int]:
    """Do the two surfaces' extents intersect, and by how many points?"""
    return _extents_meet(tifxyz.read(a), tifxyz.read(b))


def _extents_meet(sa, sb) -> tuple[bool, int]:
    """Testability is decided by whether the extents meet at all.

    NOT by a threshold on points-inside-the-other's-box. A nearly planar
    surface has a box that is degenerate in one axis, so counting points
    inside it under-reports badly and would call a genuinely testable pair
    undecidable. The census is the authority on contact; this is only a
    cheap guard against censusing two surfaces that cannot touch.
    """
    PA, PB = sa.points[sa.valid], sb.points[sb.valid]
    if len(PA) < 4 or len(PB) < 4:
        return False, 0
    lo = np.maximum(PA.min(axis=0), PB.min(axis=0))
    hi = np.minimum(PA.max(axis=0), PB.max(axis=0))
    if np.any(hi < lo):
        return False, 0
    inside = int(np.all((PA >= lo) & (PA <= hi), axis=1).sum()
                 + np.all((PB >= lo) & (PB <= hi), axis=1).sum())
    return True, inside


def n_triangles(s, maxedge: float | None = None) -> int:
    """How many triangles this surface contributes to a census.

    A surface can be entirely valid and contribute nothing: a checkerboard
    of valid cells contains no complete quad, and a coarse grid can have
    every quad dropped for edge length. Such a surface is not clean -- there
    was never anything to test -- so a pair containing one has to be refused
    rather than reported as free of contact.
    """
    from .intrinsic import retained_quads
    from .pipeline import CENSUS
    e = CENSUS["maxedge"] if maxedge is None else maxedge
    return 2 * int(retained_quads(s.points, s.valid, e).sum())


def _same_geometry(sa, sb) -> bool:
    """Are these two surfaces the same geometry under different names?

    A duplicate is not a pair. Censused as two atlas entries it would be its
    own perfect copy: contacts between the copies are coplanar, which this
    tool does not count, and any self-crossing the surface has appears
    doubled as a cross-surface transverse contact. Neither reading is a
    statement about two surfaces.
    """
    if sa.points.shape != sb.points.shape:
        return False
    if not np.array_equal(sa.valid, sb.valid):
        return False
    return bool(np.array_equal(sa.points[sa.valid], sb.points[sb.valid]))


def _untestable(res: PairResult, sa, sb, n_tri_a: int, n_tri_b: int) -> bool:
    """Refuse a pair the census cannot say anything about either way.

    These are refusals of the QUESTION, not of the answer: duplicate
    geometry, a surface with no triangles, extents that cannot meet. None
    of them can hide a contact, so unlike the sparse-overlap gate below
    they are safe to apply before the census runs.
    """
    if _same_geometry(sa, sb):
        res.reason = ("the two surfaces are the same geometry under different "
                      "names, so there is no pair to test")
        return True
    for n, name in ((n_tri_a, res.a), (n_tri_b, res.b)):
        if n == 0:
            res.reason = (f"{name} contributes no triangle to a census -- no "
                          f"four valid corners within the edge limit -- so "
                          f"there is no surface to test it against")
            return True
    meet, res.overlap_points = _extents_meet(sa, sb)
    if not meet:
        res.reason = ("the two surfaces' extents do not intersect, so they "
                      "cannot touch")
        return True
    return False


def _sparse_overlap(res: PairResult, min_overlap_points: int) -> bool:
    """Was there too little shared region to call a clean result decided?

    This one runs AFTER the census, never before. It counts vertices inside
    the shared bounding box, and two surfaces can cross through each other's
    interiors without either contributing a single vertex to that box --
    orthogonal planes are the simple case. Applied first it discards real
    conflicts the engine has already found. A detected contact outranks
    every applicability heuristic, including this one.
    """
    if res.transverse_both > 0 or res.overlap_points >= min_overlap_points:
        return False
    res.reason = (f"extents intersect but only {res.overlap_points} points "
                  f"of either surface lie in the shared region, below the "
                  f"{min_overlap_points} needed to call the absence of "
                  f"contact decided")
    return True


def _fill(res: PairResult, stats: dict, ia: int, ib: int) -> PairResult:
    """Read one pair's verdict out of a census keyed by surface index.

    The pair's count is read directly off the contacts that join the two,
    not obtained by subtracting self-censuses from a joint total. The
    subtraction was arithmetic that could go negative and had to be clamped;
    this cannot, and it means `self_a` and `self_b` are reported rather than
    load-bearing.
    """
    res.self_a = stats.get((ia, ia), {}).get("total", 0)
    res.self_b = stats.get((ib, ib), {}).get("total", 0)
    e = stats.get((ia, ib) if ia <= ib else (ib, ia), {})
    res.transverse_d0 = e.get("d0", 0)
    res.transverse_d1 = e.get("d1", 0)
    res.transverse_both = e.get("total", 0)
    res.angles_deg = e.get("angles", {})
    res.max_penetration_vx = e.get("max_penetration", 0.0)
    res.verdict = CONFLICT if res.transverse_both > 0 else NO_CONFLICT
    if res.verdict == NO_CONFLICT and (res.self_a or res.self_b):
        res.reason = ("no contact between the two; one or both self-intersect "
                      "on their own, which is reported separately")
    return res


def _pair_result(a: Path, b: Path, sa, sb, stats: dict, ia: int, ib: int,
                 min_overlap_points: int, n_tri_a: int, n_tri_b: int
                 ) -> PairResult:
    res = PairResult(a=a.name, b=b.name, verdict=NOT_TESTABLE,
                     n_valid_a=sa.n_valid, n_valid_b=sb.n_valid)
    if ia == ib:
        res.reason = ("both endpoints resolve to the same surface, so there "
                      "is no pair to test")
        return res
    if _untestable(res, sa, sb, n_tri_a, n_tri_b):
        return res
    _fill(res, stats, ia, ib)
    if _sparse_overlap(res, min_overlap_points):
        res.verdict = NOT_TESTABLE
    return res
