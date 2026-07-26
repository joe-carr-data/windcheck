"""Self-gap analysis of a traced surface, with a validity filter.

For every point on a multi-wrap trace, measure the distance to the nearest other
part of the SAME trace, excluding anything within `exclude_u` grid columns of the
point's own u. What remains is the trace's neighbouring wrap.

A correct trace reads one sheet spacing there (~17 voxels on PHerc0172). Near
zero means the trace came back onto a wrap it had already traced.

The threshold that matters is not ours. `GrowSurface.cpp:121` sets
`same_surface_th = 2.0`, and `GrowSurface.cpp:2919-2931` rejects any growth
candidate landing within 2.0 voxels of already-traced surface -- their own debug
string calls it "candidate rejected". So `frac_below_grower_th` counts points in
the published output that sit inside the pipeline's own rejection criterion.

Caveat that must travel with any number from here: that guard runs at growth
time on a coarse grid via a stochastic search, while the published tifxyz is
post-optimisation (`GrowSurface.cpp:1362` rewrites every point with no re-check
and no injectivity test). A violation in the output does not prove the guard
failed. What it does show is that nothing re-checks the criterion afterwards.

VALIDITY FILTER. Exclusion is by u-index, which assumes u advances with angle.
If a trace's flagged partners sit at |du| barely above `exclude_u`, they may be
same-wrap neighbours rather than the next wrap, and the result is an artifact.
A trace is only reported when its flagged |du| sits comfortably clear of the
cutoff. This is the one hole found in the method, so it is enforced in code
rather than eyeballed.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage

from . import atlas, tifxyz

GROWER_TH = 2.0        # GrowSurface.cpp:121  same_surface_th
FLAG_TH = 6.0          # our wider flag, for spatial clustering
DU_SAFETY = 2.0        # flagged |du| must exceed DU_SAFETY * exclude_u
AUTO_FRACTION = 0.25   # adaptive cutoff: this fraction of the measured period

# DENOMINATOR. Every reported fraction divides by the number of valid grid
# cells *submitted*, never by the number that came back with a finite distance.
#
# The published version divided by the finite ones. On a trace where almost
# nothing was measurable that denominator is tiny and the percentage stops being
# a rate: w081 published a blob_fraction of 6.67% that was 10 cells out of 150,
# ranking it second in the whole corpus behind a trace with 13,807 out of
# 196,802. Its honest unconditional figure is 0.036%.
#
# Censoring is not missing data here. A query whose nearest admissible surface
# lies beyond max_dist is a *known negative* for "nearest surface under 6 vx" --
# it cannot secretly be under 6. So the unconditional rate is well defined for
# every submitted point, and no coverage threshold is needed to compare traces.
#
# A COVERAGE_FLOOR was tried first and removed. It only hid the bad denominator,
# and it was not the untuned cut it claimed to be: the bimodal coverage gap is a
# property of the 256 vx search radius, not of the scroll. At 1024 vx the low
# group climbs from 2-4% to 24-34% and the gap fills in. Report `coverage` as an
# applicability statistic and let the reader see it; do not gate on it.


def estimate_exclude_u(path: Path, stride: int, work: Path, threads: int,
                       probe: int = 12) -> int | None:
    """Choose an exclusion width from the trace's own revolution period.

    A fixed cutoff in grid columns is not portable. The same surface published
    as tifxyz spans ~3065 columns and as OBJ ~667 parameter units, so a cutoff
    of 60 is 2% of a revolution in one and 55% in the other -- and the second
    gets rejected by the validity filter for being too close to its own cutoff.

    So: probe with a deliberately small exclusion, read the revolution period
    off where the flagged partners actually sit, and set the real cutoff to a
    quarter of it. That is small enough to keep the neighbouring wrap in view
    and large enough to exclude the local surface, in whatever units the file
    happens to use.
    """
    from . import objmesh

    if Path(path).suffix.lower() == ".obj":
        mesh = objmesh.read(path)
        if not mesh.has_param:
            return None
        pts, tags = objmesh.sample_points(mesh, stride=stride)
        span = int(mesh.tri_u.max() - mesh.tri_u.min())
        atlas.write_atlas_mesh(mesh, work / "pr_atlas.bin")
    else:
        surf = tifxyz.read(path)
        v, u = np.nonzero(surf.valid[::stride, ::stride])
        pts, tags = surf.points[v * stride, u * stride], u * stride
        span = int(surf.shape[1])
        atlas.write_atlas([_Entry(path)], work / "pr_atlas.bin")

    ex = max(span // 100, 3)
    atlas.write_queries_grouped(pts, tags, work / "pr_query.bin")
    r = atlas.run_engine(work / "pr_atlas.bin", work / "pr_query.bin",
                         work / "pr_result.bin", threads=threads, exclude_u=ex)
    d, w1 = r["d1"], r["w1"]
    m = np.isfinite(d) & (d < FLAG_TH * 4)
    if m.sum() < 50:
        return None
    period = float(np.median(np.abs(w1[m] - np.asarray(tags)[m])))
    if not np.isfinite(period) or period < 4 * ex:
        return None
    return int(max(period * AUTO_FRACTION, 3))


@dataclass
class Result:
    name: str
    u_extent: int
    n_queries: int                # valid grid cells submitted -- THE denominator
    n_measurable: int             # of those, the ones with a finite self-gap
    coverage: float               # n_measurable / n_queries; applicability only
    median_gap: float             # over measurable points, so read with coverage
    frac_below_grower_th: float   # inside their 2.0 vx rejection threshold
    frac_flagged: float           # inside our wider 6 vx flag
    largest_blob: int             # biggest contiguous flagged region, grid cells
    blob_fraction: float          # largest blob / n_queries -- USE THIS, not
                                  # largest_blob, for any cross-trace comparison
    flagged_given_neighbour: float  # flagged / n_measurable: the OLD conditional
                                    # number, kept only so published figures can
                                    # be traced. Never compare it across traces.
    largest_blob_4c: int          # same, 4-connectivity (8-connectivity above)
    top5_share: float             # concentration: top-5 blobs / all flagged
    du_p10: float
    du_median: float
    valid: bool
    reason: str

    def as_dict(self) -> dict:
        return asdict(self)


class _Entry:
    """Minimal shim so a bare path can be handed to atlas.write_atlas."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.winding = None


def analyse(path: Path | str, name: str = "", stride: int = 3,
            exclude_u: int | str = 60, workdir: Path | None = None,
            threads: int = 12) -> Result | None:
    """Run the self-gap analysis on one tifxyz surface.

    Returns None if the surface is too small to carry a previous wrap at all --
    a single-wrap patch has nothing to compare against, and reporting a number
    for it would be meaningless rather than merely noisy.
    """
    path = Path(path)
    work = Path(workdir) if workdir else Path("out")
    work.mkdir(parents=True, exist_ok=True)

    if exclude_u == "auto":
        est = estimate_exclude_u(path, stride, work, threads)
        if est is None:
            return None          # period not measurable: nothing to compare against
        exclude_u = est

    if path.suffix.lower() == ".obj":
        return _analyse_obj(path, name=name, stride=stride, exclude_u=exclude_u,
                            work=work, threads=threads)

    surf = tifxyz.read(path)
    v, u = np.nonzero(surf.valid[::stride, ::stride])
    V, U = v * stride, u * stride
    if len(V) < 5000 or surf.shape[1] < 4 * exclude_u:
        return None

    atlas.write_atlas([_Entry(path)], work / "sg_atlas.bin")
    atlas.write_queries_grouped(surf.points[V, U], U, work / "sg_query.bin")
    res = atlas.run_engine(work / "sg_atlas.bin", work / "sg_query.bin",
                           work / "sg_result.bin", threads=threads,
                           exclude_u=exclude_u)
    d, w1 = res["d1"], res["w1"]
    ok = np.isfinite(d)

    # Low coverage is a RESULT, not a failure. A surface with no neighbouring
    # wrap -- a flat sheet, a single-winding segment -- returns infinity
    # everywhere, and that is precisely the null control: the detector found
    # nothing because there is nothing. Bailing out with None here silently
    # removed the null control from any report generated through this function.
    if ok.sum() == 0:
        return Result(
            name=name or path.parent.parent.name, u_extent=int(surf.shape[1]),
            n_queries=int(d.size), n_measurable=0, coverage=0.0,
            median_gap=float("nan"),
            frac_below_grower_th=0.0, frac_flagged=0.0, largest_blob=0,
            blob_fraction=0.0, flagged_given_neighbour=0.0,
            largest_blob_4c=0, top5_share=0.0,
            du_p10=float("nan"), du_median=float("nan"), valid=True,
            reason="no neighbouring wrap within search radius (null control)",
        )

    flagged = ok & (d < FLAG_TH)
    H, W = surf.valid[::stride, ::stride].shape
    img = np.zeros((H, W), bool)
    img[v, u] = flagged
    # 8-connectivity headline, 4-connectivity reported alongside: connectivity
    # choice changes component counts, so quoting only one invites the reader to
    # assume the other agrees.
    lab, n = ndimage.label(img, structure=np.ones((3, 3)))
    if n:
        sizes = np.bincount(lab.ravel())[1:]
        order = np.argsort(-sizes)
        largest = int(sizes[order[0]])
        top5 = float(sizes[order[:5]].sum() / sizes.sum())
    else:
        largest, top5 = 0, 0.0
    lab4, n4 = ndimage.label(img)  # default structure = 4-connectivity
    largest_4c = int(np.bincount(lab4.ravel())[1:].max()) if n4 else 0

    if flagged.sum() > 20:
        du = np.abs(w1[flagged] - U[flagged]).astype(float)
        du_p10, du_med = float(np.percentile(du, 10)), float(np.median(du))
    else:
        du_p10 = du_med = float("nan")

    # Validity: flagged partners must be a real wrap away, not just past the
    # cutoff. Note the asymmetry -- this filter can only invalidate a POSITIVE.
    # A trace with no flags has nothing to be wrong about, so it passes clean;
    # calling that "rejected" would report the null control as a failure.
    n_queries = int(d.size)
    frac_flagged = flagged.sum() / n_queries
    coverage = float(ok.mean())
    if flagged.sum() == 0:
        valid, reason = True, "clean: no flagged regions"
    elif not np.isfinite(du_p10):
        valid, reason = True, "too few flagged points to characterise, treated as clean"
    elif du_p10 < DU_SAFETY * exclude_u:
        valid, reason = False, (
            f"|du| p10 = {du_p10:.0f} is within {DU_SAFETY}x the exclusion cutoff "
            f"({exclude_u}); flags may be same-wrap artifacts"
        )
    else:
        valid, reason = True, "ok"

    return Result(
        name=name or path.parent.parent.name,
        u_extent=int(surf.shape[1]),
        n_queries=n_queries,
        n_measurable=int(ok.sum()),
        coverage=coverage,
        median_gap=float(np.median(d[ok])),
        frac_below_grower_th=float((ok & (d < GROWER_TH)).sum() / n_queries),
        frac_flagged=float(frac_flagged),
        largest_blob=largest,
        blob_fraction=float(largest / n_queries),
        flagged_given_neighbour=float(flagged.sum() / ok.sum()),
        largest_blob_4c=largest_4c,
        top5_share=top5,
        du_p10=du_p10,
        du_median=du_med,
        valid=valid,
        reason=reason,
    )


def _analyse_obj(path: Path, *, name: str, stride: int, exclude_u: int,
                 work: Path, threads: int) -> Result | None:
    """Self-gap on an unstructured OBJ mesh.

    Same invariant, same kernel. The differences are that exclusion runs against
    the mesh's own `vt` u-coordinate rather than a grid column, and that the
    spatial-concentration statistics are unavailable: connected components need
    a grid, and a triangle soup has no rows to label. Those fields are reported
    as zero rather than guessed at, so an OBJ row is never silently compared
    against a tifxyz row on a statistic only one of them has.
    """
    from . import objmesh

    mesh = objmesh.read(path)
    if not mesh.has_param:
        return None                      # no vt: self-gap is undefined here
    pts, tags = objmesh.sample_points(mesh, stride=stride)
    if len(pts) < 2000:
        return None

    atlas.write_atlas_mesh(mesh, work / "sg_atlas.bin")
    atlas.write_queries_grouped(pts, tags, work / "sg_query.bin")
    res = atlas.run_engine(work / "sg_atlas.bin", work / "sg_query.bin",
                           work / "sg_result.bin", threads=threads,
                           exclude_u=exclude_u)
    d, w1 = res["d1"], res["w1"]
    ok = np.isfinite(d)
    if ok.sum() == 0:
        return Result(
            name=name or path.stem, u_extent=int(mesh.tri_u.max() - mesh.tri_u.min()),
            n_queries=int(d.size), n_measurable=0, coverage=0.0,
            median_gap=float("nan"),
            frac_below_grower_th=0.0, frac_flagged=0.0, largest_blob=0,
            blob_fraction=0.0, flagged_given_neighbour=0.0,
            largest_blob_4c=0, top5_share=0.0,
            du_p10=float("nan"), du_median=float("nan"), valid=True,
            reason="no neighbouring wrap within search radius (null control)",
        )

    flagged = ok & (d < FLAG_TH)
    n_queries = int(d.size)
    frac_flagged = flagged.sum() / n_queries
    if flagged.sum() > 20:
        du = np.abs(w1[flagged] - tags[flagged]).astype(float)
        du_p10, du_med = float(np.percentile(du, 10)), float(np.median(du))
    else:
        du_p10 = du_med = float("nan")

    if flagged.sum() == 0:
        valid, reason = True, "clean: no flagged regions"
    elif not np.isfinite(du_p10):
        valid, reason = True, "too few flagged points to characterise, treated as clean"
    elif du_p10 < DU_SAFETY * exclude_u:
        valid, reason = False, (
            f"|du| p10 = {du_p10:.0f} is within {DU_SAFETY}x the exclusion cutoff "
            f"({exclude_u}); flags may be same-wrap artifacts"
        )
    else:
        valid, reason = True, "ok (obj: concentration statistics unavailable)"

    return Result(
        name=name or path.stem,
        u_extent=int(mesh.tri_u.max() - mesh.tri_u.min()),
        n_queries=n_queries,
        n_measurable=int(ok.sum()),
        coverage=float(ok.mean()),
        median_gap=float(np.median(d[ok])),
        frac_below_grower_th=float((ok & (d < GROWER_TH)).sum() / n_queries),
        frac_flagged=float(frac_flagged),
        largest_blob=0, blob_fraction=0.0,
        flagged_given_neighbour=float(flagged.sum() / ok.sum()),
        largest_blob_4c=0, top5_share=0.0,
        du_p10=du_p10, du_median=du_med, valid=valid, reason=reason,
    )
