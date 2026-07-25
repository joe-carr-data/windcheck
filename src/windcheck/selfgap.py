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


@dataclass
class Result:
    name: str
    u_extent: int
    n_points: int
    coverage: float               # fraction with a finite self-gap
    median_gap: float
    frac_below_grower_th: float   # inside their 2.0 vx rejection threshold
    frac_flagged: float           # inside our wider 6 vx flag
    largest_blob: int             # biggest contiguous flagged region, grid cells
    blob_fraction: float          # largest blob / analysed points -- USE THIS, not
                                  # largest_blob, for any cross-trace comparison:
                                  # raw counts scale with trace size and stride, and
                                  # comparing a raw blob at one threshold against a
                                  # rate at another is how a "58x anomaly" turned
                                  # into a ~1.8x one under review
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
            exclude_u: int = 60, workdir: Path | None = None,
            threads: int = 12) -> Result | None:
    """Run the self-gap analysis on one tifxyz surface.

    Returns None if the surface is too small to carry a previous wrap at all --
    a single-wrap patch has nothing to compare against, and reporting a number
    for it would be meaningless rather than merely noisy.
    """
    path = Path(path)
    work = Path(workdir) if workdir else Path("out")
    work.mkdir(parents=True, exist_ok=True)

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
            n_points=0, coverage=0.0, median_gap=float("nan"),
            frac_below_grower_th=0.0, frac_flagged=0.0, largest_blob=0,
            blob_fraction=0.0, largest_blob_4c=0, top5_share=0.0,
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
    frac_flagged = flagged.sum() / ok.sum()
    if frac_flagged < 1e-3:
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
        n_points=int(ok.sum()),
        coverage=float(ok.mean()),
        median_gap=float(np.median(d[ok])),
        frac_below_grower_th=float((ok & (d < GROWER_TH)).sum() / ok.sum()),
        frac_flagged=float(frac_flagged),
        largest_blob=largest,
        blob_fraction=float(largest / ok.sum()),
        largest_blob_4c=largest_4c,
        top5_share=top5,
        du_p10=du_p10,
        du_median=du_med,
        valid=valid,
        reason=reason,
    )
