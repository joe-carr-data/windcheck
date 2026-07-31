"""Plant known sheet-switch defects in a clean trace, to measure detection.

Nothing in the published data says which segments are wrong, which is why three
attempts to confirm the detector against the CT volume were inconclusive. The
way out is to stop looking for ground truth and manufacture it: take a trace
that the detector reports as clean, copy a patch of its own previous wrap on top
of it, and now the answer is known exactly.

The defect injected is the "doubled back" case: the trace goes once around and
lands on the wrap it already traced, instead of advancing to the next one.
Concretely, for a patch of grid cells,

    points[v, u]  <-  points[v, u - T]

where T is the trace's revolution period in grid columns. After that the surface
at u genuinely coincides with the surface at u-T, which is exactly what the
detector is meant to find.

Two things this does NOT simulate, recorded so results are not overread:

  * A real doubled-back trace drifts onto the wrong wrap gradually and follows
    it with its own local geometry. This copies the neighbour verbatim, so the
    coincidence is exact rather than approximate. `jitter` adds noise to soften
    that, but it is still a friendlier target than reality.
  * The patch boundary is a discontinuity in the surface. Feathering blends it,
    but a mesh-quality check could in principle find the seam rather than the
    defect. The detector here uses no smoothness signal, so it cannot cheat that
    way -- but a different detector might.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile

from . import tifxyz


@dataclass
class Defect:
    """A planted defect and its exact extent, in grid coordinates."""

    v0: int
    v1: int
    u0: int
    u1: int
    period: int

    @property
    def n_cells(self) -> int:
        return (self.v1 - self.v0) * (self.u1 - self.u0)

    def mask(self, shape: tuple[int, int]) -> np.ndarray:
        m = np.zeros(shape, bool)
        m[self.v0:self.v1, self.u0:self.u1] = True
        return m


def estimate_period(surf: tifxyz.Surface, partner_u: np.ndarray,
                    query_u: np.ndarray, flagged: np.ndarray) -> int:
    """Revolution period in grid columns, from the detector's own partners."""
    if flagged.sum() < 20:
        raise ValueError("not enough flagged points to estimate the period")
    return int(np.median(np.abs(partner_u[flagged] - query_u[flagged])))


def inject(src: Path, dst: Path, defects: list[Defect], feather: int = 6,
           jitter: float = 0.0, rng: np.random.Generator | None = None) -> None:
    """Write a copy of `src` with each defect's patch replaced by its own
    previous wrap.

    `feather` blends the patch edge over that many cells so the seam is not a
    step function. `jitter` adds Gaussian noise in voxels, so the coincidence is
    approximate rather than exact.
    """
    rng = rng or np.random.default_rng(0)
    surf = tifxyz.read(src)
    P = surf.points.copy()
    valid = surf.valid.copy()
    H, W = surf.shape

    for d in defects:
        v0, v1 = max(0, d.v0), min(H, d.v1)
        u0, u1 = max(d.period, d.u0), min(W, d.u1)
        if u1 <= u0 or v1 <= v0:
            continue
        tgt = P[v0:v1, u0:u1]
        srcpatch = P[v0:v1, u0 - d.period:u1 - d.period].copy()
        ok = valid[v0:v1, u0 - d.period:u1 - d.period] & valid[v0:v1, u0:u1]

        if jitter > 0:
            srcpatch = srcpatch + rng.normal(0.0, jitter, srcpatch.shape).astype(np.float32)

        # blend weight: 1 in the interior, ramping to 0 across `feather` cells
        hh, ww = tgt.shape[:2]
        wv = np.minimum(np.arange(hh), np.arange(hh)[::-1]) / max(feather, 1)
        wu = np.minimum(np.arange(ww), np.arange(ww)[::-1]) / max(feather, 1)
        w = np.clip(np.minimum.outer(wv, wu), 0.0, 1.0)[..., None]

        blended = np.where(ok[..., None], w * srcpatch + (1 - w) * tgt, tgt)
        P[v0:v1, u0:u1] = blended

    dst.mkdir(parents=True, exist_ok=True)
    for i, axis in enumerate(("x", "y", "z")):
        arr = np.where(valid, P[..., i], -1.0).astype(np.float32)
        tifffile.imwrite(dst / f"{axis}.tif", arr)
    from .tifxyz import write_meta
    write_meta(src, dst, np.stack(
        [np.where(valid, P[..., i], -1.0).astype(np.float32) for i in range(3)],
        axis=-1), valid)


def score(detected: np.ndarray, truth: np.ndarray) -> dict:
    """Cell-level precision/recall of a detection mask against planted truth."""
    tp = int((detected & truth).sum())
    fp = int((detected & ~truth).sum())
    fn = int((~detected & truth).sum())
    prec = tp / (tp + fp) if tp + fp else float("nan")
    rec = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * prec * rec / (prec + rec) if tp and prec + rec else 0.0
    return {"tp": tp, "fp": fp, "fn": fn, "precision": prec, "recall": rec, "f1": f1}
