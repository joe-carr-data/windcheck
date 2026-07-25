"""Does a self-gap violation correspond to a real segmentation error?

The self-gap detector finds regions where a trace is coincident with itself at
exactly one revolution. That is a solid geometric fact. Whether it constitutes
an ERROR is a separate claim, and this module is the attempt to settle it
against the CT volume.

The test is the *missed-sheet* signature. Along the surface normal at a point,
find every bright band (each is a papyrus sheet) and ask which of them the trace
occupies. Occupancy is measured as distance from the ray point to the trace
surface, computed by the C++ engine -- not as "is there a grid sample nearby",
which fails because the grid samples every ~20 voxels laterally and so misses
even the sheet the trace is sitting on.

If a trace doubled back onto a wrap it had already traced, it never covered the
wrap it should have advanced to, so some nearby band should be orphaned. If
flagged and clean regions orphan bands at the same rate, the self-gap signal is
geometric only and does not evidence an error.

Sanity check that must pass before any number here is read: distance-to-surface
at t=0 must be ~0, since the trace is by construction sitting there.
"""

from __future__ import annotations

import numpy as np

from . import atlas

TS = np.arange(-28, 28.5, 1.0)      # sample offsets along the normal, voxels
BAND_FRAC = 0.55                    # band threshold, fraction of profile range
OCCUPIED_VX = 4.0                   # a band is occupied if the trace is this close
LATERAL_BOX = 70                    # half-size of the CT box fetched per seed


def normal_at(surf, vv: int, uu: int):
    """Unit surface normal from the tifxyz grid tangents."""
    P = surf.points
    du = P[vv, min(uu + 2, P.shape[1] - 1)] - P[vv, max(uu - 2, 0)]
    dv = P[min(vv + 2, P.shape[0] - 1), uu] - P[max(vv - 2, 0), uu]
    n = np.cross(du, dv)
    m = np.linalg.norm(n)
    return n / m if m > 1e-6 else None


def band_centres(profile: np.ndarray) -> list[int]:
    """Indices of bright bands in a CT intensity profile."""
    sm = np.convolve(np.nan_to_num(profile, nan=float(np.nanmean(profile))),
                     np.ones(5) / 5, "same")
    above = sm > sm.min() + BAND_FRAC * max(np.ptp(sm), 1.0)
    out, i = [], 0
    while i < len(above):
        if above[i]:
            k = i
            while k + 1 < len(above) and above[k + 1]:
                k += 1
            out.append((i + k) // 2)
            i = k + 1
        else:
            i += 1
    return out


def sample_profiles(surf, z0, mask, V, U, rng, n_boxes: int, per_box: int):
    """Fetch CT profiles along the normal at points drawn from `mask`."""
    idxs = np.nonzero(mask)[0]
    if not len(idxs):
        return np.zeros((0, len(TS))), np.zeros((0, 3))
    seeds = rng.choice(idxs, size=min(n_boxes, len(idxs)), replace=False)
    profs, rays = [], []
    for seed in seeds:
        c = surf.points[V[seed], U[seed]].astype(float)
        lo = np.floor(c - LATERAL_BOX).astype(int)[::-1]
        hi = lo + 2 * LATERAL_BOX + 1
        try:
            box = np.asarray(z0[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]).astype(np.float32)
        except (OSError, ValueError, KeyError, IndexError):
            # Streaming a chunk from S3 fails in several unrelated ways (network,
            # missing chunk, out-of-range slice). Any of them means "no data for
            # this seed", and skipping is correct -- but the sample is dropped,
            # not silently zero-filled, so it cannot bias the comparison.
            continue
        near = idxs[np.all(np.abs(surf.points[V[idxs], U[idxs]] - c) < 45, axis=1)]
        if not len(near):
            continue
        for j in rng.choice(near, size=min(per_box, len(near)), replace=False):
            vv, uu = V[j], U[j]
            p = surf.points[vv, uu].astype(float)
            n = normal_at(surf, vv, uu)
            if n is None:
                continue
            pts = p[None, :] + TS[:, None] * n[None, :]
            q = np.round(pts[:, ::-1]).astype(int) - lo
            good = np.all((q >= 0) & (q < box.shape), axis=1)
            if good.mean() < 0.8:
                continue
            pr = np.full(len(TS), np.nan)
            pr[good] = box[q[good, 0], q[good, 1], q[good, 2]]
            profs.append(pr)
            rays.append(pts)
    if not profs:
        return np.zeros((0, len(TS))), np.zeros((0, 3))
    return np.array(profs), np.concatenate(rays)


def occupancy(rays: np.ndarray, atlas_bin, qpath, rpath, threads: int = 12) -> np.ndarray:
    """Distance from every ray point to the trace surface, shaped (n, len(TS))."""
    atlas.write_queries(rays, qpath)
    r = atlas.run_engine(atlas_bin, qpath, rpath, threads=threads, max_dist=64.0)
    return r["d1"].reshape(-1, len(TS))


def score(profiles: np.ndarray, dists: np.ndarray) -> dict:
    """Band and orphan statistics for a set of sampled points."""
    nb, orp = [], []
    for pr, dd in zip(profiles, dists):
        cen = band_centres(pr)
        if not cen:
            continue
        nb.append(len(cen))
        orp.append(sum(1 for ci in cen if dd[ci] > OCCUPIED_VX))
    nb, orp = np.array(nb), np.array(orp)
    zi = len(TS) // 2
    return {
        "n": len(nb),
        "sanity_d_at_t0": float(np.median(dists[:, zi])) if len(dists) else float("nan"),
        "bands": float(nb.mean()) if len(nb) else float("nan"),
        "orphans": float(orp.mean()) if len(orp) else float("nan"),
        "orphan_rate": float(orp.sum() / max(nb.sum(), 1)),
        "frac_ge1": float(np.mean(orp >= 1)) if len(orp) else float("nan"),
        "_orp": orp,
    }
