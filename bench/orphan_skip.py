"""GATE CONDITION 1: are there zero sheets between known-adjacent wraps?

The 44 labelled windings of PHerc0172 are consecutive sheets, w_k and w_{k+1}.
Between them there is, by definition, no third sheet. Walking the CT from a
point on w_k to the nearest point on w_{k+1} must therefore cross no bright
band strictly in between.

If it does, the band detector is picking up fibre texture or noise rather than
sheets, and no orphan test built on it can work. This is the cheapest way to
find that out.
"""
import numpy as np, contextlib, io
from pathlib import Path
from windcheck import atlas, tifxyz
import zarr, s3fs

fs = s3fs.S3FileSystem(anon=True)
Z = zarr.open(s3fs.S3Map("vesuvius-challenge-open-data/PHerc0172/volumes/"
     "20241024131839-7.910um-53keV-masked.zarr", s3=fs, check=False), mode="r")["0"]
WORK = Path("out/orphan"); WORK.mkdir(parents=True, exist_ok=True)

def normal_at(s, vv, uu):
    P = s.points
    du = P[vv, min(uu+2, P.shape[1]-1)] - P[vv, max(uu-2, 0)]
    dv = P[min(vv+2, P.shape[0]-1), uu] - P[max(vv-2, 0), uu]
    n = np.cross(du, dv); m = np.linalg.norm(n)
    return n/m if m > 1e-6 else None

def bands_between(box, lo, p, direction, d_stop, margin=3.0, frac=0.55):
    """Count bright bands strictly between the surface and its neighbour."""
    ts = np.arange(0.0, d_stop + 0.5, 0.5)
    if len(ts) < 8: return None
    prof = []
    for t in ts:
        q = np.round((p + t*direction)[::-1]).astype(int) - lo
        prof.append(box[tuple(q)] if np.all((q >= 0) & (q < box.shape)) else np.nan)
    prof = np.array(prof, float)
    if np.isnan(prof).mean() > 0.2: return None
    sm = np.convolve(np.nan_to_num(prof, nan=float(np.nanmean(prof))), np.ones(5)/5, "same")
    rng_ = max(np.ptp(sm), 1.0)
    above = sm > sm.min() + frac*rng_
    inner = (ts > margin) & (ts < d_stop - margin)
    a = above & inner
    return int((np.diff(a.astype(int)) == 1).sum() + (a[0] if a.size else 0))

ents = atlas.discover(Path("data/scroll5_tifxyz"))
refs = {e.winding: e for e in ents if e.winding is not None}
rng = np.random.default_rng(2)
counts, gaps = [], []

for w in (56, 62, 70, 78, 84, 90):
    if w not in refs or w+1 not in refs: continue
    s = tifxyz.read(refs[w].path)
    atlas.write_atlas([refs[w+1]], WORK/"a.bin")           # the NEXT sheet only
    v, u = np.nonzero(s.valid[::6, ::6]); V, U = v*6, u*6
    pts = s.points[V, U]
    atlas.write_queries(pts, WORK/"q.bin")
    with contextlib.redirect_stderr(io.StringIO()):
        r = atlas.run_engine(WORK/"a.bin", WORK/"q.bin", WORK/"r.bin", threads=12)
    d = r["d1"]
    good = np.isfinite(d) & (d > 8) & (d < 30)             # a plausible one-sheet gap
    idx = np.nonzero(good)[0]
    if len(idx) < 30: continue
    for seed in rng.choice(idx, size=min(4, len(idx)), replace=False):
        c = pts[seed]
        lo = np.floor(c - 45).astype(int)[::-1]; hi = lo + 91
        try: box = np.asarray(Z[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]).astype(np.float32)
        except Exception: continue
        near = idx[np.all(np.abs(pts[idx] - c) < 30, axis=1)]
        for j in rng.choice(near, size=min(25, len(near)), replace=False):
            n = normal_at(s, V[j], U[j])
            if n is None: continue
            # point the normal toward the neighbouring sheet
            for sign in (1.0, -1.0):
                probe = pts[j] + sign*n*d[j]
                atlas.write_queries(probe[None, :].astype(np.float32), WORK/"q1.bin")
                with contextlib.redirect_stderr(io.StringIO()):
                    rr = atlas.run_engine(WORK/"a.bin", WORK/"q1.bin", WORK/"r1.bin", threads=2)
                if rr["d1"][0] < 4.0: break
            else:
                continue
            nb = bands_between(box, lo, pts[j], sign*n, float(d[j]))
            if nb is not None:
                counts.append(nb); gaps.append(float(d[j]))

counts = np.array(counts)
print(f"pairs of KNOWN-ADJACENT sheets probed: {len(counts)}")
print(f"  median gap between them: {np.median(gaps):.1f} vx")
print(f"  bands found strictly between: median {np.median(counts):.1f}, "
      f"mean {counts.mean():.2f}")
print(f"  fraction reporting >=1 band: {np.mean(counts >= 1):.1%}")
print(f"\nCONDITION 1 (median 0 and <20% with >=1): "
      f"{'PASS' if np.median(counts) == 0 and np.mean(counts >= 1) < 0.20 else 'FAIL'}")
