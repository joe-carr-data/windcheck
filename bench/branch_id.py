"""Given a flagged region, can we tell WHICH branch is the wrong one?

Detection localises a defect. A corrector needs more: of the two coincident
branches, which one should be moved? On planted defects the answer is known --
the injected patch is wrong, the wrap it was copied from is right -- so this is
directly measurable.

A caution that shapes the reading: once two branches are coincident they occupy
the same 3D location and therefore see the same CT neighbourhood. Anything
purely local is expected to fail. Three signals are measured:

  ct       mean CT intensity under each branch. Expected to fail, and included
           precisely so that failure is on the record rather than assumed.
  seam     size of the 3D position jump across each branch's boundary in u.
           A spliced patch has a discontinuity where it was inserted; the
           original does not. This is expected to work HERE and to generalise
           badly, because a real doubled-back trace drifts in smoothly.
  normal   agreement between each branch's surface normal and the local trend
           of the trace either side of it.

PASS if some signal picks the wrong branch in >=80% of trials. Below that a
corrector cannot be driven by it.
"""
import numpy as np
from pathlib import Path
from windcheck import atlas, tifxyz, inject
import zarr, s3fs

STR = 3
WORK = Path("out/branch"); WORK.mkdir(parents=True, exist_ok=True)
fs = s3fs.S3FileSystem(anon=True)
Z = zarr.open(s3fs.S3Map("vesuvius-challenge-open-data/PHerc0172/volumes/"
     "20241024131839-7.910um-53keV-masked.zarr", s3=fs, check=False), mode="r")["0"]

def ct_mean(pts, n_boxes=6, rng=None):
    """CT intensity under a branch, via small boxes around sampled seeds.

    A branch spans thousands of voxels, so it cannot be fetched as one box.
    Sample seeds, pull a modest box around each, and average the intensity at
    whichever branch points fall inside it.
    """
    rng = rng or np.random.default_rng(0)
    if len(pts) == 0: return np.nan
    vals = []
    for seed in rng.choice(len(pts), size=min(n_boxes, len(pts)), replace=False):
        c = pts[seed]
        lo = np.floor(c - 40).astype(int)[::-1]; hi = lo + 81
        try:
            box = np.asarray(Z[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]).astype(np.float32)
        except (OSError, ValueError, KeyError, IndexError):
            continue
        q = np.round(pts[:, ::-1]).astype(int) - lo
        ok = np.all((q >= 0) & (q < box.shape), axis=1)
        if ok.any():
            vals.append(box[q[ok,0], q[ok,1], q[ok,2]].mean())
    return float(np.mean(vals)) if vals else np.nan

def seam(P, valid, v0, v1, u, w=4):
    """3D jump across grid column boundary u."""
    a, b = u - 1, u
    if a < 0 or b >= P.shape[1]: return np.nan
    m = valid[v0:v1, a] & valid[v0:v1, b]
    if m.sum() < 10: return np.nan
    step = np.linalg.norm(P[v0:v1, b][m] - P[v0:v1, a][m], axis=1)
    ref = []
    for c in (a - w, b + w):
        if 0 <= c < P.shape[1]-1:
            mm = valid[v0:v1, c] & valid[v0:v1, c+1]
            if mm.sum() > 10:
                ref.append(np.median(np.linalg.norm(P[v0:v1,c+1][mm]-P[v0:v1,c][mm],axis=1)))
    base = np.median(ref) if ref else 1.0
    return float(np.median(step) / max(base, 1e-6))

ents = atlas.discover(Path("data/scroll5_tifxyz"))
hosts = {e.long_id.split('_')[-2]: e for e in ents if e.is_auto_grown}
rng = np.random.default_rng(5)
periods = {"0": 271, "1": 313, "3": 436}

print(f"{'host':>5} {'jit':>5} | {'ct wrong':>9} {'ct right':>9} {'ct picks':>9} "
      f"| {'seam wrong':>11} {'seam right':>11} {'seam picks':>11}")
print("-"*84)
votes = {"ct": [], "seam": []}
for tag in ("0", "1", "3"):
    host = hosts[tag]; T = periods[tag]
    s0 = tifxyz.read(host.path); H, W = s0.shape
    dv, du = 100, min(220, T - 40)
    v0, u0 = H//2 - dv//2, int(W*0.55)
    for jit in (0.4, 5.0):
        d = inject.Defect(v0, v0+dv, u0, u0+du, T)
        out = WORK/f"h{tag}_j{jit}"
        inject.inject(host.path, out, [d], feather=8, jitter=jit, rng=rng)
        s = tifxyz.read(out); P, V = s.points, s.valid
        def sample(uu0, uu1):
            m = V[v0:v1_, uu0:uu1] if False else V[v0:v0+dv, uu0:uu1]
            pts = P[v0:v0+dv, uu0:uu1][m]
            return pts[::7]
        wrong = sample(u0, u0+du)                 # the injected patch
        right = sample(u0-T, u0-T+du)             # the wrap it copied
        cw, cr = ct_mean(wrong, rng=rng), ct_mean(right, rng=rng)
        sw = seam(P, V, v0, v0+dv, u0)
        sr = seam(P, V, v0, v0+dv, u0-T)
        ct_ok = cw < cr                            # guess: wrong branch darker
        se_ok = (sw > sr) if np.isfinite(sw) and np.isfinite(sr) else False
        votes["ct"].append(ct_ok); votes["seam"].append(se_ok)
        print(f"{tag:>5} {jit:>5} | {cw:>9.1f} {cr:>9.1f} {str(ct_ok):>9} "
              f"| {sw:>11.2f} {sr:>11.2f} {str(se_ok):>11}")

print()
for k, v in votes.items():
    r = np.mean(v)
    print(f"  {k:>5}: picks the wrong branch in {r:.0%} of {len(v)} trials"
          f"   {'PASS' if r >= 0.8 else 'FAIL'}")
