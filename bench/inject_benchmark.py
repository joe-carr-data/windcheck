"""Blind benchmark: plant defects of known extent, measure detection.

PRE-REGISTERED GATE (written before running):
  PASS  recall >= 0.70 AND precision >= 0.50 on planted defects at the two
        larger sizes, AND the clean host's false-positive rate stays under 2%.
  FAIL  anything less -> the detector cannot localise a defect it was handed,
        so a corrector built on it has nothing to act on. Drop the corrector.
"""
import numpy as np
from pathlib import Path
from scipy import ndimage
from windcheck import atlas, tifxyz, inject

STR, EXU, FLAG = 3, 60, 6.0
WORK = Path("out/bench"); WORK.mkdir(parents=True, exist_ok=True)

def detect(path):
    """flagged mask on the full grid, plus partner-u and query-u arrays"""
    s = tifxyz.read(path)
    class E: pass
    e = E(); e.path = path; e.winding = None
    atlas.write_atlas([e], WORK/"a.bin")
    v,u = np.nonzero(s.valid[::STR,::STR]); V,U = v*STR, u*STR
    atlas.write_queries_grouped(s.points[V,U], U, WORK/"q.bin")
    r = atlas.run_engine(WORK/"a.bin", WORK/"q.bin", WORK/"r.bin",
                         threads=12, exclude_u=EXU)
    d, w1 = r["d1"], r["w1"]; ok = np.isfinite(d)
    H,W = s.valid[::STR,::STR].shape
    m = np.zeros((H,W), bool); m[v,u] = ok & (d < FLAG)
    return m, w1, U, ok, s

ents = atlas.discover(Path("data/scroll5_tifxyz"))
host = [e for e in ents if e.is_auto_grown and e.long_id.endswith("_0_flatboi")][0]
base_mask, w1, U, ok, s = detect(host.path)
H, W = s.shape
period = inject.estimate_period(s, w1, U, (base_mask[::1,::1][np.nonzero(s.valid[::STR,::STR])]))
print(f"host: {host.long_id[-12:]}   grid {H}x{W}   period ~{period} columns")
print(f"      baseline flagged: {base_mask.mean():.3%} of grid cells\n")

sizes = [(40, 120), (80, 200), (150, 240)]      # u extent must stay under the period
rng = np.random.default_rng(3)
print(f"{'defect (v x u)':>16} {'cells':>8} {'recall':>8} {'precision':>10} {'f1':>7}")
results = []
for (dv, du) in sizes:
    v0 = H//2 - dv//2
    u0 = int(W*0.55)
    d = inject.Defect(v0, v0+dv, u0, u0+du, period)
    out = WORK/f"inj_{dv}x{du}"
    inject.inject(host.path, out, [d], feather=6, jitter=0.4, rng=rng)
    m, *_ = detect(out)
    # score on the subsampled grid the detector actually works on
    # A doubled-back trace has TWO coincident regions: the patch, and the wrap
    # it was copied from. Flagging both is correct, so both are ground truth.
    truth_full = d.mask((H, W))
    src = inject.Defect(d.v0, d.v1, d.u0 - period, d.u1 - period, period)
    truth_full |= src.mask((H, W))
    truth = truth_full[::STR, ::STR]
    new = m & ~base_mask                      # detections not present before injection
    sc = inject.score(new, truth)
    results.append((dv, du, sc))
    print(f"{f'{dv} x {du}':>16} {d.n_cells:>8,} {sc['recall']:>8.2f} "
          f"{sc['precision']:>10.2f} {sc['f1']:>7.2f}")

print(f"\nfalse positives outside any defect, clean host: {base_mask.mean():.3%}")
big = [r for r in results if r[0] >= 80]
gate = all(r[2]['recall'] >= 0.70 and r[2]['precision'] >= 0.50 for r in big) and base_mask.mean() < 0.02
print(f"\nPRE-REGISTERED GATE: {'PASS' if gate else 'FAIL'}")
