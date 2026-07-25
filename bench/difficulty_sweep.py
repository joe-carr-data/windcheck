"""How hard can a planted defect get before the check stops finding it?

The first benchmark planted a verbatim copy of the neighbouring wrap. A real
doubled-back trace is not verbatim: it lands on the wrong sheet and then follows
it with its own local geometry, so the two surfaces are close but not identical.

`jitter` is the knob for that. At 0 the copy is exact; at the sheet spacing
(~17 vx) the "defect" is no longer meaningfully on the wrong wrap at all. Recall
must fall somewhere in between, and where it falls is the honest statement of
what this check can see.

Reported as a curve rather than a single number, because a single number from
the easiest setting would overstate it.
"""
import numpy as np
from pathlib import Path
from windcheck import atlas, tifxyz, inject

STR, EXU, FLAG = 3, 60, 6.0
WORK = Path("out/sweep"); WORK.mkdir(parents=True, exist_ok=True)

def detect(path):
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
hosts = {e.long_id.split('_')[-2]: e for e in ents if e.is_auto_grown}
rng = np.random.default_rng(11)

print("recall vs injection difficulty (precision in brackets)\n")
print(f"{'host':>6} {'period':>7} | " + " ".join(f"{j:>11}" for j in
      ["jitter 0.4","jitter 2","jitter 5","jitter 9","jitter 14"]))
print("-"*84)
for tag in ("0", "1", "3"):
    host = hosts[tag]
    base, w1, U, ok, s = detect(host.path)
    H, W = s.shape
    flagged_sub = base[np.nonzero(s.valid[::STR,::STR])]
    try:
        period = inject.estimate_period(s, w1, U, flagged_sub)
    except ValueError:
        print(f"{tag:>6}   (period not measurable)"); continue
    dv, du = 100, min(220, period - 40)
    v0, u0 = H//2 - dv//2, int(W*0.55)
    row = []
    for jit in (0.4, 2.0, 5.0, 9.0, 14.0):
        d = inject.Defect(v0, v0+dv, u0, u0+du, period)
        out = WORK/f"h{tag}_j{jit}"
        inject.inject(host.path, out, [d], feather=8, jitter=jit, rng=rng)
        m, *_ = detect(out)
        truth = (d.mask((H,W)) |
                 inject.Defect(d.v0,d.v1,d.u0-period,d.u1-period,period).mask((H,W)))[::STR,::STR]
        sc = inject.score(m & ~base, truth)
        row.append(f"{sc['recall']:.2f} [{sc['precision']:.2f}]")
    print(f"{tag:>6} {period:>7} | " + " ".join(f"{c:>11}" for c in row))
print("\nsheet spacing is ~17 vx: at jitter 14 the copy is barely on the wrong wrap.")
