"""Calibration for test C: what AUC can raw vertices reach on a KNOWN defect?

Test C asks whether flagged regions separate from unflagged ones on
vertex-to-vertex distance, with no interpolation. It returned a median AUC of
0.816 on real traces, which the pre-registration calls "partial".

That number is uninterpretable on its own. At 20 vx sampling a genuinely
touching surface still has its nearest vertex up to half a grid pitch away
tangentially, so there is a ceiling on how well raw vertices *can* separate a
real contact. This measures that ceiling on planted defects, where the answer is
known exactly.

Read it honestly in both directions:
  ceiling ~0.82  -> real traces are at the ceiling; test C's "partial" is a
                    limit of the sampling, not evidence against the flags
  ceiling ~0.98  -> raw vertices could have separated a real defect cleanly,
                    and the real traces' 0.816 indicates genuine contamination

The injection is a FRIENDLIER target than reality (it copies the neighbouring
wrap verbatim, so coincidence is exact), which biases the ceiling UP. Jitter
softens that. A ceiling measured this way is therefore an upper bound, and the
gap it implies is an upper bound too.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from windcheck import atlas, inject, selfgap, tifxyz
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from interp_support import auc, vertex_distance

STRIDE = 3
WORK = Path("out/interp")
WORK.mkdir(parents=True, exist_ok=True)


def measure(path: Path, exclude_u: int):
    surf = tifxyz.read(path)
    v, u = np.nonzero(surf.valid[::STRIDE, ::STRIDE])
    V, U = v * STRIDE, u * STRIDE
    pts = surf.points[V, U]
    atlas.write_atlas([selfgap._Entry(path)], WORK / "c_a.bin")
    atlas.write_queries_grouped(pts, U, WORK / "c_q.bin")
    res = atlas.run_engine(WORK / "c_a.bin", WORK / "c_q.bin", WORK / "c_r.bin",
                           threads=12, exclude_u=exclude_u)
    d = res["d1"].astype(np.float64)
    return surf, v, u, V, U, pts, d, res["w1"]


if __name__ == "__main__":
    ents = atlas.discover(Path("data/scroll5_tifxyz"))
    host = [e for e in ents if e.is_auto_grown
            and e.long_id.endswith("_0_flatboi")][0]
    ex = selfgap.estimate_exclude_u(host.path, STRIDE, WORK, 12)
    print(f"host {host.long_id[-12:]}   exclude_u {ex}")

    surf, v, u, V, U, pts, d0, w1 = measure(host.path, ex)
    base_flag = np.isfinite(d0) & (d0 < selfgap.FLAG_TH)
    period = inject.estimate_period(surf, w1, U, base_flag)
    H, W = surf.shape
    print(f"grid {H}x{W}  period {period}  baseline flagged {base_flag.mean():.2%}")

    rows = []
    for jit in (0.4, 5.0, 9.0):
        dv, du = 150, 240
        v0 = H // 2 - dv // 2
        u0 = max(period + 10, W // 2 - du // 2)
        defect = inject.Defect(v0, v0 + dv, u0, u0 + du, period)
        dst = WORK / f"inj_j{jit}.tifxyz"
        inject.inject(host.path, dst, [defect], jitter=jit,
                      rng=np.random.default_rng(3))

        _, vi, ui, Vi, Ui, ptsi, di, _ = measure(dst, ex)
        truth = defect.mask(surf.shape)[Vi, Ui]
        dvert = vertex_distance(ptsi, Ui, ex)
        fin = np.isfinite(dvert)

        a_truth = auc(dvert[truth & fin], dvert[(~truth) & fin])
        det = np.isfinite(di) & (di < selfgap.FLAG_TH)
        recall = float(det[truth].mean())
        rows.append({"jitter": jit, "auc_vertex_vs_truth": a_truth,
                     "detector_recall_on_patch": recall,
                     "n_truth": int(truth.sum())})
        print(f"  jitter {jit:4.1f}   AUC(d_vertex, TRUE defect) = {a_truth:.3f}"
              f"   detector recall {recall:.2f}")

    print("\nreal traces (test C): median AUC 0.816, range 0.779-0.870")
    best = max(r["auc_vertex_vs_truth"] for r in rows)
    print(f"planted-defect ceiling: {best:.3f} (best over jitter levels)")
    print("ceiling is an UPPER bound: injections copy the wrap verbatim, "
          "which is friendlier than a real drift.")
    (WORK / "ceiling.json").write_text(json.dumps(rows, indent=2))
