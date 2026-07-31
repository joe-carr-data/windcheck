"""Route A / test A2: does surface texture break at a flag boundary?

Gate pre-registered in notes/PREREG-routeA.md. The descriptor (16 magnitude-
weighted gradient-orientation bins + 8 radial power-spectrum bins, L2
normalised) was fixed in that file before any number below was computed.

A2 is the test the corrector needs. The mesh is smooth across a doubling-back
fault -- the grower walked across a contact point, it did not jump -- so every
geometric seam statistic has scored chance. Texture is the channel where the
two sheets are different physical papyrus and cannot agree.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
from scipy import ndimage

from windcheck import atlas, selfgap, tifxyz
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from texture_a1 import open_surface_volume, LAYERS, SCALE

WIN = 20         # grid columns per side of the cut
N_ORI, N_RAD = 16, 8


def descriptor(img: np.ndarray) -> np.ndarray:
    """Orientation energy + radial spectrum, L2-normalised. Fixed pre-hoc."""
    gy, gx = np.gradient(img)
    mag = np.hypot(gx, gy)
    ang = np.arctan2(gy, gx) % np.pi                    # undirected
    ori = np.histogram(ang, bins=N_ORI, range=(0, np.pi), weights=mag)[0]

    f = np.abs(np.fft.fftshift(np.fft.fft2(img - img.mean())))
    h, w = f.shape
    yy, xx = np.mgrid[:h, :w]
    rr = np.hypot(yy - h / 2, xx - w / 2)
    rmax = min(h, w) / 2
    bins = np.clip((rr / rmax * N_RAD).astype(int), 0, N_RAD - 1)
    rad = np.array([f[bins == k].mean() for k in range(N_RAD)])

    d = np.concatenate([ori / (ori.sum() + 1e-9), rad / (rad.sum() + 1e-9)])
    return d / (np.linalg.norm(d) + 1e-9)


def window(arr, v0: int, v1: int, u0: int, u1: int):
    blk = arr[LAYERS, v0 * SCALE:v1 * SCALE, u0 * SCALE:u1 * SCALE]
    if blk.size == 0:
        return None
    img = blk.astype(np.float32).mean(0)
    if img.std() < 1e-3:
        return None
    img -= img.mean(1, keepdims=True)
    return img / (img.std(1, keepdims=True) + 1e-6)


def cut_distance(arr, v0: int, v1: int, x: int):
    L = window(arr, v0, v1, x - WIN, x)
    R = window(arr, v0, v1, x, x + WIN)
    if L is None or R is None:
        return None
    a, b = descriptor(L), descriptor(R)
    return float(1.0 - a @ b)


def analyse(seg: str, sample: str, n_null: int, seed: int) -> dict:
    d = glob.glob(f"data/{'scroll5' if sample == 'PHerc0172' else sample}_tifxyz/"
                  f"{seg}/mesh/*-on-*-7.91um.tifxyz")
    path = Path(sorted(d)[-1])
    stride = 3
    work = Path("out/routeA")
    ex = selfgap.estimate_exclude_u(path, stride, work, 12)
    s = tifxyz.read(path)
    v, u = np.nonzero(s.valid[::stride, ::stride])
    V, U = v * stride, u * stride
    atlas.write_atlas([selfgap._Entry(path)], work / "a2a.bin")
    atlas.write_queries_grouped(s.points[V, U], U, work / "a2q.bin")
    res = atlas.run_engine(work / "a2a.bin", work / "a2q.bin", work / "a2r.bin",
                           threads=12, exclude_u=ex)
    dd = res["d1"]
    ok = np.isfinite(dd)
    flag = ok & (dd < selfgap.FLAG_TH)
    H, W = s.valid[::stride, ::stride].shape
    img = np.zeros((H, W), bool)
    img[v, u] = flag
    lab, n = ndimage.label(img, structure=np.ones((3, 3)))
    if n == 0:
        return {"segment": seg, "skip": "no flagged region"}
    k = int(np.argmax(np.bincount(lab.ravel())[1:])) + 1
    blob = lab == k

    # densest v band, then the sharpest u onset inside it
    rows = blob.sum(1)
    band = max(range(0, H - 14), key=lambda a: rows[a:a + 14].sum())
    v0, v1 = band * stride, (band + 14) * stride
    cols = blob[band:band + 14].sum(0)
    onset = int(np.nonzero(cols)[0].min()) * stride

    arr = open_surface_volume(seg, sample)
    D_fault = cut_distance(arr, v0, v1, onset)
    if D_fault is None:
        return {"segment": seg, "skip": "fault window off-surface"}

    # null: cuts inside u ranges with no flagged cell anywhere in the band
    anyflag = img[band:band + 14].sum(0) > 0
    clean = np.nonzero(~anyflag)[0] * stride
    clean = clean[(clean > WIN + 5) & (clean < s.shape[1] - WIN - 5)]
    clean = clean[np.abs(clean - onset) > 3 * WIN]
    rng = np.random.default_rng(seed)
    rng.shuffle(clean)
    nulls = []
    for x in clean:
        if len(nulls) >= n_null:
            break
        val = cut_distance(arr, v0, v1, int(x))
        if val is not None:
            nulls.append(val)
    if len(nulls) < 20:
        return {"segment": seg, "skip": f"only {len(nulls)} null cuts"}

    nulls = np.array(nulls)
    pct = float((nulls < D_fault).mean() * 100)
    out = {"segment": seg, "exclude_u": int(ex), "v_band": [v0, v1],
           "fault_u": onset, "D_fault": D_fault,
           "null_n": len(nulls), "null_median": float(np.median(nulls)),
           "null_p95": float(np.percentile(nulls, 95)),
           "fault_percentile": pct,
           "verdict": "PASS" if pct >= 95 else "FAIL"}
    print(f"{seg[:46]:48s} v{v0}-{v1} u={onset:5d}  "
          f"D_fault {D_fault:.4f}  null med {out['null_median']:.4f} "
          f"p95 {out['null_p95']:.4f}  -> pct {pct:5.1f}  {out['verdict']}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("segments", nargs="+")
    ap.add_argument("--nulls", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--sample", default="PHerc0172")
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    rs = [analyse(s, a.sample, a.nulls, a.seed) for s in a.segments]
    good = [r for r in rs if "verdict" in r]
    n_pass = sum(r["verdict"] == "PASS" for r in good)
    print(f"\nA2 GATE: {n_pass}/{len(good)} boundaries at or above null p95 "
          f"(pre-registered threshold: 3 of 4)")
    if a.json:
        a.json.write_text(json.dumps(rs, indent=2))
