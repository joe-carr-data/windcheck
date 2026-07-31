"""Route A / test A1: are flagged regions duplicated *surface*, not just
duplicated geometry?

Gate is pre-registered in notes/PREREG-routeA.md. Nothing here is fitted: the
patch size, layer range and statistic were fixed before the first number.

The detector says two parts of a trace occupy the same place. That claim has
only ever been checked against the mesh geometry it was derived from. Here the
surface-volume channel -- the CT resampled along the sheet, 33 layers thick, at
20x the grid -- is asked to agree, and it did not see the geometry.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import s3fs
import zarr

from windcheck import atlas, selfgap, tifxyz

BUCKET = "vesuvius-challenge-open-data"
LEVEL = "1"          # grid x10
LAYERS = slice(10, 23)
PATCH = 12           # grid cells per side -> 120 px at level 1
SCALE = 10           # level-1 px per grid cell


def open_surface_volume(seg: str, sample: str = "PHerc0172"):
    fs = s3fs.S3FileSystem(anon=True)
    base = f"{BUCKET}/{sample}/segments/{seg}/surface-volumes/"
    sub = [p for p in fs.ls(base) if p.endswith(".zarr")]
    if not sub:
        raise FileNotFoundError(f"no surface-volume zarr under {base}")
    g = zarr.open(s3fs.S3Map(root=sub[0], s3=fs, check=False), mode="r")
    return g[LEVEL]


def flagged_with_partners(path: Path, stride: int, threads: int = 12):
    """Flagged full-grid cells and the partner column the engine matched."""
    work = Path("out/routeA")
    work.mkdir(parents=True, exist_ok=True)
    ex = selfgap.estimate_exclude_u(path, stride, work, threads)
    if ex is None:
        raise RuntimeError("revolution period not measurable")
    s = tifxyz.read(path)
    v, u = np.nonzero(s.valid[::stride, ::stride])
    V, U = v * stride, u * stride
    atlas.write_atlas([selfgap._Entry(path)], work / "a.bin")
    atlas.write_queries_grouped(s.points[V, U], U, work / "q.bin")
    res = atlas.run_engine(work / "a.bin", work / "q.bin", work / "r.bin",
                           threads=threads, exclude_u=ex)
    d, w1 = res["d1"], res["w1"]
    ok = np.isfinite(d)
    flag = ok & (d < selfgap.FLAG_TH)
    return s, V, U, flag, ok, w1, ex


def locate_partner(s, V0: int, U0: int, ucol: int, half: int = 2):
    """Grid cell nearest in 3D to (V0,U0), searched in columns ucol +- half."""
    lo, hi = max(0, ucol - half), min(s.shape[1], ucol + half + 1)
    sub = s.points[:, lo:hi]
    val = s.valid[:, lo:hi]
    if not val.any():
        return None
    dd = np.linalg.norm(sub - s.points[V0, U0], axis=-1)
    dd[~val] = np.inf
    vi, ui = np.unravel_index(np.argmin(dd), dd.shape)
    if not np.isfinite(dd[vi, ui]):
        return None
    return int(vi), int(lo + ui), float(dd[vi, ui])


def patch(arr, vg: int, ug: int):
    """Standardised texture patch centred on a grid cell, or None if clipped."""
    y0, x0 = (vg - PATCH // 2) * SCALE, (ug - PATCH // 2) * SCALE
    n = PATCH * SCALE
    if y0 < 0 or x0 < 0 or y0 + n > arr.shape[1] or x0 + n > arr.shape[2]:
        return None
    blk = arr[LAYERS, y0:y0 + n, x0:x0 + n]
    if blk.size == 0:
        return None
    img = blk.astype(np.float32).mean(0)
    img -= img.mean(1, keepdims=True)
    sd = img.std(1, keepdims=True)
    if (sd < 1e-3).mean() > 0.5:
        return None                      # empty / off-surface region
    return img / (sd + 1e-6)


def corr(a, b) -> float:
    a = a.ravel() - a.mean()
    b = b.ravel() - b.mean()
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-9))


def run(seg: str, n_pairs: int, seed: int, sample: str) -> dict:
    d = glob.glob(f"data/{'scroll5' if sample == 'PHerc0172' else sample}_tifxyz/"
                  f"{seg}/mesh/*-on-*-7.91um.tifxyz")
    path = Path(sorted(d)[-1])
    s, V, U, flag, ok, w1, ex = flagged_with_partners(path, stride=3)
    print(f"{seg}\n  grid {s.shape}  exclude_u {ex}  "
          f"flagged {flag.sum()}/{ok.sum()} ({flag.mean()*100:.2f}%)")

    arr = open_surface_volume(seg, sample)
    print(f"  surface volume {arr.shape}  ratio "
          f"{arr.shape[1]/s.shape[0]:.1f} x {arr.shape[2]/s.shape[1]:.1f}")

    rng = np.random.default_rng(seed)
    fi = np.nonzero(flag)[0]
    ni = np.nonzero(ok & ~flag)[0]
    rng.shuffle(fi)
    rng.shuffle(ni)

    # The null must be matched in column offset, otherwise it tests distance
    # rather than flag status. Each null cell is paired with the cell |du| away
    # that the flagged population actually uses.
    du_pop = np.abs(w1[flag] - U[flag])
    du_pop = du_pop[du_pop > ex]

    out = {"flagged": [], "null": []}
    for label, idx in (("flagged", fi), ("null", ni)):
        got = 0
        for i in idx:
            if got >= n_pairs:
                break
            v0, u0 = int(V[i]), int(U[i])
            ucol = (int(w1[i]) if label == "flagged"
                    else u0 + int(rng.choice(du_pop)) * rng.choice([-1, 1]))
            if not (0 <= ucol < s.shape[1]):
                continue
            p = locate_partner(s, v0, u0, ucol)
            if p is None:
                continue
            v1, u1, dist = p
            A, B = patch(arr, v0, u0), patch(arr, v1, u1)
            if A is None or B is None:
                continue
            out[label].append({"v0": v0, "u0": u0, "v1": v1, "u1": u1,
                               "d3d": dist, "r": corr(A, B)})
            got += 1
        print(f"  {label}: {got} pairs")

    rf = np.array([x["r"] for x in out["flagged"]])
    rn = np.array([x["r"] for x in out["null"]])
    if len(rf) == 0 or len(rn) == 0:
        print("  INSUFFICIENT PAIRS")
        return {"segment": seg, "n_flagged": len(rf), "n_null": len(rn)}

    p95 = float(np.percentile(rn, 95))
    med = float(np.median(rf))
    verdict = "PASS" if med > p95 else "FAIL"
    print(f"\n  flagged  median r = {med:+.3f}   "
          f"[p10 {np.percentile(rf,10):+.3f}  p90 {np.percentile(rf,90):+.3f}]")
    print(f"  null     median r = {np.median(rn):+.3f}   p95 = {p95:+.3f}")
    print(f"  d3d flagged median {np.median([x['d3d'] for x in out['flagged']]):.2f} vx"
          f"   null {np.median([x['d3d'] for x in out['null']]):.2f} vx")
    print(f"  GATE (median flagged r > null p95): {verdict}")
    return {"segment": seg, "exclude_u": int(ex), "n_flagged": len(rf),
            "n_null": len(rn), "median_r_flagged": med,
            "median_r_null": float(np.median(rn)), "null_p95": p95,
            "verdict": verdict, "pairs": out}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("segments", nargs="+")
    ap.add_argument("--pairs", type=int, default=60)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--sample", default="PHerc0172")
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    results = [run(s, a.pairs, a.seed, a.sample) for s in a.segments]
    if a.json:
        a.json.write_text(json.dumps(results, indent=2))
        print(f"\nwrote {a.json}")
