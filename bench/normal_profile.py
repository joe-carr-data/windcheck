"""Tight packing or doubled back? Sample the CT along the surface normal.

This is the discriminator Paul's question needs, and it does not depend on any
resampled per-segment artifact -- it reads the scroll volume directly.

At a flagged point the trace is within 6 vx of another part of itself. Two
mutually exclusive explanations:

  tight packing   two DISTINCT sheets happen to lie ~5 vx apart, the trace
                  correctly follows one on each pass. Then a second sheet
                  really is there, and the intensity profile along the normal
                  shows a SECOND PEAK a few voxels out.

  doubled back    the trace visited ONE sheet twice. Then there is no
                  neighbouring sheet that close, and the profile shows a single
                  peak with the next sheet at the normal spacing (~17 vx).

So the test is simply: do flagged points have a near neighbour peak that
unflagged points lack? Averaged over thousands of points, this is robust to the
CT noise that sank earlier single-location attempts.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import s3fs
import zarr

from windcheck import atlas, selfgap, tifxyz

BUCKET = "vesuvius-challenge-open-data"
VOLUME = "20241024131839-7.910um-53keV-masked.zarr"
SEG = ("data/scroll5_tifxyz/20251115002745-auto_grown_20251115002740308_5_"
       "flatboi/mesh/20251115002745-on-20241024131839-7.91um.tifxyz")
REACH = 26          # vx each side of the surface
STEP = 1.0


def normals(P: np.ndarray, V: np.ndarray) -> np.ndarray:
    """Unit normals from central differences on the grid."""
    dv = np.zeros_like(P)
    du = np.zeros_like(P)
    dv[1:-1] = P[2:] - P[:-2]
    du[:, 1:-1] = P[:, 2:] - P[:, :-2]
    n = np.cross(du, dv)
    ln = np.linalg.norm(n, axis=-1, keepdims=True)
    return n / np.maximum(ln, 1e-6)


def main(n_boxes: int, per_box: int, box: int, seed: int, out: Path):
    s = tifxyz.read(Path(SEG))
    P, V = s.points, s.valid
    stride = 3
    work = Path("out/interp")
    ex = selfgap.estimate_exclude_u(Path(SEG), stride, work, 12)
    v, u = np.nonzero(V[::stride, ::stride])
    Vi, Ui = v * stride, u * stride
    pts = P[Vi, Ui]
    atlas.write_atlas([selfgap._Entry(Path(SEG))], work / "n_a.bin")
    atlas.write_queries_grouped(pts, Ui, work / "n_q.bin")
    res = atlas.run_engine(work / "n_a.bin", work / "n_q.bin", work / "n_r.bin",
                           threads=12, exclude_u=ex)
    d = res["d1"].astype(np.float64)
    flagged = np.isfinite(d) & (d < selfgap.FLAG_TH)
    print(f"exclude_u {ex}   flagged {flagged.sum()} / {len(d)}")

    N = normals(P, V)
    nrm = N[Vi, Ui]
    good = np.isfinite(nrm).all(1) & (np.linalg.norm(nrm, axis=1) > 0.5)
    rng = np.random.default_rng(seed)

    fs = s3fs.S3FileSystem(anon=True)
    vol = zarr.open(s3fs.S3Map(root=f"{BUCKET}/PHerc0172/volumes/{VOLUME}",
                               s3=fs, check=False), mode="r")["0"]
    offs = np.arange(-REACH, REACH + STEP, STEP)

    # One S3 read per BOX, not per point. Per-point reads pulled several 128^3
    # chunks each and ran at ~20 profiles/minute.
    #
    # This also fixes a confound: flagged and unflagged points are now drawn
    # from the SAME neighbourhood, so local packing, damage and contrast are
    # controlled. Comparing flagged points against unflagged ones drawn from
    # anywhere in the scroll would confound "flagged" with "wherever flags
    # happen to live".
    fidx = np.nonzero(flagged & good)[0]
    rng.shuffle(fidx)
    seeds = fidx[:n_boxes]
    results = {"flagged": [], "unflagged": []}

    for b, si in enumerate(seeds):
        c = pts[si]
        half = box // 2
        lo = np.array([c[2] - half, c[1] - half, c[0] - half])
        hi = lo + box
        if (lo < REACH).any() or (hi[0] >= vol.shape[0]
                                  or hi[1] >= vol.shape[1]
                                  or hi[2] >= vol.shape[2]):
            continue
        sub = vol[int(lo[0]):int(hi[0]), int(lo[1]):int(hi[1]),
                  int(lo[2]):int(hi[2])]

        inside = ((np.abs(pts[:, 0] - c[0]) < half - REACH)
                  & (np.abs(pts[:, 1] - c[1]) < half - REACH)
                  & (np.abs(pts[:, 2] - c[2]) < half - REACH) & good)
        nf = int((inside & flagged).sum())
        nu = int((inside & ~flagged).sum())
        for label, mask in (("flagged", inside & flagged),
                            ("unflagged", inside & ~flagged)):
            idx = np.nonzero(mask)[0]
            rng.shuffle(idx)
            for i in idx[:per_box]:
                xyz = pts[i][None, :] + offs[:, None] * nrm[i][None, :]
                zz = np.round(xyz[:, 2] - lo[0]).astype(int)
                yy = np.round(xyz[:, 1] - lo[1]).astype(int)
                xx = np.round(xyz[:, 0] - lo[2]).astype(int)
                if (zz.min() < 0 or yy.min() < 0 or xx.min() < 0
                        or zz.max() >= sub.shape[0]
                        or yy.max() >= sub.shape[1]
                        or xx.max() >= sub.shape[2]):
                    continue
                results[label].append(sub[zz, yy, xx].astype(np.float32))
        print(f"  box {b+1}/{len(seeds)}  in-box flagged {nf} unflagged {nu}"
              f"  totals f={len(results['flagged'])} "
              f"u={len(results['unflagged'])}", flush=True)

    results = {k: np.array(v) for k, v in results.items()}
    for k, v in results.items():
        print(f"{k}: {len(v)} profiles")

    summary = {}
    for label, prof in results.items():
        mean = prof.mean(0)
        summary[label] = {"n": int(len(prof)), "offsets": offs.tolist(),
                          "mean_profile": mean.tolist()}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))
    np.save(out.with_suffix(".flagged.npy"), results["flagged"])
    np.save(out.with_suffix(".unflagged.npy"), results["unflagged"])
    print("wrote", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--boxes", type=int, default=12)
    ap.add_argument("--per-box", type=int, default=60)
    ap.add_argument("--box", type=int, default=320)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("out/routeA/normal_profile.json"))
    a = ap.parse_args()
    main(a.boxes, a.per_box, a.box, a.seed, a.out)
