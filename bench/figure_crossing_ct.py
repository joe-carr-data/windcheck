"""The picture: a CT cross-section with the trace drawn on it, jumping wraps.

Every other figure in this project argues. This one shows. A reader who knows
what a rolled scroll looks like can see a sheet switch without reading a number,
and can disagree with us using their own eyes -- which is the honest form of
persuasion when, as here, the statistical claim that a crossing IS an error has
not been established.

What is drawn:

  * one z-slice of the published CT volume, at full resolution, so individual
    papyrus sheets are 12-17 pixels apart and clearly separate;
  * every point of the trace that lies within a couple of voxels of that plane,
    coloured by its position ALONG THE TRACE (its `u` column);
  * the flagged crossing marked.

The colour is the whole point. On a correct trace the colours sweep smoothly
outwards, one wrap at a time. Where the trace has switched sheets, two very
different colours appear at the same place -- the trace has returned to a wrap it
already traced, thousands of columns later.

The volume is read directly from the public S3 zarr, uncompressed uint8 in
128^3 chunks, so only the chunks covering the window are fetched (~50 MB) rather
than the 7.6 GB a full slice would need.

    uv run --extra viz python bench/figure_crossing_ct.py
"""
from __future__ import annotations

import argparse
import io
import json
import subprocess
from pathlib import Path

import numpy as np

BUCKET = "s3://vesuvius-challenge-open-data"
ZARR = (f"{BUCKET}/PHerc0172/volumes/"
        "20241024131839-7.910um-53keV-masked.zarr")
CHUNK = 128
SHAPE = (20820, 6700, 9100)          # z, y, x at level 0


def fetch_chunk(level: int, zc: int, yc: int, xc: int) -> np.ndarray | None:
    """One 128^3 uint8 chunk, or None where the volume is empty there."""
    key = f"{ZARR}/{level}/{zc}/{yc}/{xc}"
    r = subprocess.run(["aws", "s3", "cp", key, "-", "--no-sign-request"],
                       capture_output=True)
    if r.returncode != 0 or len(r.stdout) != CHUNK ** 3:
        return None
    return np.frombuffer(r.stdout, dtype=np.uint8).reshape(CHUNK, CHUNK, CHUNK)


def slice_window(z: int, y0: int, y1: int, x0: int, x1: int,
                 level: int = 0) -> np.ndarray:
    """Assemble one z-plane over [y0:y1, x0:x1] from the chunks covering it."""
    out = np.zeros((y1 - y0, x1 - x0), dtype=np.uint8)
    zc, zi = divmod(z, CHUNK)
    for yc in range(y0 // CHUNK, (y1 - 1) // CHUNK + 1):
        for xc in range(x0 // CHUNK, (x1 - 1) // CHUNK + 1):
            blk = fetch_chunk(level, zc, yc, xc)
            if blk is None:
                continue
            gy0, gx0 = yc * CHUNK, xc * CHUNK
            ty0, tx0 = max(y0, gy0), max(x0, gx0)
            ty1, tx1 = min(y1, gy0 + CHUNK), min(x1, gx0 + CHUNK)
            if ty1 <= ty0 or tx1 <= tx0:
                continue
            out[ty0 - y0:ty1 - y0, tx0 - x0:tx1 - x0] = \
                blk[zi, ty0 - gy0:ty1 - gy0, tx0 - gx0:tx1 - gx0]
    return out


def plane_curve(P: np.ndarray, V: np.ndarray, z: float):
    """Where the surface crosses the plane z = const, one point per column.

    Sampling mesh vertices that happen to lie near the plane gives scattered
    dots, because the grid pitch is ~20 voxels and almost nothing lands exactly
    on a slice. The surface's actual intersection with the plane is a curve, so
    it is computed as one: within each column, find the row where the z
    coordinate crosses the plane and interpolate between the two neighbouring
    vertices.

    Returns (x, y, u) arrays. Columns whose rows never straddle the plane
    contribute nothing, which is correct -- the surface is not there.
    """
    xs, ys, us = [], [], []
    Z = P[:, :, 2]
    for u in range(P.shape[1]):
        col = V[:, u]
        if col.sum() < 2:
            continue
        rows = np.nonzero(col)[0]
        zc = Z[rows, u]
        d = zc - z
        sign_change = np.nonzero(np.diff(np.sign(d)) != 0)[0]
        for i in sign_change:
            r0, r1 = rows[i], rows[i + 1]
            if r1 - r0 > 2:              # a hole in the grid, not a crossing
                continue
            d0, d1 = d[i], d[i + 1]
            t = 0.0 if d1 == d0 else float(-d0 / (d1 - d0))
            xs.append(P[r0, u, 0] + t * (P[r1, u, 0] - P[r0, u, 0]))
            ys.append(P[r0, u, 1] + t * (P[r1, u, 1] - P[r0, u, 1]))
            us.append(u)
    return np.array(xs), np.array(ys), np.array(us)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--segment", default="data/scroll5_tifxyz/"
                    "20251115002745-auto_grown_20251115002740308_5_flatboi")
    ap.add_argument("--volume", default="20241024131839")
    ap.add_argument("--pairs", default="out/check/"
                    "20251115002745-auto_grown_20251115002740_pairs.csv")
    ap.add_argument("--half", type=int, default=320,
                    help="half-width of the CT window, in voxels")
    ap.add_argument("--zthick", type=float, default=2.5,
                    help="mesh points within this many voxels of the plane")
    ap.add_argument("--out", type=Path, default=Path("out/fig_crossing_ct.png"))
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from windcheck import tifxyz

    from windcheck.check import load_pairs

    mesh = sorted(Path(a.segment).glob(f"mesh/*{a.volume}*.tifxyz"))[0]
    s = tifxyz.read(mesh)
    P, V = s.points, s.valid

    # Centre on the MOST separated crossing, not the busiest slice. A busy slice
    # shows a healthy spiral -- several wraps at several radii, which is what a
    # correct trace looks like and proves nothing. The picture has to show the
    # case the measurement is about: two parts of the trace thousands of columns
    # apart meeting at one place.
    rec = load_pairs(Path(a.pairs))
    du = np.abs(rec["u1"].astype(np.int64) - rec["u2"].astype(np.int64))
    i = int(np.argmax(du))
    centre = P[int(rec["v1"][i]), int(rec["u1"][i])]
    z = int(round(float(centre[2])))
    print(f"widest crossing: u={rec['u1'][i]} and u={rec['u2'][i]}, "
          f"{du[i]} columns apart")
    x0, x1 = int(centre[0]) - a.half, int(centre[0]) + a.half
    y0, y1 = int(centre[1]) - a.half, int(centre[1]) + a.half
    x0, y0 = max(0, x0), max(0, y0)
    x1, y1 = min(SHAPE[2], x1), min(SHAPE[1], y1)
    print(f"crossing at x={centre[0]:.0f} y={centre[1]:.0f} z={z}")
    print(f"fetching CT window y[{y0}:{y1}] x[{x0}:{x1}] at z={z} ...")

    ct = slice_window(z, y0, y1, x0, x1)
    print(f"  window {ct.shape}, nonzero {100 * (ct > 0).mean():.1f}%")

    px, py, cu = plane_curve(P, V, z)
    inside = (px >= x0) & (px < x1) & (py >= y0) & (py < y1)
    px, py, cu = px[inside], py[inside], cu[inside]

    fig, ax = plt.subplots(figsize=(9.2, 8.4))
    ax.imshow(ct, cmap="gray", origin="upper",
              extent=[x0, x1, y1, y0], interpolation="nearest")
    sc = ax.scatter(px, py, c=cu, s=13, cmap="turbo", alpha=0.95, linewidths=0)
    ax.plot(centre[0], centre[1], "o", mfc="none", mec="white", mew=2.4, ms=26)
    ax.plot(centre[0], centre[1], "o", mfc="none", mec="red", mew=1.4, ms=26)

    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("position along the trace  (grid column u)")
    ax.set_xlabel("x (voxels)")
    ax.set_ylabel("y (voxels)")
    ax.set_title(
        f"One papyrus sheet, traced over and over\n"
        f"CT slice z={z}. Each dot is the trace, coloured by how far along "
        f"itself it is.\n"
        f"A correct trace crosses a given sheet once, so one sheet should carry "
        f"one colour.\n"
        f"Here the whole colour range lies on a single sheet: at the marked "
        f"point the trace\nreturns to it {int(du[i])} columns "
        f"({du[i] / 518.6:.1f} revolutions) later.", fontsize=10)
    fig.tight_layout()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=165)
    print(f"wrote {a.out}")

    if len(cu):
        print(f"  trace columns present in this plane: "
              f"{cu.min()} - {cu.max()}  (span {cu.max() - cu.min()})")


if __name__ == "__main__":
    main()
