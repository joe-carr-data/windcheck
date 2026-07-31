"""The inspection the maintainers asked for: look at a flagged region in the CT.

Paul's question was whether the flagged regions are actually defective, since
wraps can genuinely sit very close together. That cannot be answered from the
mesh, and it cannot be answered from the per-segment surface volume either --
that is resampled along the mesh's own normal, so it inherits whatever the mesh
believes and shows chevrons rather than sheets.

This cuts a real axial slice from the scroll volume and draws the trace on top,
which is the view VC3D gives. Sheets appear as spiral arcs. The question then
becomes visual and direct: at a flagged location, does the trace sit on one arc
and come back onto the *same* arc a revolution later, or are there two distinct
arcs that happen to be very close?
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import s3fs
import zarr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from windcheck import tifxyz

BUCKET = "vesuvius-challenge-open-data"
VOLUME = "20241024131839-7.910um-53keV-masked.zarr"
SEG = ("data/scroll5_tifxyz/20251115002745-auto_grown_20251115002740308_5_"
       "flatboi/mesh/20251115002745-on-20241024131839-7.91um.tifxyz")


def main(level: int, half: int, slab: float, out: Path):
    s = tifxyz.read(Path(SEG))
    P, V = s.points, s.valid
    m = np.load("out/routeA/blob_full.npy")
    big = np.zeros(V.shape, bool)
    vv, uu = np.nonzero(m)
    big[np.minimum(vv * 3, V.shape[0] - 1),
        np.minimum(uu * 3, V.shape[1] - 1)] = True

    # Pick the z that carries the most flagged points, then centre on where
    # those flagged points actually are. The blob spans most of the scroll, so
    # its overall median lands in empty space -- an earlier version of this
    # script cropped a window containing zero flagged points.
    fp = P[big & V]
    zs = fp[:, 2]
    edges = np.arange(zs.min(), zs.max() + 2 * slab, 2 * slab)
    hist, _ = np.histogram(zs, bins=edges)
    z0 = float(edges[int(np.argmax(hist))] + slab)
    near = fp[np.abs(zs - z0) < slab]
    # Centre on an actual flagged point, never on the centroid: a scroll
    # cross-section is a spiral, so the mean/median of points along it lands in
    # the empty middle. Two earlier runs cropped the hole and found no trace.
    d2 = ((near[:, None, 0] - near[None, :, 0]) ** 2
          + (near[:, None, 1] - near[None, :, 1]) ** 2)
    dense = int(np.argmax((d2 < 300.0 ** 2).sum(1)))
    cx, cy = float(near[dense, 0]), float(near[dense, 1])
    print(f"densest flagged z={z0:.0f} holds {len(near)} flagged points; "
          f"centred on flagged point x {cx:.0f} y {cy:.0f} "
          f"({int((d2[dense] < 300.0 ** 2).sum())} flagged within 300 vx)")

    sel = (np.abs(P[:, :, 2] - z0) < slab) & V
    pts, fl = P[sel], big[sel]
    print(f"trace points in +-{slab} vx slab: {sel.sum()}, flagged {fl.sum()}")

    sc = 2 ** level
    x0, x1 = int(cx - half), int(cx + half)
    y0, y1 = int(cy - half), int(cy + half)
    fs = s3fs.S3FileSystem(anon=True)
    vol = zarr.open(s3fs.S3Map(root=f"{BUCKET}/PHerc0172/volumes/{VOLUME}",
                               s3=fs, check=False), mode="r")[str(level)]
    print(f"reading level {level} {vol.shape} ...")
    img = vol[int(round(z0 / sc)),
              max(0, y0 // sc):y1 // sc,
              max(0, x0 // sc):x1 // sc]
    print("slice", img.shape)

    keep = ((pts[:, 0] > x0) & (pts[:, 0] < x1)
            & (pts[:, 1] > y0) & (pts[:, 1] < y1))
    pts, fl = pts[keep], fl[keep]
    print(f"points inside crop: {len(pts)}, flagged {fl.sum()}")

    fig, ax = plt.subplots(1, 2, figsize=(19, 9.5))
    for k in (0, 1):
        ax[k].imshow(img, cmap="gray", origin="upper",
                     extent=[x0, x1, y1, y0])
        ax[k].set_xlabel("x (vx)")
    ax[0].set_title(f"scroll volume, axial z={z0:.0f}, level {level} "
                    f"({sc * 7.91:.1f} um/px)")
    ax[0].set_ylabel("y (vx)")
    ax[1].scatter(pts[~fl, 0], pts[~fl, 1], s=1.4, c="#33bbff",
                  label="trace, not flagged", linewidths=0)
    ax[1].scatter(pts[fl, 0], pts[fl, 1], s=1.4, c="#ff3333",
                  label="trace, flagged (<6 vx self-gap)", linewidths=0)
    ax[1].legend(loc="upper right", framealpha=0.85, markerscale=6)
    ax[1].set_title("same slice with the trace drawn on it")
    plt.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, default=1)
    ap.add_argument("--half", type=int, default=600)
    ap.add_argument("--slab", type=float, default=6.0)
    ap.add_argument("--out", type=Path, default=Path("out/routeA/inspect.png"))
    a = ap.parse_args()
    main(a.level, a.half, a.slab, a.out)
