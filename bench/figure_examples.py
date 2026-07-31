"""Render one wrap-scale and one local self-overlap in the CT volume.

Numbers do not make this legible. An axial slice with the trace drawn on it does:
the sheets are visible as bright arcs, and where the trace passes through itself
you can see the same trace occupying one arc twice, coloured by how far apart the
two visits are along its own parameter.

Left panel is the volume alone so the reader can see the sheets without our
interpretation painted over them. Right panel adds the trace.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import s3fs
import zarr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crossing_analyse import load                             # noqa: E402

from windcheck import tifxyz                                  # noqa: E402

BUCKET = "vesuvius-challenge-open-data"
VOLUME = "20241024131839-7.910um-53keV-masked.zarr"

CASES = [
    ("wrap-scale overlap",
     "20251115002745-auto_grown_20251115002740308_5_flatboi", 300, None),
    ("local overlap",
     "20251109232817-w065_20251109232817724_flatboi", None, 30),
]


def pick(seg: str, gt: int | None, le: int | None):
    rec = load(Path("out/crossing") / f"{seg[:40]}_d0.csv")
    sep = np.maximum(np.abs(rec["v1"] - rec["v2"]),
                     np.abs(rec["u1"] - rec["u2"]))
    sub = rec[sep > gt] if gt else rec[sep <= le]
    r = sub[int(np.argmax(sub["pen"]))]
    return ((int(r["v1"]), int(r["u1"])), (int(r["v2"]), int(r["u2"])),
            float(r["pen"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--half", type=int, default=260)
    ap.add_argument("--slab", type=float, default=7.0)
    ap.add_argument("--level", type=int, default=0)
    ap.add_argument("--out", type=Path,
                    default=Path("out/routeA/examples.png"))
    a = ap.parse_args()

    fs = s3fs.S3FileSystem(anon=True)
    vol = zarr.open(s3fs.S3Map(root=f"{BUCKET}/PHerc0172/volumes/{VOLUME}",
                               s3=fs, check=False), mode="r")[str(a.level)]
    sc = 2 ** a.level

    fig, axes = plt.subplots(2, 2, figsize=(17, 16))
    for row, (title, seg, gt, le) in enumerate(CASES):
        q1, q2, pen = pick(seg, gt, le)
        s = tifxyz.read(sorted((Path("data/scroll5_tifxyz") / seg)
                               .glob("mesh/*20241024131839*.tifxyz"))[0])
        P, V = s.points, s.valid
        c = P[q1]
        gsep = max(abs(q1[0] - q2[0]), abs(q1[1] - q2[1]))

        sel = (np.abs(P[:, :, 2] - c[2]) < a.slab) & V
        gi, gj = np.nonzero(sel)
        pts, ucoord = P[sel], gj.astype(float)
        keep = (np.abs(pts[:, 0] - c[0]) < a.half) & (np.abs(pts[:, 1] - c[1]) < a.half)
        pts, ucoord = pts[keep], ucoord[keep]

        x0, x1 = int(c[0] - a.half), int(c[0] + a.half)
        y0, y1 = int(c[1] - a.half), int(c[1] + a.half)
        img = vol[int(round(c[2] / sc)), y0 // sc:y1 // sc, x0 // sc:x1 // sc]
        print(f"{title}: quads {q1} vs {q2}, grid separation {gsep}, "
              f"{len(pts)} trace points in view", flush=True)

        for col in (0, 1):
            ax = axes[row, col]
            ax.imshow(img, cmap="gray", extent=[x0, x1, y1, y0])
            ax.set_xlabel("x (voxels)")
            if col == 0:
                ax.set_ylabel("y (voxels)")
                ax.set_title(f"{title} — CT only, z={c[2]:.0f}", loc="left")
            else:
                sct = ax.scatter(pts[:, 0], pts[:, 1], s=13, c=ucoord,
                                 cmap="autumn", linewidths=0)
                # The two marks land almost on top of each other -- that IS the
                # finding -- so the parameter distance has to be written out or
                # the reader cannot see that these are two different visits.
                for q, mk, dy in ((q1, "o", -34), (q2, "s", 26)):
                    ax.scatter([P[q][0]], [P[q][1]], s=210, marker=mk,
                               facecolors="none", edgecolors="#00e5ff", lw=2.0)
                    ax.annotate(f"u = {q[1]}", (P[q][0], P[q][1]),
                                textcoords="offset points", xytext=(16, dy),
                                color="#00e5ff", fontsize=11, weight="bold",
                                arrowprops=dict(arrowstyle="-", color="#00e5ff",
                                                lw=1.0, alpha=0.8))
                plt.colorbar(sct, ax=ax, label="u = position along trace",
                             fraction=0.046)
                ax.set_title(f"same trace, {gsep} grid cells apart, "
                             f"penetration {pen:.1f} vx", loc="left")
    fig.suptitle("The trace passing through itself: one wrap-scale case, "
                 "one local case\ncircle and square mark the two crossing quads",
                 fontsize=13)
    fig.tight_layout()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=105)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
