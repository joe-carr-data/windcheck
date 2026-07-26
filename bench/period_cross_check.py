"""Two independent estimates of columns-per-revolution, compared.

The turning estimator in `revolution_diag.py` integrates how far the surface's
own centreline rotates. It is model-free but it needs a centre, and a segment
covering only part of a crushed spiral can put that centre outside the arc.

This second estimator uses no centre and no axis at all: for each grid point it
asks which other part of the SAME surface, at least `gap` columns away, is
physically nearest. On a rolled sheet that nearest part is the adjacent wrap,
so the column offset to it is one revolution. It is the same machinery the
self-gap code already uses (`selfgap.estimate_exclude_u`), reported raw instead
of scaled to a quarter.

Neither estimator is trusted alone. Where they agree the period is real; where
they disagree the disagreement is the finding, and gets looked at rather than
averaged away.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from revolution_diag import revolution_period                  # noqa: E402

from windcheck import atlas, tifxyz                            # noqa: E402


class _Entry:
    def __init__(self, p):
        self.path, self.name, self.winding = p, p.name, None


def period_from_neighbour(path: Path, work: Path, threads: int = 0,
                          stride: int = 4) -> tuple[float, int]:
    """Median column offset to the physically nearest non-local surface.

    Returns (period, n_supporting_samples). A period supported by few samples
    is not a period; the caller decides, so the count is returned rather than
    folded into a None.
    """
    s = tifxyz.read(path)
    v, u = np.nonzero(s.valid[::stride, ::stride])
    if len(v) < 2000:
        return float("nan"), 0
    V, U = v * stride, u * stride
    pts = s.points[V, U]
    ncol = s.shape[1]
    gap = max(ncol // 100, 3)

    work.mkdir(parents=True, exist_ok=True)
    atlas.write_atlas([_Entry(path)], work / "pc_atlas.bin")
    atlas.write_queries_grouped(pts, U, work / "pc_query.bin")
    r = atlas.run_engine(work / "pc_atlas.bin", work / "pc_query.bin",
                         work / "pc_result.bin", threads=threads, exclude_u=gap)
    d, w1 = r["d1"], r["w1"]
    # Keep only points whose nearest non-local partner is close enough to be a
    # neighbouring sheet rather than a fold elsewhere in the volume. 24 vx is
    # ~2x the 12-17 vx sheet spacing measured on this scroll.
    m = np.isfinite(d) & (d < 24.0)
    if m.sum() < 200:
        return float("nan"), int(m.sum())
    off = np.abs(w1[m] - U[m])
    return float(np.median(off)), int(m.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--volume", default="",
                    help="MUST be the volume the census ran on; the offsets "
                         "here are in that grid's columns")
    ap.add_argument("--work", type=Path, default=Path("out/period"))
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()

    rows = []
    print(f"{'segment':40s} {'turn':>9s} {'neigh':>9s} {'ratio':>7s} {'n':>7s}")
    print("-" * 76)
    for d in sorted(a.root.iterdir()):
        if not d.is_dir():
            continue
        m = sorted(d.glob(f"mesh/*{a.volume}*.tifxyz"))
        if not m:
            continue
        s = tifxyz.read(m[0])
        turn = revolution_period(s.points, s.valid)
        neigh, n = period_from_neighbour(m[0], a.work)
        ratio = turn / neigh if (turn == turn and neigh == neigh and neigh) else float("nan")
        rows.append({"segment": d.name, "turn": turn, "neigh": neigh,
                     "ratio": ratio, "n_support": n})
        print(f"{d.name[:38]:40s} {turn:9.1f} {neigh:9.1f} {ratio:7.3f} {n:7d}",
              flush=True)

    if a.json:
        a.json.write_text(json.dumps(rows, indent=2, default=float))
        print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
