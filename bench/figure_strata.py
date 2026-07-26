"""The whole distribution as a scatter, not two extreme values.

A gap between a maximum and a minimum is two numbers, and two numbers can be
produced by luck. This plots every segment, so a reader sees the shape and can
disagree with the reading.

The axes are chosen to remove a circularity in the headline table. There, both
the grouping variable (covering span) and the outcome (separation) are divided
by the same estimated period, so an error in that period moves a point along
the diagonal and could manufacture a trend. Here:

    x  =  u_span / estimated_period        how far round the scroll it reaches
    y  =  du_max / u_span                  separation as a fraction of the
                                           segment's OWN column extent

The outcome no longer contains the estimated period at all. What remains is the
exposure relationship, stated honestly: a longer surface has more room to
contain a widely separated pair, and the plot shows how much of the effect that
accounts for.

    uv run python bench/figure_strata.py --out out/fig_strata.png
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

CORPORA = [
    ("Scroll 5", "out/crossing/revdiag.json", "tab:blue", "o"),
    ("PHerc0814", "out/crossing_0814/revdiag.json", "tab:orange", "s"),
    ("PHerc0139", "out/crossing_0139/revdiag.json", "tab:green", "^"),
    ("PHerc1667", "out/crossing_1667/revdiag.json", "tab:red", "D"),
]
FLOOR = 3e-4          # so "no crossing at all" is visible on a log axis


def points(path: Path):
    x, y, clean = [], [], []
    for r in json.loads(path.read_text()):
        span, sep = r.get("span_per_rev"), r.get("sep_max_per_span")
        if span is None or sep is None or span < 0.2:
            continue
        x.append(span)
        y.append(sep if sep > 0 else FLOOR)
        clean.append(sep == 0)
    return np.array(x), np.array(y), np.array(clean, dtype=bool)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("out/fig_strata.png"))
    a = ap.parse_args()

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    total = 0
    for name, path, colour, marker in CORPORA:
        p = Path(path)
        if not p.exists():
            continue
        x, y, clean = points(p)
        total += len(x)
        ax.scatter(x[~clean], y[~clean], c=colour, marker=marker, s=34,
                   alpha=0.8, edgecolors="none", label=f"{name} ({len(x)})")
        # Segments with no crossing at all are drawn hollow on the floor line:
        # dropping them would quietly remove the cleanest evidence there is.
        if clean.any():
            ax.scatter(x[clean], y[clean], facecolors="none", edgecolors=colour,
                       marker=marker, s=34, linewidths=0.9)

    ax.axvline(1.0, color="0.55", lw=0.9, ls="--")
    ax.axvline(2.0, color="0.55", lw=0.9, ls="--")
    ax.text(1.0, 1.35, " covers one\n revolution", fontsize=8, color="0.35",
            va="top")
    ax.text(2.0, 1.35, " two", fontsize=8, color="0.35", va="top")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("covering span  =  column span / estimated revolution period")
    ax.set_ylabel("separation / own column span")
    ax.set_title("Self-overlap separation against how far a segment reaches\n"
                 f"{total} segments, four scrolls; hollow = no crossing found",
                 fontsize=11)
    ax.grid(True, which="major", lw=0.4, alpha=0.4)
    ax.legend(fontsize=8, loc="lower right", framealpha=0.9)
    fig.tight_layout()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=170)
    print(f"wrote {a.out}  ({total} segments)")


if __name__ == "__main__":
    main()
