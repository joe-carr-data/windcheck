"""Regenerate docs/img/triage_ranking.png from the published sample outputs.

The figure ranks every multi-wrap trace by the size of its largest contiguous
flagged region. Ranking is all it does: the labelled and auto-grown populations
overlap at the low end, so a high bar means "open this one first", not "this
trace is wrong".

    uv run python bench/make_triage_figure.py

Reads sample_outputs/*_report.json, so the figure can never drift from the
numbers quoted in the README and the paper.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "img" / "triage_ranking.png"


def short(sample: str, name: str) -> str:
    tag = sample.replace("PHerc", "")
    if "auto_grown_" in name:
        rest = name.split("auto_grown_")[-1].replace("_flatboi", "")
        return f"{tag}  auto_grown {rest}"
    return f"{tag}  {name[:28]}"


def main() -> None:
    rows = []
    for sample in ("PHerc0172", "PHerc0814"):
        p = ROOT / "sample_outputs" / f"{sample}_report.json"
        if not p.exists():
            continue
        for r in json.loads(p.read_text()):
            # multi-wrap traces only; single windings have no previous wrap
            if "auto_grown" not in r["name"] and "auto_trace" not in r["name"]:
                continue
            rows.append((short(sample, r["name"]), r["blob_fraction"] * 100.0,
                         bool(r["valid"])))
    if not rows:
        raise SystemExit("no sample_outputs found; run the selfgap command first")

    rows.sort(key=lambda t: -t[1])
    labels = [f"{n}  (rejected)" if not ok else n for n, _, ok in rows]
    vals = [v for _, v, _ in rows]
    colors = ["#c0392b" if i == 0 and rows[0][2] else
              ("#cccccc" if not ok else "#2b6cb0")
              for i, (_, _, ok) in enumerate(rows)]

    fig, ax = plt.subplots(figsize=(10.5, 0.32 * len(rows) + 1.6))
    y = range(len(rows))
    ax.barh(list(y), vals, color=colors)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    for i, v in enumerate(vals):
        ax.text(v + max(vals) * 0.012, i, f"{v:.3f}", va="center", fontsize=8.5)
    ax.set_xlim(0, max(vals) * 1.16)
    ax.set_xlabel("largest contiguous flagged region, % of valid queries "
                  "submitted  (blob%)")
    ax.set_title("Triage order: which traces a human should look at first\n"
                 "grey = rejected by the validity filter, not a finding",
                 loc="left", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT, dpi=130)
    print(f"wrote {OUT}  ({len(rows)} traces, top {vals[0]:.3f}%)")


if __name__ == "__main__":
    main()
