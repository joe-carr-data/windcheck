"""Is the spectrum's shape an artifact of the diagonal choice?

The two tessellations (d0/d1) produce different crossing-pair sets, so events
do not correspond one-to-one; the comparison is at the levels that are
well-defined under both: pooled survival fractions, per-segment maxima and
medians, and the bimodality fractions themselves.

    uv run python bench/spectrum_compare.py \
        --a out/spectrum_full_d0.json --b out/spectrum_full_d1.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

THRESH = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0)


def seps(rows):
    return np.array([e["separation_mm"] for r in rows for e in r["events"]
                     if e["separation_mm"] is not None])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=Path, default=Path("out/spectrum_full_d0.json"))
    ap.add_argument("--b", type=Path, default=Path("out/spectrum_full_d1.json"))
    p = ap.parse_args()
    A = json.loads(p.a.read_text())
    B = json.loads(p.b.read_text())
    sa, sb = seps(A), seps(B)

    print(f"{p.a.name}: {len(A)} segments, {len(sa)} finite events")
    print(f"{p.b.name}: {len(B)} segments, {len(sb)} finite events\n")

    print("pooled survival P(sep >= s):        d0        d1     delta")
    for s in THRESH:
        fa, fb = float((sa >= s).mean()), float((sb >= s).mean())
        print(f"   {s:6.1f} mm                 {fa:9.4f} {fb:9.4f}  {fb-fa:+8.4f}")

    ia = {r["segment"]: r for r in A}
    matched = [(ia[r["segment"]], r) for r in B if r["segment"] in ia]
    dmax, dmed = [], []
    for ra, rb in matched:
        if ra["sep_mm_max"] and rb["sep_mm_max"]:
            dmax.append(rb["sep_mm_max"] / ra["sep_mm_max"])
        if ra["sep_mm_median"] and rb["sep_mm_median"]:
            dmed.append(rb["sep_mm_median"] / ra["sep_mm_median"])
    dmax, dmed = np.array(dmax), np.array(dmed)
    print(f"\nper-segment ratio d1/d0 over {len(matched)} shared segments:")
    print(f"   max separation: median {np.median(dmax):.3f}  "
          f"IQR [{np.percentile(dmax,25):.3f}, {np.percentile(dmax,75):.3f}]  "
          f"range [{dmax.min():.3f}, {dmax.max():.3f}]")
    print(f"   median sep:     median {np.median(dmed):.3f}  "
          f"IQR [{np.percentile(dmed,25):.3f}, {np.percentile(dmed,75):.3f}]  "
          f"range [{dmed.min():.3f}, {dmed.max():.3f}]")

    # the two decades that carry the bimodality claim
    for name, s in (("d0", sa), ("d1", sb)):
        mid = float(((s >= 0.3) & (s < 30)).mean())
        low = float((s < 0.1).mean())
        print(f"\n{name}: P(<0.1mm)={low:.4f}   P(0.3-30mm)={mid:.4f}   "
              f"P(>=30mm)={float((s >= 30).mean()):.4f}")

    inter_a = sum(r["events_inter_component"] for r in A)
    inter_b = sum(r["events_inter_component"] for r in B)
    print(f"\ninter-component events: d0 {inter_a}, d1 {inter_b}")


if __name__ == "__main__":
    main()
