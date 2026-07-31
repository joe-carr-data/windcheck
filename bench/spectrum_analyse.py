"""Does the event-level intrinsic distribution have structure?

The round-11 question. The per-segment maximum could not distinguish a
continuous distribution from structure destroyed by aggregation; this reads the
full event spectrum and answers with the distribution itself.

    uv run python bench/spectrum_analyse.py --json out/spectrum_d0.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", type=Path, default=Path("out/spectrum_d0.json"))
    a = ap.parse_args()
    rows = json.loads(a.json.read_text())

    ev = [dict(e, corpus=r["corpus"], segment=r["segment"])
          for r in rows for e in r["events"]]
    seps = np.array([e["separation_mm"] for e in ev
                     if e["separation_mm"] is not None])
    inter = [e for e in ev if e["same_component"] is False]
    capped = [r for r in rows if r["events_measured"] < r["events_total"]]

    print(f"{len(rows)} segments, {len(ev)} events measured, "
          f"{len(inter)} inter-component")
    print(f"capped segments: {len(capped)} "
          f"(events dropped: {sum(r['events_total'] - r['events_measured'] for r in capped)})"
          " -- caps keep the LARGEST events, so upper tail is complete\n")

    # The distribution itself, log-spaced survival
    print("pooled survival  P(separation >= s):")
    for s in (0.1, 0.3, 1, 3, 10, 30, 100, 300):
        f = float((seps >= s).mean())
        print(f"   {s:6.1f} mm   {f:７.4f}   n={int((seps >= s).sum())}"
              .replace("７", ""))

    # largest multiplicative gaps in the EVENT distribution
    v = np.sort(seps[seps > 0.05])
    g = v[1:] / np.maximum(v[:-1], 1e-9)
    wide = (v[1:] - v[:-1]) >= 0.5
    print("\nlargest event-level gaps (>=0.5 mm wide):")
    order = np.argsort(-(g * wide))[:5]
    for i in order:
        if wide[i]:
            print(f"   {v[i]:9.2f} -> {v[i+1]:9.2f} mm   {g[i]:5.2f}x")

    # per-segment bimodality: fraction of segments whose events span >30x
    span = []
    for r in rows:
        s = [e["separation_mm"] for e in r["events"]
             if e["separation_mm"] and e["separation_mm"] > 0.05]
        if len(s) >= 5:
            span.append((max(s) / min(s), r["segment"], len(s)))
    span.sort(reverse=True)
    wide_frac = sum(1 for x in span if x[0] > 30) / max(len(span), 1)
    print(f"\nsegments with >=5 events: {len(span)}; "
          f"fraction whose events span >30x: {wide_frac:.2f}")

    print("\ninter-component events by segment:")
    from collections import Counter
    for (c, s), n in Counter((e["corpus"], e["segment"]) for e in inter).most_common(8):
        print(f"   {c:10s} {s[:40]:42s} {n:4d}")


if __name__ == "__main__":
    main()
