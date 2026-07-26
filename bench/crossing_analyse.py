"""Turn raw crossing pairs into the units a reader can actually interpret.

Three things the raw CSV cannot tell you:

**Pairs are not defects.** One crossing patch, where a sheet passes through
itself over some area, produces thousands of intersecting triangle pairs. Quoting
the pair count makes a single event look like a catastrophe. Participating quads
are clustered by grid adjacency into regions, and a crossing EVENT is a pair of
regions. That is the unit to report.

**"Nonlocal" is a threshold, not a fact.** Existence of an intersection needs no
threshold; calling it wrap-scale does. So the separation distribution is reported
whole, as a survival curve, rather than collapsed to a count past one cutoff.

**Diagonal agreement is setwise.** Equal counts under the two quad
triangulations do not show the same regions crossed; two disjoint sets can have
identical sizes. Events are matched by canonical quad-pair key.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def load(csv: Path, verdict: str = "transverse") -> np.ndarray:
    if not csv.exists():
        return np.empty(0, dtype=[("v1", "i4"), ("u1", "i4"), ("v2", "i4"),
                                  ("u2", "i4"), ("pen", "f8"), ("ang", "f8")])
    rows = []
    with csv.open() as fh:
        header = fh.readline().strip().split(",")
        has_margin = "penetration" in header
        for line in fh:
            p = line.rstrip("\n").split(",")
            if p[4] != verdict:
                continue
            rows.append((int(p[0]), int(p[1]), int(p[2]), int(p[3]),
                         float(p[5]) if has_margin else np.nan,
                         float(p[6]) if has_margin else np.nan))
    return np.array(rows, dtype=[("v1", "i4"), ("u1", "i4"), ("v2", "i4"),
                                 ("u2", "i4"), ("pen", "f8"), ("ang", "f8")])


class DSU:
    def __init__(self):
        self.p: dict = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def components(rec: np.ndarray) -> tuple[dict, int]:
    """Cluster participating quads into grid-connected regions."""
    quads = set()
    for r in rec:
        quads.add((int(r["v1"]), int(r["u1"])))
        quads.add((int(r["v2"]), int(r["u2"])))
    d = DSU()
    for q in quads:
        d.find(q)
    for (v, u) in quads:
        for dv in (-1, 0, 1):
            for du in (-1, 0, 1):
                if dv == 0 and du == 0:
                    continue
                n = (v + dv, u + du)
                if n in quads:
                    d.union((v, u), n)
    label = {q: d.find(q) for q in quads}
    return label, len(set(label.values()))


def events(rec: np.ndarray, label: dict) -> set:
    """A crossing event is an unordered pair of distinct regions."""
    ev = set()
    for r in rec:
        a = label[(int(r["v1"]), int(r["u1"]))]
        b = label[(int(r["v2"]), int(r["u2"]))]
        ev.add((a, b) if a <= b else (b, a))
    return ev


def canonical(rec: np.ndarray) -> set:
    out = set()
    for r in rec:
        a = (int(r["v1"]), int(r["u1"]))
        b = (int(r["v2"]), int(r["u2"]))
        out.add((a, b) if a <= b else (b, a))
    return out


def summarise(seg: str, d0: Path, d1: Path) -> dict:
    r0, r1 = load(d0), load(d1)
    if len(r0) == 0:
        return {"segment": seg, "pairs": 0, "regions": 0, "events": 0,
                "events_both_diagonals": 0, "max_separation": 0}
    sep0 = np.maximum(np.abs(r0["v1"] - r0["v2"]), np.abs(r0["u1"] - r0["u2"]))
    lab0, nreg = components(r0)
    ev0 = events(r0, lab0)

    # setwise, not by count
    both = len(canonical(r0) & canonical(r1)) if len(r1) else 0

    surv = {int(t): int((sep0 > t).sum()) for t in (0, 50, 100, 200, 300, 500,
                                                   750, 1000, 1500)}
    return {
        "segment": seg,
        "pairs": int(len(r0)),
        "regions": nreg,
        "events": len(ev0),
        "quad_pairs_shared_with_other_diagonal": both,
        "quad_pairs_d0": len(canonical(r0)),
        "quad_pairs_d1": len(canonical(r1)) if len(r1) else 0,
        "max_separation": int(sep0.max()),
        "median_separation": float(np.median(sep0)),
        "separation_survival": surv,
        "penetration_median": float(np.nanmedian(r0["pen"])),
        "penetration_min": float(np.nanmin(r0["pen"])),
        "angle_median_deg": float(np.nanmedian(r0["ang"])),
    }


def corpus_filter(root: Path | None):
    """Accept only CSVs belonging to the segments of one corpus.

    Runs over different scrolls wrote into a shared output directory before
    per-run work dirs existed, so a bare glob picks up PHerc0139 and PHerc1667
    files alongside Scroll 5's and silently reports a mixture as one corpus.
    Naming the root makes the population explicit at the point of use, which is
    where it has to be right.
    """
    if root is None:
        return lambda stem: True
    known = {d.name[:40] for d in root.iterdir() if d.is_dir()}
    return lambda stem: stem in known


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, default=Path("out/crossing"))
    ap.add_argument("--root", type=Path, default=Path("data/scroll5_tifxyz"),
                    help="corpus whose segments to include; pass an empty "
                         "string to take every CSV in --dir")
    ap.add_argument("--json", type=Path, default=Path("out/crossing/analysis.json"))
    a = ap.parse_args()

    keep = corpus_filter(a.root if a.root and str(a.root) else None)
    rows = []
    for d0 in sorted(a.dir.glob("*_d0.csv")):
        seg = d0.name[:-7]
        if not keep(seg):
            continue
        rows.append(summarise(seg, d0, a.dir / f"{seg}_d1.csv"))

    lab = [r for r in rows if "-w0" in r["segment"]]
    ag = [r for r in rows if "auto_grown" in r["segment"]
          or "auto_trace" in r["segment"]]

    print(f"{'segment':40s} {'pairs':>9s} {'regions':>8s} {'events':>7s} "
          f"{'maxsep':>7s} {'pen med':>8s}")
    print("-" * 84)
    for r in sorted(rows, key=lambda x: -x["events"])[:24]:
        print(f"{r['segment'][:38]:40s} {r['pairs']:9,d} {r['regions']:8,d} "
              f"{r['events']:7,d} {r['max_separation']:7,d} "
              f"{r['penetration_median']:8.3f}")

    print(f"\nlabelled single-wrap ({len(lab)}): "
          f"{sum(r['events'] for r in lab):,} events, "
          f"max separation {max((r['max_separation'] for r in lab), default=0)}")
    print(f"multi-wrap ({len(ag)}): {sum(r['events'] for r in ag):,} events, "
          f"max separation {max((r['max_separation'] for r in ag), default=0)}")

    a.json.write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
