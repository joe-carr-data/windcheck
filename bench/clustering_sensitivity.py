"""How much does the event count depend on the clustering rule?

Two defensible definitions:

  region-pair    cluster participating quads by grid adjacency, then call an
                 unordered pair of regions one event
  product-space  connect two crossing records only when BOTH endpoints are
                 adjacent (allowing the endpoint order to swap), which tracks
                 continuity of the contact relation rather than connectivity of
                 either surface patch alone

The second is the better definition. Reported together so the reader can see the
size of the choice rather than take one number on trust.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from crossing_analyse import DSU, components, corpus_filter, events, load


def product_events(rec) -> int:
    recs = [((int(r["v1"]), int(r["u1"])), (int(r["v2"]), int(r["u2"])))
            for r in rec]
    idx = {}
    for i, (a, b) in enumerate(recs):
        for key in ((a, b), (b, a)):
            idx.setdefault(key, []).append(i)
    d = DSU()
    for i in range(len(recs)):
        d.find(i)
    for i, (a, b) in enumerate(recs):
        for da in (-1, 0, 1):
            for du in (-1, 0, 1):
                for db in (-1, 0, 1):
                    for dv in (-1, 0, 1):
                        k = ((a[0] + da, a[1] + du), (b[0] + db, b[1] + dv))
                        for j in idx.get(k, ()):
                            d.union(i, j)
    return len({d.find(i) for i in range(len(recs))})


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, default=Path("out/crossing"))
    ap.add_argument("--root", type=Path,
                    default=Path("data/scroll5_tifxyz"))
    args = ap.parse_args()
    keep = corpus_filter(args.root if str(args.root) else None)

    cut = 200
    CAP = 40000            # the O(n^2) product-space union is what this bounds
    rows = []
    skipped, empty = [], 0
    for f in sorted(args.dir.glob("*_d0.csv")):
        seg = f.name[:-7]
        if not keep(seg):
            continue
        rec = load(f)
        if len(rec) == 0:
            empty += 1
            continue
        if len(rec) > CAP:
            # Counted and named, not dropped in silence. A sensitivity analysis
            # that quietly excludes its largest cases reads as covering
            # everything while covering only the easy half, and the excluded
            # traces here are exactly the wrap-scale ones the result is about.
            skipped.append((seg, len(rec)))
            continue
        sep = np.maximum(np.abs(rec["v1"] - rec["v2"]),
                         np.abs(rec["u1"] - rec["u2"]))
        for tag, sub in (("local", rec[sep <= cut]), ("nonlocal", rec[sep > cut])):
            if len(sub) == 0:
                continue
            lab, _ = components(sub)
            rows.append((seg, tag, len(events(sub, lab)), product_events(sub)))
    print(f"{'class':10s} {'n segs':>7s} {'region-pair':>12s} {'product-space':>14s} {'ratio':>7s}")
    print("-" * 56)
    for tag in ("local", "nonlocal"):
        r = [x for x in rows if x[1] == tag]
        if not r:
            continue
        a = sum(x[2] for x in r); b = sum(x[3] for x in r)
        print(f"{tag:10s} {len(r):7d} {a:12,d} {b:14,d} {b/max(a,1):7.2f}x")
    print(f"\nScope: {empty} file(s) had no crossings; {len(skipped)} exceeded "
          f"the {CAP:,}-pair cap on the product-space union and were NOT "
          f"analysed.")
    for seg, n in sorted(skipped, key=lambda x: -x[1]):
        print(f"    not analysed: {seg[:44]:46s} {n:9,d} pairs")
    print("Every ratio above is therefore 'on the analysed subset'.")
    print("\nThe SEGMENT-level result is unaffected: which traces have wrap-scale")
    print("overlap does not depend on how events within them are grouped.")
