"""Run `windcheck check` over every segment we hold, and bundle the results.

The point is that nobody should have to download 18 GB to find out whether their
own trace self-overlaps. This writes one certificate per segment plus a single
index, small enough to commit, so the answer is a file lookup.

Volumes are named per corpus and must match the one the census used, because the
crossing indices are indices into that grid -- reading a different resolution of
the same segment silently reinterprets them, which has happened here once and
produced an entirely plausible fiction.

    uv run python bench/precompute_all.py --out out/bundle
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from windcheck import check as C
from windcheck import classify

CORPORA = [
    ("Scroll 1 (PHercParis4)", "data/scroll1_tifxyz", "20230205180739"),
    ("Scroll 5 (PHerc0172)", "data/scroll5_tifxyz", "20241024131839"),
    ("PHerc0139", "data/PHerc0139_tifxyz", ""),
    ("PHerc0814", "data/PHerc0814_tifxyz", "20250804134230"),
    ("PHerc1667", "data/PHerc1667_tifxyz", "20231117161658"),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("out/bundle"))
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    certs = a.out / "certificates"
    certs.mkdir(exist_ok=True)

    index = []
    for corpus, root, volume in CORPORA:
        rootp = Path(root)
        if not rootp.exists():
            print(f"skip {corpus}: {root} not present")
            continue
        for d in sorted(p for p in rootp.iterdir() if p.is_dir()):
            try:
                r = C.analyse(d, certs, volume=volume, threads=0)
            except Exception as e:                      # noqa: BLE001
                print(f"  {d.name[:40]:42s} FAILED {type(e).__name__}: {e}")
                continue
            if r is None:
                continue
            index.append({
                "corpus": corpus, "segment": r["name"],
                "grid": r["grid"], "triangles": r["triangles"],
                "pairs": r["pairs"], "events": r["events"],
                "events_beyond_cut": r["events_beyond_cut"],
                "covering_span_rev": (round(r["span_rev"], 3)
                                      if r["span_rev"] else None),
                "separation_rev": (round(r["sep_rev"], 3)
                                   if r["sep_rev"] else None),
                "band": r["band"],
                "crossing_status": classify.crossing_status(r["pairs"]),
                "period_status": r["period_status"],
                "verdict": r["verdict"],
            })
            print(f"  {corpus[:12]:14s} {r['name'][:38]:40s} "
                  f"{r['band']:16s} period {r['period_status']:12s} "
                  f"sep {str(round(r['sep_rev'], 3)) if r['sep_rev'] else '-':>7s}",
                  flush=True)

    (a.out / "index.json").write_text(json.dumps(index, indent=1))
    # Intermediates are large and are not part of the bundle: the per-pair CSVs
    # and the engine's serialised atlas/query/result buffers.
    for pat in ("*_pairs.csv", "_atlas.bin", "_np_atlas.bin", "_np_query.bin",
                "_np_result.bin"):
        for f in certs.glob(pat):
            f.unlink()

    n = len(index)
    from collections import Counter
    print(f"\n{n} segments -> {a.out}/index.json")
    for label, key in (("separation band", "band"),
                       ("period status", "period_status"),
                       ("crossing status", "crossing_status")):
        print(f"  {label}:")
        for k, v in sorted(Counter(r[key] for r in index).items()):
            print(f"     {k:22s} {v:4d}")


if __name__ == "__main__":
    main()
