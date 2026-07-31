"""The self-intersection spectrum of every censused segment.

Per segment: every crossing EVENT with its intrinsic set-to-set separation
along the surface, its component status, and its size. This is the object
round 11 said must exist before any claim about the shape of the intrinsic
distribution -- the per-segment maximum could not distinguish "continuous"
from "structure destroyed by taking a maximum".

Events are processed largest-first and capped per segment; the cap and how much
it dropped are recorded on the row, never silent.

    uv run python bench/spectrum_census.py --diagonal 0 --out out/spectrum_d0.json
"""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np

from windcheck import tifxyz
from windcheck.check import load_pairs
from windcheck.intrinsic import oriented_events, segment_spectrum

CORPORA = [
    ("Scroll 1", "data/scroll1_tifxyz", "20230205180739", "out/crossing_s1"),
    ("Scroll 5", "data/scroll5_tifxyz", "20241024131839", "out/crossing"),
    ("PHerc0139", "data/PHerc0139_tifxyz", "20250728140407", "out/crossing_0139"),
    ("PHerc0814", "data/PHerc0814_tifxyz", "20250804134230", "out/crossing_0814"),
    ("PHerc1667", "data/PHerc1667_tifxyz", "20231117161658", "out/crossing_1667"),
]
RES = re.compile(r"-(\d+\.?\d*)um\.tifxyz$")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--diagonal", type=int, default=0, choices=(0, 1))
    ap.add_argument("--max-events", type=int, default=150,
                    help="0 = uncapped (needs the early-exit engine)")
    ap.add_argument("--out", type=Path, default=Path("out/spectrum_d0.json"))
    ap.add_argument("--resume", action="store_true",
                    help="skip segments already in the sidecar .jsonl")
    ap.add_argument("--shard", default="",
                    help="'i/N': process only segments with index%%N==i; "
                         "checkpoints go to a per-shard jsonl, and every "
                         "shard's jsonl (plus the shared one) counts as done")
    a = ap.parse_args()
    shard_i, shard_n = (map(int, a.shard.split("/")) if a.shard else (0, 1))

    # Checkpoint sidecar: one row per line, appended as each segment lands,
    # so a killed run costs nothing. The final JSON is assembled from it.
    base = a.out.with_suffix(".jsonl")
    jsonl = (base if not a.shard
             else a.out.with_suffix(f".s{shard_i}of{shard_n}.jsonl"))
    jsonl.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    if a.resume:
        siblings = {base, jsonl} | set(base.parent.glob(
            base.stem + ".s*.jsonl"))
        for f in sorted(siblings):
            if f.exists():
                rows += [json.loads(x) for x in f.read_text().splitlines() if x]
        print(f"resuming: {len(rows)} segments already done", flush=True)
        if not jsonl.exists():
            jsonl.write_text("")
    else:
        jsonl.write_text("")
    done = {(r["corpus"], r["segment"]) for r in rows}
    seg_index = 0
    t0 = time.time()
    for corpus, root, volume, work in CORPORA:
        rootp = Path(root)
        if not rootp.exists():
            continue
        for d in sorted(p for p in rootp.iterdir() if p.is_dir()):
            seg_index += 1
            if (seg_index - 1) % shard_n != shard_i or (corpus, d.name) in done:
                continue
            m = sorted(d.glob(f"mesh/*{volume}*.tifxyz"))
            if not m or not (m[0] / "x.tif").exists():
                continue
            res = RES.search(m[0].name)
            if not res:
                continue
            vx = float(res.group(1))
            s = tifxyz.read(m[0])
            rec = load_pairs(Path(work) / f"{d.name[:40]}_d{a.diagonal}.csv")
            nv, nu = s.points.shape[:2]
            if len(rec):
                keep = ((rec["v1"] < nv - 1) & (rec["v2"] < nv - 1)
                        & (rec["u1"] < nu - 1) & (rec["u2"] < nu - 1))
                rec = rec[keep]
            t = time.time()
            spec = segment_spectrum(s.points, s.valid, rec, vx,
                                    diagonal=a.diagonal,
                                    max_events=(a.max_events or None))
            n_groups = 0 if len(rec) == 0 else len(oriented_events(rec))
            n_amb = sum(1 for e in spec if e.get("ambiguous"))
            seps = [e["separation_mm"] for e in spec
                    if e["separation_mm"] is not None]
            inter = sum(1 for e in spec if e["same_component"] is False)
            rows.append({
                "corpus": corpus, "segment": d.name, "voxel_um": vx,
                "pairs": int(len(rec)), "events_total": int(n_groups),
                "events_measured": len(spec), "events_inter_component": inter,
                "events_ambiguous": n_amb,
                "sep_mm_max": (round(max(seps), 3) if seps else None),
                "sep_mm_median": (round(float(np.median(seps)), 3)
                                  if seps else None),
                "events": spec,
            })
            with jsonl.open("a") as f:
                f.write(json.dumps(rows[-1]) + "\n")
            print(f"{corpus:10s} {d.name[:36]:38s} ev {len(spec):4d}/{n_groups:4d} "
                  f"inter {inter:3d}  max "
                  f"{('%9.2f' % max(seps)) if seps else '        -'} mm "
                  f"[{time.time() - t:5.1f}s]", flush=True)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    if a.shard:
        print(f"shard {a.shard} done; merge all *.jsonl into {a.out} "
              "when every shard has finished")
    else:
        a.out.write_text(json.dumps(rows, indent=1))
    ev = [e for r in rows for e in r["events"]]
    seps = [e["separation_mm"] for e in ev if e["separation_mm"] is not None]
    print(f"\n{len(rows)} segments, {len(ev)} events measured "
          f"({sum(1 for e in ev if e['same_component'] is False)} inter-component)"
          f" in {(time.time() - t0) / 60:.1f} min")
    if seps:
        q = np.percentile(seps, [50, 90, 99])
        print(f"pooled separations mm: median {q[0]:.2f}  p90 {q[1]:.2f}  "
              f"p99 {q[2]:.2f}  max {max(seps):.2f}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
