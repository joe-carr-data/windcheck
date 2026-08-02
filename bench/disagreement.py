"""Does inter-patch disagreement predict error where no crossing exists?

Our census has one structural blind spot, stated publicly since July: a
trace that leaves the correct sheet and never returns need not
self-intersect, so we cannot see it. Every error we detect folds back.

The patch set is 84,316 independent reconstructions of one scroll. Where
two overlap they disagree about where the sheet is, and that disagreement
is measurable EVERYWHERE they overlap, not only where they cross. If it
predicts self-intersection in published traces of the same region, then
it is a quality signal that works where the census is blind.

Pre-registered in notes/PREREG-DISAGREEMENT.md before any statistic was
computed against any label. The gate is AUC > 0.75 to stand, <= 0.65 to
die. Read that file before changing anything here.

Stages:

    control    the statistic must recover a known 0.25 vx displacement
    field      per-pair disagreement, from patches only, no labels
    evaluate   AUC against census labels, with the confound checks
"""
from __future__ import annotations

import argparse
import json
import concurrent.futures as cf
import random
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from windcheck import tifxyz  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from patch_audit import download, have_patch  # noqa: E402

OUT = Path("out/disagree")


def disagreement(a: Path, b: Path, max_samples: int = 4000) -> dict | None:
    """Median distance from A's points inside B's box to B's surface.

    Median rather than mean: one bad corner must not carry the statistic.
    Returns None when the overlap holds too few points to be meaningful,
    and that refusal is recorded rather than silently treated as zero.
    """
    sa, sb = tifxyz.read(a), tifxyz.read(b)
    PA, PB = sa.points[sa.valid], sb.points[sb.valid]
    if len(PA) < 50 or len(PB) < 50:
        return None
    lo, hi = PB.min(axis=0), PB.max(axis=0)
    inside = np.all((PA >= lo) & (PA <= hi), axis=1)
    n_in = int(inside.sum())
    if n_in < 50:
        return None
    P = PA[inside]
    if len(P) > max_samples:
        idx = np.linspace(0, len(P) - 1, max_samples).astype(int)
        P = P[idx]
    d, _ = cKDTree(PB).query(P, k=1)
    return {"n_points": n_in, "median_vx": float(np.median(d)),
            "p90_vx": float(np.percentile(d, 90)),
            "mean_vx": float(d.mean()),
            "centroid": [float(x) for x in P.mean(axis=0)]}


def stage_control(data: Path) -> int:
    """A known displacement must come back out. Otherwise nothing else counts."""
    names = Path("out/patches/index.txt").read_text().split()
    ok = True
    for name in names[1000:1003]:
        if not (data / name).is_dir() and download(name, data) is None:
            continue
        s = tifxyz.read(data / name)
        for shift in (0.0, 0.25, 1.0):
            d = data / f"_ctrl_shift{shift}.tifxyz"
            if d.exists():
                shutil.rmtree(d)
            d.mkdir(parents=True)
            import tifffile
            for i, ax in enumerate("xyz"):
                arr = np.asarray(s.points[..., i], np.float32).copy()
                if i == 0:
                    arr = arr + np.float32(shift)
                arr[~s.valid] = -1.0
                tifffile.imwrite(d / f"{ax}.tif", arr)
            (d / "meta.json").write_text('{"format":"tifxyz","scale":[1,1]}')
            r = disagreement(data / name, d)
            shutil.rmtree(d, ignore_errors=True)
            if r is None:
                print(f"  {name[:34]} shift {shift}: REFUSED (too few points)")
                continue
            err = abs(r["median_vx"] - shift)
            good = err <= 0.05
            ok &= good
            print(f"  {name[:34]:36s} shift {shift:4.2f} -> "
                  f"measured {r['median_vx']:6.3f}  err {err:5.3f}  "
                  f"{'PASS' if good else 'FAIL'}")
    print(f"\ncontrol {'PASSED' if ok else 'FAILED'}")
    if not ok:
        print("the statistic cannot recover a known displacement; "
              "per the pre-registration the study stops here")
    return 0 if ok else 1


def stage_field(pairs_file: Path, out: Path, sample: int, seed: int,
                data: Path, jobs: int = 24) -> int:
    """Per-pair disagreement. No crossing information is read here."""
    rng = random.Random(seed)
    total, reservoir = 0, []
    with pairs_file.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            total += 1
            if len(reservoir) < sample:
                reservoir.append(json.loads(line))
            else:
                j = rng.randrange(total)
                if j < sample:
                    reservoir[j] = json.loads(line)
    print(f"streamed {total} pairs, sampled {len(reservoir)}")

    done = set()
    if out.exists():
        done = {(r["a"], r["b"]) for r in
                (json.loads(l) for l in out.read_text().splitlines() if l.strip())}
    todo = [p for p in reservoir if (p["a"], p["b"]) not in done]
    t0 = time.time()

    def one(p: dict) -> dict:
        rec = dict(p)
        try:
            fetched = all(have_patch(data / p[k]) or download(p[k], data)
                          for k in ("a", "b"))
            if not fetched:
                rec["status"] = "download_failed"
                return rec
            d = disagreement(data / p["a"], data / p["b"])
            if d is None:
                rec["status"] = "insufficient_overlap"
            else:
                rec.update(d)
                rec["status"] = "ok"
        except Exception as e:                              # noqa: BLE001
            # One unreadable pair must not end a run of tens of thousands.
            rec["status"] = f"error: {type(e).__name__}: {e}"
        return rec

    # Fetching two patches dominates each pair, so this is latency-bound in
    # exactly the way the census was. Single-threaded it ran at 0.5 pairs/s,
    # which is too slow to populate cubes densely enough to evaluate.
    with out.open("a") as sink, cf.ThreadPoolExecutor(jobs) as pool:
        for i, rec in enumerate(pool.map(one, todo), 1):
            sink.write(json.dumps(rec) + "\n")
            if i % 500 == 0:
                sink.flush()
                el = time.time() - t0
                print(f"  {i}/{len(todo)}  {i/el:.1f}/s  "
                      f"eta {(len(todo)-i)/(i/el)/60:.0f}m", flush=True)
    print(f"wrote {out}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["control", "field"])
    ap.add_argument("--pairs-file", default="out/patches/pairs.jsonl")
    ap.add_argument("--out", default=str(OUT / "field.jsonl"))
    ap.add_argument("--data", default="out/patches/pairdata")
    ap.add_argument("--sample", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--jobs", type=int, default=24)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    Path(a.data).mkdir(parents=True, exist_ok=True)
    if a.stage == "control":
        return stage_control(Path(a.data))
    return stage_field(Path(a.pairs_file), Path(a.out), a.sample, a.seed,
                       Path(a.data), a.jobs)


if __name__ == "__main__":
    raise SystemExit(main())
