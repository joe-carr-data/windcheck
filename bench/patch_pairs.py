"""Does patch cleanliness compose?

Every one of the 84,316 published patches was censused individually and
84,311 came back clean. That says nothing about what happens when they
are assembled, and assembly is the point of the dataset: two surfaces
that are each free of self-intersection can still pass through each
other, and any merge of such a pair either self-intersects or must throw
part of one away.

The engine censuses one surface at a time, so a pair is tested by
stitching both into a single grid separated by a band of invalid rows.
The band is far wider than the Chebyshev adjacency exclusion, so no
genuine contact between the two is suppressed, and since each patch is
independently clean, every contact the census reports is necessarily
BETWEEN them. That is not a merge -- there is no reparametrisation and
no seam -- but it answers the load-bearing question, which is whether
the two occupy the same space transversally.

Three stages, each resumable:

    bboxes   fetch meta.json for every patch (1 request each, not 4)
    pairs    index the boxes and enumerate overlapping candidates
    census   fetch, stitch and census a sample of those pairs
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import itertools
import json
import random
import shutil
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import tifffile  # noqa: E402

from windcheck import pipeline, tifxyz  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from patch_audit import BASE, PREFIX, download, get, sha256  # noqa: E402

INDEX = Path("out/patches/index.txt")
BOXES = Path("out/patches/bboxes.jsonl")

# Wider than the Chebyshev exclusion of 1 by a comfortable margin, so a
# real contact between the two patches can never be mistaken for adjacency.
GAP_ROWS = 8


# ------------------------------------------------------------------ bboxes

def stage_bboxes(names: list[str], jobs: int) -> None:
    done = set()
    if BOXES.exists():
        done = {json.loads(l)["patch"] for l in BOXES.read_text().splitlines()
                if l.strip()}
    todo = [n for n in names if n not in done]
    print(f"bboxes: {len(done)} cached, {len(todo)} to fetch")
    if not todo:
        return
    t0 = time.time()
    with BOXES.open("a") as sink, cf.ThreadPoolExecutor(jobs) as pool:
        def one(n: str) -> dict:
            body = get(f"{PREFIX}/{n}/meta.json")
            if body is None:
                return {"patch": n, "status": "fetch_failed"}
            try:
                m = json.loads(body)
                return {"patch": n, "status": "ok", "bbox": m["bbox"],
                        "area_vx2": m.get("area_vx2"),
                        "surface_points": m.get("surface_points")}
            except Exception as e:                          # noqa: BLE001
                return {"patch": n, "status": f"parse_failed: {e}"}

        for i, rec in enumerate(pool.map(one, todo), 1):
            sink.write(json.dumps(rec) + "\n")
            if i % 2000 == 0:
                el = time.time() - t0
                print(f"  {i}/{len(todo)}  {i/el:.0f}/s  "
                      f"eta {(len(todo)-i)/(i/el)/60:.0f}m", flush=True)
    print(f"bboxes done in {time.time()-t0:.0f}s")


# ------------------------------------------------------------------- pairs

def stage_pairs(out: Path, min_overlap_vx: float) -> None:
    recs = [json.loads(l) for l in BOXES.read_text().splitlines() if l.strip()]
    ok = [r for r in recs if r.get("status") == "ok"]
    print(f"boxes: {len(ok)} usable of {len(recs)}")

    lo = np.array([r["bbox"][0] for r in ok], float)
    hi = np.array([r["bbox"][1] for r in ok], float)
    names = [r["patch"] for r in ok]
    ext = hi - lo
    print(f"box extent, median per axis: {np.median(ext, axis=0).round(1)}")

    # Uniform grid over the median box size: a pair can only overlap if it
    # shares a cell, which turns an 84k^2 comparison into a local one.
    cell = float(np.median(ext.max(axis=1)))
    buckets: dict[tuple, list[int]] = defaultdict(list)
    for i in range(len(ok)):
        a = np.floor(lo[i] / cell).astype(int)
        b = np.floor(hi[i] / cell).astype(int)
        for ix in range(a[0], b[0] + 1):
            for iy in range(a[1], b[1] + 1):
                for iz in range(a[2], b[2] + 1):
                    buckets[(ix, iy, iz)].append(i)

    seen, pairs = set(), []
    for idxs in buckets.values():
        for i, j in itertools.combinations(sorted(idxs), 2):
            if (i, j) in seen:
                continue
            seen.add((i, j))
            o_lo = np.maximum(lo[i], lo[j])
            o_hi = np.minimum(hi[i], hi[j])
            d = o_hi - o_lo
            if np.all(d > 0):
                vol = float(np.prod(d))
                if vol >= min_overlap_vx:
                    pairs.append({"a": names[i], "b": names[j],
                                  "overlap_vx3": vol,
                                  "overlap_extent": [float(x) for x in d]})
    pairs.sort(key=lambda p: -p["overlap_vx3"])
    out.write_text("\n".join(json.dumps(p) for p in pairs) + "\n")
    print(f"overlapping pairs (bbox volume >= {min_overlap_vx:g}): {len(pairs)}")
    print(f"wrote {out}")


# ------------------------------------------------------------------ census

def stitch(a: Path, b: Path, dst: Path) -> tuple[int, int]:
    """Both patches in one grid, separated by GAP_ROWS invalid rows."""
    sa, sb = tifxyz.read(a), tifxyz.read(b)
    va, ua = sa.shape
    vb, ub = sb.shape
    V, U = va + GAP_ROWS + vb, max(ua, ub)
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    for i, ax in enumerate("xyz"):
        out = np.full((V, U), -1.0, np.float32)
        pa = np.asarray(sa.points[..., i], np.float32).copy()
        pa[~sa.valid] = -1.0
        out[:va, :ua] = pa
        pb = np.asarray(sb.points[..., i], np.float32).copy()
        pb[~sb.valid] = -1.0
        out[va + GAP_ROWS:, :ub] = pb
        tifffile.imwrite(dst / f"{ax}.tif", out)
    (dst / "meta.json").write_text(json.dumps(
        {"format": "tifxyz", "scale": [1, 1],
         "note": "two patches stitched for a pairwise census; not a merge"}))
    return sa.n_valid, sb.n_valid


def census_pair(a: Path, b: Path, work: Path, tag: str) -> dict:
    dst = work / f"stitch_{tag}.tifxyz"
    na, nb = stitch(a, b, dst)
    total, per = 0, {}
    for d in (0, 1):
        _, c = pipeline.run_engine(dst, f"pair{tag}", work, d, {"threads": 1})
        per[f"d{d}"] = int(c["transverse"])
        total += int(c["transverse"])
    shutil.rmtree(dst, ignore_errors=True)
    return {"valid_a": na, "valid_b": nb, "transverse_d0": per["d0"],
            "transverse_d1": per["d1"], "transverse_both": total,
            "cross": total > 0}


def stage_census(pairs_file: Path, out: Path, sample: int, seed: int,
                 jobs: int, data: Path, work: Path) -> None:
    # The pair file is ~6 GB at 27.8M rows. Reading it whole to draw a few
    # hundred rows costs tens of gigabytes of RAM and risks the run being
    # killed halfway. Reservoir-sample in one streaming pass instead: memory
    # is O(sample), the draw is still uniform, and it is reproducible for a
    # given seed because the reservoir depends only on the row order.
    rng = random.Random(seed)
    total = 0
    reservoir: list[dict] = []
    with pairs_file.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            total += 1
            if sample <= 0:
                reservoir.append(json.loads(line))
                continue
            if len(reservoir) < sample:
                reservoir.append(json.loads(line))
            else:
                j = rng.randrange(total)
                if j < sample:
                    reservoir[j] = json.loads(line)
    pairs, chosen = reservoir, reservoir
    print(f"streamed {total} pairs, reservoir {len(reservoir)}")
    done = set()
    if out.exists():
        done = {(r["a"], r["b"]) for r in
                (json.loads(l) for l in out.read_text().splitlines() if l.strip())}
    todo = [p for p in chosen if (p["a"], p["b"]) not in done]
    print(f"pairs to census: {len(todo)}")

    work.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    n_cross = n_ok = 0

    def one(p: dict, i: int) -> dict:
        rec = dict(p)
        for key in ("a", "b"):
            if not (data / p[key]).is_dir() and download(p[key], data) is None:
                rec["status"] = "download_failed"
                return rec
        try:
            rec.update(census_pair(data / p["a"], data / p["b"], work,
                                   f"{i:06d}"))
            rec["status"] = "ok"
        except Exception as e:                              # noqa: BLE001
            rec["status"] = f"census_failed: {type(e).__name__}: {e}"
        return rec

    with out.open("a") as sink, cf.ThreadPoolExecutor(jobs) as pool:
        futs = [pool.submit(one, p, i) for i, p in enumerate(todo)]
        for k, f in enumerate(cf.as_completed(futs), 1):
            rec = f.result()
            sink.write(json.dumps(rec) + "\n")
            sink.flush()
            if rec.get("status") == "ok":
                n_ok += 1
                n_cross += bool(rec["cross"])
            if k % 50 == 0:
                el = time.time() - t0
                print(f"  {k}/{len(todo)}  crossing={n_cross}/{n_ok}  "
                      f"{k/el:.1f}/s", flush=True)

    print(f"\npairs censused    {n_ok}")
    print(f"pairs that cross  {n_cross}"
          + (f"  ({100*n_cross/n_ok:.1f}%)" if n_ok else ""))
    print(f"wrote {out}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["bboxes", "pairs", "census"])
    ap.add_argument("--jobs", type=int, default=32)
    ap.add_argument("--sample", type=int, default=300)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-overlap", type=float, default=1e6)
    ap.add_argument("--pairs-file", default="out/patches/pairs.jsonl")
    ap.add_argument("--out", default="out/patches/pair_census.jsonl")
    ap.add_argument("--data", default="out/patches/pairdata")
    ap.add_argument("--work", default="out/patches/pairwork")
    a = ap.parse_args()

    BOXES.parent.mkdir(parents=True, exist_ok=True)
    if a.stage == "bboxes":
        stage_bboxes(INDEX.read_text().split(), a.jobs)
    elif a.stage == "pairs":
        stage_pairs(Path(a.pairs_file), a.min_overlap)
    else:
        stage_census(Path(a.pairs_file), Path(a.out), a.sample, a.seed,
                     a.jobs, Path(a.data), Path(a.work))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
