"""Evaluate the pre-registered claim: does disagreement predict error?

Gate, fixed in notes/PREREG-DISAGREEMENT.md before any statistic met any
label: AUC > 0.75 stands, 0.65 < AUC <= 0.75 is suggestive only, AUC
<= 0.65 is dead with no reframing.

Labels come from the pinned corpus census and enter only here. The
disagreement field was computed from patches alone.

Two confound checks are run whether or not the headline passes:

* geometry-only baseline. If patch area, overlap volume and covering-pair
  count reach within 0.05 AUC of the disagreement statistic, the signal
  is explained by how the patches are shaped rather than by what they
  say, and the result is void.
* coverage strata. If high-disagreement cubes are simply cubes with more
  pairs, the effect is an artifact of coverage, so AUC is recomputed
  within strata of covering-pair count.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import tifffile  # noqa: E402

CROSSING_DIRS = ["out/crossing", "out/crossing_s1"]


def auc(score: np.ndarray, label: np.ndarray) -> float:
    """Rank-based AUC, ties averaged."""
    pos, neg = score[label == 1], score[label == 0]
    if not len(pos) or not len(neg):
        return float("nan")
    order = np.argsort(np.concatenate([pos, neg]), kind="mergesort")
    ranks = np.empty(len(order), float)
    ranks[order] = np.arange(1, len(order) + 1)
    s = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    # average ties
    allv = np.concatenate([pos, neg])
    for v in np.unique(allv):
        m = allv == v
        if m.sum() > 1:
            ranks[m] = ranks[m].mean()
    return float((ranks[s == 1].sum() - len(pos) * (len(pos) + 1) / 2)
                 / (len(pos) * len(neg)))


def crossing_points(index: dict, cube: int) -> tuple[set, dict]:
    """Cubes containing a published self-intersection, and trace coverage."""
    positives, covered = set(), {}
    for seg, v in index["segments"].items():
        if v["scroll"] != "Scroll 1":
            continue
        mesh = Path(v["original_mesh"])
        if not mesh.is_dir():
            continue
        x = tifffile.imread(mesh / "x.tif").astype(np.float32)
        y = tifffile.imread(mesh / "y.tif").astype(np.float32)
        z = tifffile.imread(mesh / "z.tif").astype(np.float32)
        valid = ~((x == -1) & (y == -1) & (z == -1))
        P = np.stack([x, y, z], -1)
        # every cube this trace passes through is evaluable
        pts = P[valid]
        for c in np.unique((pts // cube).astype(np.int32), axis=0):
            covered[tuple(c)] = covered.get(tuple(c), 0) + 1
        stem = seg
        for d in CROSSING_DIRS:
            for f in Path(d).glob(f"{stem}_d*.csv"):
                with f.open() as fh:
                    for r in csv.DictReader(fh):
                        if r["verdict"] != "transverse":
                            continue
                        vv, uu = int(r["v1"]), int(r["u1"])
                        if vv < P.shape[0] and uu < P.shape[1] and valid[vv, uu]:
                            positives.add(tuple((P[vv, uu] // cube).astype(int)))
    return positives, covered


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--field", default="out/disagree/field.jsonl")
    ap.add_argument("--cube", type=int, default=256)
    ap.add_argument("--min-pairs", type=int, default=3)
    ap.add_argument("--out", default="out/disagree/evaluation.json")
    a = ap.parse_args()

    recs = [json.loads(l) for l in Path(a.field).read_text().splitlines()
            if l.strip()]
    ok = [r for r in recs if r.get("status") == "ok"]
    print(f"field records {len(recs)}, usable {len(ok)}")

    by_cube: dict[tuple, list] = {}
    for r in ok:
        c = tuple((np.array(r["centroid"]) // a.cube).astype(int))
        by_cube.setdefault(c, []).append(r)

    index = json.loads(Path("out/release/index.json").read_text())
    positives, covered = crossing_points(index, a.cube)
    print(f"cubes with disagreement data {len(by_cube)}; "
          f"cubes touched by a published trace {len(covered)}; "
          f"cubes containing a published crossing {len(positives)}")

    rows = []
    for c, rs in by_cube.items():
        if len(rs) < a.min_pairs or c not in covered:
            continue
        rows.append({
            "cube": list(map(int, c)),
            "n_pairs": len(rs),
            "disagreement": float(np.median([r["median_vx"] for r in rs])),
            "mean_overlap": float(np.mean([r["overlap_vx3"] for r in rs])),
            "mean_points": float(np.mean([r["n_points"] for r in rs])),
            "label": int(c in positives)})

    if len(rows) < 20:
        print(f"\nONLY {len(rows)} evaluable cubes: too few to evaluate. "
              "Increase --sample in the field stage or widen --cube.")
        Path(a.out).write_text(json.dumps(
            {"status": "insufficient_data", "evaluable_cubes": len(rows)},
            indent=1) + "\n")
        return 2

    score = np.array([r["disagreement"] for r in rows])
    label = np.array([r["label"] for r in rows])
    npos, nneg = int(label.sum()), int((1 - label).sum())
    A = auc(score, label)

    rng = np.random.default_rng(0)
    null = np.array([auc(score, rng.permutation(label)) for _ in range(1000)])
    p = float(np.mean(null >= A))

    geo = {}
    for k in ("n_pairs", "mean_overlap", "mean_points"):
        geo[k] = auc(np.array([r[k] for r in rows]), label)
    best_geo = max(abs(v - 0.5) for v in geo.values()) + 0.5

    strata = {}
    q = np.quantile([r["n_pairs"] for r in rows], [0.5])
    for name, m in (("low_coverage", np.array([r["n_pairs"] for r in rows]) <= q[0]),
                    ("high_coverage", np.array([r["n_pairs"] for r in rows]) > q[0])):
        if m.sum() > 10 and 0 < label[m].sum() < m.sum():
            strata[name] = auc(score[m], label[m])

    verdict = ("STANDS" if A > 0.75 else
               "SUGGESTIVE" if A > 0.65 else "DEAD")
    confounded = best_geo >= A - 0.05

    print(f"\nevaluable cubes {len(rows)}  (positive {npos}, negative {nneg})")
    print(f"AUC                     {A:.3f}")
    print(f"permutation null        mean {null.mean():.3f}, p = {p:.4f}")
    print(f"geometry-only baselines { {k: round(v, 3) for k, v in geo.items()} }")
    print(f"within-coverage strata  { {k: round(v, 3) for k, v in strata.items()} }")
    print(f"\nGATE: AUC > 0.75 stands, <= 0.65 dead   ->  {verdict}")
    if confounded:
        print("CONFOUNDED: geometry alone reaches within 0.05 of the "
              "disagreement statistic. Per the pre-registration the result "
              "is void regardless of AUC.")

    Path(a.out).write_text(json.dumps({
        "cube_vx": a.cube, "min_pairs": a.min_pairs,
        "evaluable_cubes": len(rows), "positives": npos, "negatives": nneg,
        "auc": A, "permutation_p": p, "permutation_null_mean": float(null.mean()),
        "geometry_only_auc": geo, "coverage_strata_auc": strata,
        "verdict": verdict, "confounded": bool(confounded),
        "gate": "AUC>0.75 stands, 0.65-0.75 suggestive, <=0.65 dead",
        "rows": rows}, indent=1) + "\n")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
