"""Rigid-infeasibility lower bounds for crossing events, corpus-wide.

For each event, the certificate-grade bound (round 21): actual crossing
triangle pairs (the census rows' 2x2 combos that geometrically intersect,
recomputed on the measured grid), per-pair minimum-norm Minkowski exit with
explicit facet normalization, L_safe = L - stated numerical allowance, and
the maximizing pair/facet as witness. Verdict per event against the
admissible relative budget (2 vx at lambda=0.5, strict per-vertex tier):

  certified_infeasible  (L_safe > admissible: no rigid relative core
                         translation within the tier can clear it)
  inconclusive          (bound below the budget -- says nothing)

Model scope, stated on every record: rules out RIGID RELATIVE CORE
TRANSLATION only; not non-rigid deformation, remeshing, or cutting.

    uv run python bench/min_exit_sweep.py                 # whole corpus
    uv run python bench/min_exit_sweep.py --segment NAME  # one certificate
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
from windcheck import tifxyz                                    # noqa: E402
from windcheck.check import load_pairs                          # noqa: E402
from windcheck.clearance import min_exit_from_pairs             # noqa: E402
from windcheck.intrinsic import (SurfaceGraph, _tri_tri_segment,  # noqa: E402
                                 oriented_events)

CORPORA = [
    ("Scroll 1", "data/scroll1_tifxyz", "20230205180739", "out/crossing_s1"),
    ("Scroll 5", "data/scroll5_tifxyz", "20241024131839", "out/crossing"),
    ("PHerc0139", "data/PHerc0139_tifxyz", "20250728140407", "out/crossing_0139"),
    ("PHerc0814", "data/PHerc0814_tifxyz", "20250804134230", "out/crossing_0814"),
    ("PHerc1667", "data/PHerc1667_tifxyz", "20231117161658", "out/crossing_1667"),
]
ADMISSIBLE_REL_VX = 2.0


def event_bounds(g: SurfaceGraph, X: np.ndarray, rec, ev) -> tuple:
    """Actual intersecting triangle pairs of the event's census rows."""
    pairs = []
    for i in ev["rows"]:
        r = rec[i]
        for t1 in g.quad_triangles(int(r["v1"]), int(r["u1"])):
            for t2 in g.quad_triangles(int(r["v2"]), int(r["u2"])):
                if (t1 >= 0).all() and (t2 >= 0).all() and \
                   _tri_tri_segment(X[t1], X[t2]) is not None:
                    pairs.append((X[t1], X[t2]))
    if not pairs:
        return 0.0, 0.0, None, 0
    coords = np.abs(np.concatenate([np.concatenate(p) for p in pairs[:8]]))
    ulp = float(np.max(np.spacing(coords.astype(np.float32))))
    allowance = max(8 * ulp, 1e-3)
    L, L_safe, wit = min_exit_from_pairs(pairs, allowance)
    return L, L_safe, {"allowance_vx": allowance, **(wit or {})}, len(pairs)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--segment", default="")
    ap.add_argument("--out", type=Path, default=Path("out/min_exit_sweep.jsonl"))
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if a.out.exists() and not a.segment:
        for x in a.out.read_text().splitlines():
            if x:
                done.add(json.loads(x)["segment"])
    commit = subprocess.run(["git", "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    t0 = time.time()
    for corpus, root, volume, work in CORPORA:
        rootp = Path(root)
        if not rootp.exists():
            continue
        for d in sorted(p for p in rootp.iterdir() if p.is_dir()):
            if a.segment and a.segment not in d.name:
                continue
            if d.name in done:
                continue
            m = sorted(d.glob(f"mesh/*{volume}*.tifxyz"))
            if not m or not (m[0] / "x.tif").exists():
                continue
            s = tifxyz.read(m[0])
            nv, nu = s.shape
            rows_out = []
            for diag in (0, 1):
                csv = Path(work) / f"{d.name[:40]}_d{diag}.csv"
                rec = load_pairs(csv)
                if len(rec):
                    rec = rec[(rec["v1"] < nv-1) & (rec["v2"] < nv-1)
                              & (rec["u1"] < nu-1) & (rec["u2"] < nu-1)]
                if not len(rec):
                    continue
                g = SurfaceGraph(s.points, s.valid, diag)
                gv, gu = np.nonzero(g.idx >= 0)
                X = np.empty((g.n, 3))
                X[g.idx[gv, gu]] = np.asarray(s.points, np.float64)[gv, gu]
                for k, ev in enumerate(oriented_events(rec)):
                    if ev["ambiguous"]:
                        continue
                    L, L_safe, wit, npairs = event_bounds(g, X, rec, ev)
                    rows_out.append({
                        "diagonal": diag, "event_index": k,
                        "n_rows": len(ev["rows"]), "n_tri_pairs": npairs,
                        "L_vx": round(L, 6), "L_safe_vx": round(L_safe, 6),
                        "witness": wit,
                        "verdict": ("certified_infeasible"
                                    if L_safe > ADMISSIBLE_REL_VX
                                    else "inconclusive"),
                    })
            row = {"corpus": corpus, "segment": d.name,
                   "admissible_rel_vx": ADMISSIBLE_REL_VX,
                   "model_scope": ("rules out rigid relative core "
                                   "translation <= admissible only; not "
                                   "non-rigid deformation, remeshing, or "
                                   "cutting"),
                   "code_commit": commit,
                   "mesh_x_sha256": hashlib.sha256(
                       (m[0] / "x.tif").read_bytes()).hexdigest(),
                   "events": rows_out}
            with a.out.open("a") as f:
                f.write(json.dumps(row) + "\n")
            ninf = sum(1 for e in rows_out
                       if e["verdict"] == "certified_infeasible")
            print(f"{corpus:10s} {d.name[:36]:38s} events {len(rows_out):5d} "
                  f"certified-infeasible {ninf:4d} [{time.time()-t0:6.0f}s]",
                  flush=True)


if __name__ == "__main__":
    main()
