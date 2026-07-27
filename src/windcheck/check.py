"""One command, one segment, one answer.

Everything else in this repository is a bench script that reproduces a figure in
the write-up. This is the part someone else runs:

    windcheck check data/scroll5_tifxyz/20251115002745-auto_grown_..._flatboi

    20251115002745-auto_grown_20251115002740308_5_flatboi
      grid                637 x 3065          triangles  3,535,554
      covering span       5.91 revolutions
      widest separation   4.91 revolutions
      crossing events     563   (181 separated by 1.6 revolutions or more)
      VERDICT             self-intersection present; widest separation 4.91
                          revolutions along the trace's own parameter

      certificate  out/check/<name>_certificate.json
      overlay      out/check/<name>_points.json    <- opens in VC3D

The verdict is about the surface, not about the tracer. This measures that a
surface meets itself and how far apart along its own parameter; it does not
establish that any crossing is a tracing error, and the certificate carries that
caveat on its face so the number cannot travel without it.

Self-contained on purpose: it imports nothing from `bench/`, because those
scripts exist to reproduce figures and should be free to change.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from . import atlas, classify, tifxyz
from .certificate import (WRAP_SCALE_CUT_REV, certificate, write_collection)

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / "engines" / "selfcross"

PAIR_DTYPE = [("v1", "i4"), ("u1", "i4"), ("v2", "i4"), ("u2", "i4"),
              ("pen", "f8"), ("ang", "f8")]


class _Entry:
    def __init__(self, path: Path):
        self.path, self.winding = Path(path), None


def find_mesh(target: Path, volume: str = "") -> Path | None:
    """Accept a .tifxyz directly, or a segment directory to search inside."""
    target = Path(target)
    if (target / "meta.json").exists() and (target / "x.tif").exists():
        return target
    for pat in (f"mesh/*{volume}*.tifxyz", f"*{volume}*.tifxyz"):
        hits = sorted(target.glob(pat))
        if hits:
            return hits[0]
    return None


def load_pairs(csv: Path) -> np.ndarray:
    if not csv.exists():
        return np.empty(0, dtype=PAIR_DTYPE)
    rows = []
    with csv.open() as fh:
        has_margin = "penetration" in fh.readline()
        for line in fh:
            p = line.rstrip("\n").split(",")
            if p[4] != "transverse":
                continue
            rows.append((int(p[0]), int(p[1]), int(p[2]), int(p[3]),
                         float(p[5]) if has_margin else np.nan,
                         float(p[6]) if has_margin else np.nan))
    return np.array(rows, dtype=PAIR_DTYPE)


def revolution_period(P: np.ndarray, V: np.ndarray) -> float:
    """Columns per revolution, from the turning of the surface's own centreline.

    No axis is fitted; axis fitting failed its own positive control on this data.
    Validated externally against published winding counts on 31 of 33 Scroll 1
    segments at r = 0.9999 (`bench/winding_control.py`).
    """
    cols, have = [], []
    for u in range(P.shape[1]):
        m = V[:, u]
        if m.sum() >= 3:
            cols.append(np.median(P[m, u], axis=0))
            have.append(u)
    if len(have) < 8:
        return float("nan")
    c = np.asarray(cols) - np.mean(cols, axis=0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    ang = np.unwrap(np.arctan2(c @ vt[1], c @ vt[0]))
    turn = float(abs(ang[-1] - ang[0]))
    if turn < 0.35:                       # under ~20 degrees: no period
        return float("nan")
    return (have[-1] - have[0]) * 2.0 * np.pi / turn


def neighbour_period(mesh: Path, out: Path, stride: int = 4,
                     threads: int = 0) -> float:
    """Second, centre-free period estimate: columns to the adjacent wrap.

    For each grid point, the physically nearest part of the SAME surface at
    least a few columns away is the neighbouring wrap, so the column offset to
    it is one revolution. Returns nan when there is no neighbouring wrap inside
    the segment -- a single winding has none, which is a property of the segment
    and not a failure.
    """
    eng = ROOT / "engines" / "atlas_query"
    if not eng.exists():
        return float("nan")
    s = tifxyz.read(mesh)
    v, u = np.nonzero(s.valid[::stride, ::stride])
    if len(v) < 2000:
        return float("nan")
    V, U = v * stride, u * stride
    gap = max(s.shape[1] // 100, 3)
    out.mkdir(parents=True, exist_ok=True)
    atlas.write_atlas([_Entry(mesh)], out / "_np_atlas.bin")
    atlas.write_queries_grouped(s.points[V, U], U, out / "_np_query.bin")
    r = atlas.run_engine(out / "_np_atlas.bin", out / "_np_query.bin",
                         out / "_np_result.bin", threads=threads, exclude_u=gap)
    d, w1 = r["d1"], r["w1"]
    m = np.isfinite(d) & (d < 24.0)
    if m.sum() < 200:
        return float("nan")
    return float(np.median(np.abs(w1[m] - U[m])))


def components(rec: np.ndarray) -> dict:
    """Grid-connected regions of participating quads, by union-find."""
    quads = set()
    for r in rec:
        quads.add((int(r["v1"]), int(r["u1"])))
        quads.add((int(r["v2"]), int(r["u2"])))
    parent: dict = {q: q for q in quads}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (v, u) in quads:
        for dv in (-1, 0, 1):
            for du in (-1, 0, 1):
                n = (v + dv, u + du)
                if n in parent:
                    a, b = find((v, u)), find(n)
                    if a != b:
                        parent[a] = b
    return {q: find(q) for q in quads}


def events(rec: np.ndarray, lab: dict) -> set:
    """A crossing EVENT is a pair of regions, not a pair of triangles."""
    return {tuple(sorted((lab[(int(r["v1"]), int(r["u1"]))],
                          lab[(int(r["v2"]), int(r["u2"]))])))
            for r in rec}


def quad_area_mm2(P, V, vx_um) -> np.ndarray:
    p00, p10, p01, p11 = P[:-1, :-1], P[1:, :-1], P[:-1, 1:], P[1:, 1:]
    ok = V[:-1, :-1] & V[1:, :-1] & V[:-1, 1:] & V[1:, 1:]
    a = 0.5 * np.linalg.norm(np.cross(p01 - p00, p11 - p00), axis=-1) \
        + 0.5 * np.linalg.norm(np.cross(p11 - p00, p10 - p00), axis=-1)
    a[~ok] = 0.0
    return a * (vx_um / 1000.0) ** 2


def analyse(target: Path, out: Path, volume: str = "", threads: int = 0,
            cell: float = 40.0, maxedge: float = 60.0) -> dict | None:
    mesh = find_mesh(target, volume)
    if mesh is None:
        print(f"windcheck: no .tifxyz found under {target}", file=sys.stderr)
        return None
    if not ENGINE.exists():
        print(f"windcheck: engine not built: {ENGINE}\n"
              "  clang++ -O3 -std=c++17 -pthread -o engines/selfcross "
              "engines/selfcross.cpp", file=sys.stderr)
        return None

    name = mesh.parent.parent.name or mesh.stem
    out.mkdir(parents=True, exist_ok=True)
    s = tifxyz.read(mesh)
    P, V = s.points, s.valid

    m = __import__("re").search(r"-(\d+\.?\d*)um\.tifxyz$", mesh.name)
    vx_um = float(m.group(1)) if m else float("nan")

    atlas.write_atlas([_Entry(mesh)], out / "_atlas.bin")
    csv = out / f"{name[:40]}_pairs.csv"
    proc = subprocess.run(
        [str(ENGINE), str(out / "_atlas.bin"), str(csv), str(threads),
         str(cell), "1", "0", str(maxedge)],
        capture_output=True, text=True, check=True)
    stats = json.loads(proc.stdout.strip().splitlines()[-1])

    rec = load_pairs(csv)
    period = revolution_period(P, V)
    ucols = np.nonzero(V.any(0))[0]
    span_rev = (float(ucols[-1] - ucols[0]) / period) if period == period else None

    sep_rev = n_events = n_far = 0
    pts: list = []
    if len(rec) and period == period:
        du = np.abs(rec["u1"].astype(np.int64) - rec["u2"].astype(np.int64))
        sep_rev = float(du.max()) / period
        n_events = len(events(rec, components(rec)))
        far = rec[du > WRAP_SCALE_CUT_REV * period]
        if len(far):
            l2 = components(far)
            n_far = len(events(far, l2))
            best: dict = {}
            for r in far:
                k = tuple(sorted((l2[(int(r["v1"]), int(r["u1"]))],
                                  l2[(int(r["v2"]), int(r["u2"]))])))
                if k not in best or r["pen"] > best[k]["pen"]:
                    best[k] = r
            pts = [P[int(r["v1"]), int(r["u1"])].tolist() for r in best.values()]
            assert len(pts) == n_far

    neigh = neighbour_period(mesh, out, threads=threads)
    pstatus = classify.period_status(period if period == period else None, neigh)

    area = quad_area_mm2(P, V, vx_um) if vx_um == vx_um else np.zeros((1, 1))
    part = 0.0
    if len(rec) and vx_um == vx_um:
        uq = np.unique(np.stack([np.concatenate([rec["v1"], rec["v2"]]),
                                 np.concatenate([rec["u1"], rec["u2"]])], 1),
                       axis=0)
        inb = (uq[:, 0] < area.shape[0]) & (uq[:, 1] < area.shape[1])
        part = float(area[uq[inb, 0], uq[inb, 1]].sum())

    if pts:
        write_collection(out / f"{name[:40]}_points.json",
                         f"windcheck-overlaps-{name[:22]}", pts,
                         metadata={"tool": "windcheck",
                                   "cut_revolutions": WRAP_SCALE_CUT_REV})

    cert = certificate(
        segment=name, mesh_path=mesh, voxel_um=vx_um,
        params={"exclude": 1, "cell": cell, "maxedge": maxedge,
                "touch_tol": 1e-3, "diagonal": 0},
        total_area_mm2=float(area.sum()),
        participating_quad_area_mm2=part,
        parameter_span_mm_estimate=0.0,
        separation_revolutions=(sep_rev or None),
        covering_span_revolutions=span_rev,
        revolution_period_columns={
            "turning_estimate": round(period, 1) if period == period else None,
            "neighbour_estimate": round(neigh, 1) if neigh == neigh else None,
            "method": "centreline turning, cross-checked against the column "
                      "offset to the neighbouring wrap"},
        period_status=pstatus,
        n_pairs=int(len(rec)), n_events=n_events, events_beyond_cut=n_far,
        median_penetration_vx=(float(np.nanmedian(rec["pen"])) if len(rec) else 0.0))
    (out / f"{name[:40]}_certificate.json").write_text(json.dumps(cert, indent=2))

    return {"name": name, "mesh": mesh, "grid": list(P.shape[:2]),
            "triangles": stats["triangles"], "pairs": int(len(rec)),
            "period": period, "neighbour": neigh, "period_status": pstatus,
            "span_rev": span_rev, "sep_rev": sep_rev or None,
            "band": classify.separation_band(sep_rev or None),
            "events": n_events, "events_beyond_cut": n_far,
            "verdict": cert["verdict"], "note": cert["note"],
            "out": out, "points": len(pts)}


def report(r: dict) -> None:
    print(f"\n{r['name']}")
    print(f"  grid                {r['grid'][0]} x {r['grid'][1]}"
          f"{'':6s}triangles  {r['triangles']:,}")
    print(f"  covering span       "
          f"{('%.2f revolutions' % r['span_rev']) if r['span_rev'] else 'not measurable'}")
    if r["pairs"] == 0:
        print("  self-overlap        none found")
    else:
        sep = (f"{r['sep_rev']:.2f} revolutions" if r["sep_rev"] else "n/a")
        if r["period_status"] != "agreed":
            sep += f"    (period {r['period_status']}; treat with caution)"
        print(f"  widest separation   {sep}")
        print(f"  crossing events     {r['events']:,}"
              f"   ({r['events_beyond_cut']:,} separated by "
              f"{classify.SEP_WIDE} revolutions or more)")
    print(f"  VERDICT             {r['verdict']}")
    print(f"\n  certificate  {r['out']}/{r['name'][:40]}_certificate.json")
    if r["points"]:
        print(f"  overlay      {r['out']}/{r['name'][:40]}_points.json"
              f"   <- open in VC3D ({r['points']} points)")
    print(f"\n  {r['note']}")
    print("  This measures that the surface meets itself, and how far apart "
          "along its own\n  parameter. It does not establish that any crossing "
          "is a tracing error.")
