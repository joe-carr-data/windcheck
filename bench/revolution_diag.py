"""Diagnose WHY the absolute-mm separation criterion inverts across scrolls.

Step 0 of the handoff. The 17.2 mm / 135.7 mm gap is a Scroll 5 observation:
PHerc0139's labelled windings reach 89 mm and PHerc1667's reach 1237 mm. Before
choosing a normalisation, look at the shape of the offending crossings rather
than at their headline number -- the lesson from the first census, which spent an
hour measuring local mesh roughness because the count looked plausible.

For each segment this prints:

  shape           grid rows x columns
  u_rev           columns per revolution, measured from the geometry itself
                  (see revolution_period below)
  span/rev        how many revolutions the segment covers
  du_max          largest column separation among crossing pairs
  du_max/u_rev    that separation as a fraction of a revolution
  end             fraction of far pairs whose two ends sit at OPPOSITE ends of
                  the segment's own column range -- the signature of a winding
                  that closes on itself rather than of a doubled trace
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crossing_analyse import load                              # noqa: E402

from windcheck import tifxyz                                   # noqa: E402


def turning_profile(P: np.ndarray, V: np.ndarray) -> np.ndarray | None:
    """Unwrapped angle, in radians, at each grid column.

    Dividing a column offset by a period assumes columns are spaced uniformly in
    angle. They are not: these grids come out of a flattening step, so the
    parameterisation is uniform in the flattened space rather than around the
    scroll. Reading the angle off directly removes the assumption instead of
    hoping it holds.

    Returns one value per column of `P`, with columns that carry no valid
    surface filled by interpolation from their neighbours. `None` when the
    centreline is too short or turns too little to carry an angle.
    """
    cols, have = [], []
    for u in range(P.shape[1]):
        m = V[:, u]
        if m.sum() >= 3:
            cols.append(np.median(P[m, u], axis=0))
            have.append(u)
    if len(have) < 8:
        return None
    c = np.asarray(cols) - np.mean(cols, axis=0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    x, y = c @ vt[0], c @ vt[1]
    ang = np.unwrap(np.arctan2(y, x))
    if abs(ang[-1] - ang[0]) < 0.35:          # under ~20 degrees: no angle
        return None
    return np.interp(np.arange(P.shape[1]), np.asarray(have), ang)


def revolution_period(P: np.ndarray, V: np.ndarray) -> float:
    """Columns per revolution: the PCA-centroid polar-angle period.

    Not a tangent-turning integral, which an earlier docstring claimed. What it
    actually does, stated so the assumption is visible: take one 3D point per
    grid column, project onto the two leading principal directions of that
    centreline, unwrap the polar angle about the centreline's own centroid, and
    divide the column range by the total angle swept. No axis is fitted and no
    circle assumed -- both failed their positive control on this scroll -- but
    there IS an implicit centre, the PCA centroid, and it is well-conditioned
    only when the centreline wraps around it. That is why short arcs are
    rejected and why the angular measure derived from the same profile is
    unusable at short range.

    The column *range* is used, not the count of retained columns: internal
    missing columns would otherwise shorten the numerator and shrink the period.

    Returns nan when the surface turns too little to carry a period, which is
    itself informative -- such a segment cannot close on itself.
    """
    theta = turning_profile(P, V)
    if theta is None:
        return float("nan")
    ucols = np.nonzero(V.any(0))[0]
    if len(ucols) < 8:
        return float("nan")
    u_lo, u_hi = int(ucols[0]), int(ucols[-1])
    turn = float(np.abs(theta[u_hi] - theta[u_lo]))
    if turn < 0.35:                      # under ~20 degrees: not measurable
        return float("nan")
    return (u_hi - u_lo) * 2.0 * np.pi / turn


def diag(seg_dir: Path, csv: Path, volume: str,
         census_grid: dict[str, list] | None = None) -> dict | None:
    m = sorted(seg_dir.glob(f"mesh/*{volume}*.tifxyz"))
    if not m:
        return None
    s = tifxyz.read(m[0])
    P, V = s.points, s.valid

    # The CSV holds (v, u) indices into the grid the CENSUS measured. Reading a
    # different resolution of the same segment here silently reinterprets those
    # indices against the wrong grid, and the result looks entirely plausible --
    # it produced an 11-segment "half a revolution" cluster on PHerc1667 that
    # survived a second normalisation before the shape mismatch was spotted.
    # Nothing downstream can detect it, so it is caught here or not at all.
    if census_grid is not None:
        want = census_grid.get(seg_dir.name)
        if want is None:
            return None
        if list(P.shape[:2]) != list(want["grid"]):
            raise SystemExit(
                f"{seg_dir.name}: census measured a {want['grid']} grid but "
                f"{m[0].name} is {list(P.shape[:2])}. The crossing indices "
                f"belong to a different resolution -- pass the --volume that "
                f"the census used ({want.get('mesh', {}).get('volume')}).")
        cm = want.get("mesh") or {}
        if cm.get("name") and cm["name"] != m[0].name:
            raise SystemExit(
                f"{seg_dir.name}: census ran on {cm['name']} but this is "
                f"{m[0].name}. Same grid shape is not the same mesh.")
    u_rev = revolution_period(P, V)
    ucols = np.nonzero(V.any(0))[0]
    if len(ucols) == 0:
        return None
    u_lo, u_hi = int(ucols[0]), int(ucols[-1])
    span = u_hi - u_lo + 1

    rec = load(csv)
    out = {"segment": seg_dir.name, "shape": list(P.shape[:2]),
           "u_span": span, "u_rev": round(u_rev, 1) if u_rev == u_rev else None,
           "span_per_rev": round(span / u_rev, 3) if u_rev == u_rev else None,
           "pairs": int(len(rec))}

    # The angular span is a property of the SURFACE, so it is measured before
    # looking at crossings. A segment with no crossings at all is the most
    # informative case there is -- it is the null control -- and computing the
    # span only when there is something to divide would drop exactly those.
    theta = turning_profile(P, V)
    span_rev = (float(np.abs(theta[u_hi] - theta[u_lo])) / (2.0 * np.pi)
                if theta is not None else None)

    if len(rec) == 0:
        # Same keys as the populated branch. An absent key and a zero are the
        # same fact here -- no crossings -- and a consumer that has to tell them
        # apart with .get() will eventually get it wrong.
        out.update(du_max=0, dv_max=0, sep_max=0, du_max_per_rev=0.0,
                   sep_max_per_span=0.0, end_frac=0.0,
                   ang_sep_max_rev=0.0, ang_sep_p99_rev=0.0,
                   ang_span_rev=(round(span_rev, 3) if span_rev is not None
                                 else None))
        return out

    # Every index must address the grid actually loaded. The shape check above
    # catches a wholesale resolution swap; this catches the subtler case of two
    # meshes with the same shape, and costs nothing.
    nv, nu = P.shape[:2]
    for f, n in (("v1", nv), ("v2", nv), ("u1", nu), ("u2", nu)):
        if rec[f].size and (rec[f].min() < 0 or rec[f].max() >= n):
            raise SystemExit(
                f"{seg_dir.name}: CSV column {f} ranges "
                f"{rec[f].min()}..{rec[f].max()}, outside this {nv}x{nu} grid.")

    du = np.abs(rec["u1"].astype(np.int64) - rec["u2"].astype(np.int64))
    dv = np.abs(rec["v1"].astype(np.int64) - rec["v2"].astype(np.int64))
    sep = np.maximum(du, dv)

    # Angular separation read straight off the turning profile: no assumption
    # that a column offset is proportional to an angle.
    out["ang_sep_max_rev"] = out["ang_sep_p99_rev"] = out["ang_span_rev"] = None
    if theta is not None:
        n = len(theta)
        u1 = np.clip(rec["u1"].astype(np.int64), 0, n - 1)
        u2 = np.clip(rec["u2"].astype(np.int64), 0, n - 1)
        dth = np.abs(theta[u1] - theta[u2]) / (2.0 * np.pi)
        # Self-consistency: two parts of one segment cannot be further apart in
        # angle than the segment itself spans. When they are, the unwrapped
        # angle has jumped -- the centreline wandered near its own projected
        # centre, where arctan2 is unstable -- and the whole profile is
        # unusable. This gate independently rejects exactly the two segments
        # the two period estimators already disagreed 16x about, which is why
        # it is trusted: two unrelated checks naming the same cases.
        if dth.max() <= span_rev * 1.02 + 1e-6:
            out["ang_sep_max_rev"] = round(float(dth.max()), 3)
            out["ang_sep_p99_rev"] = round(float(np.percentile(dth, 99)), 3)
            out["ang_span_rev"] = round(span_rev, 3)
        else:
            out["angle_rejected"] = (
                f"separation {dth.max():.3f} exceeds span {span_rev:.3f} rev")
    # "at opposite ends": both quads within 15% of the span of their own end.
    edge = max(1, int(0.15 * span))
    lo = np.minimum(rec["u1"], rec["u2"])
    hi = np.maximum(rec["u1"], rec["u2"])
    far = sep >= np.percentile(sep, 90)
    endish = (lo <= u_lo + edge) & (hi >= u_hi - edge)
    out.update(
        du_max=int(du.max()), dv_max=int(dv.max()), sep_max=int(sep.max()),
        du_max_per_rev=(round(float(du.max()) / u_rev, 3)
                        if u_rev == u_rev else None),
        sep_max_per_span=round(float(sep.max()) / span, 3),
        end_frac=round(float(endish[far].mean()), 3) if far.any() else 0.0,
    )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--volume", required=True)
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--json", type=Path)
    ap.add_argument("--census", type=Path, default=None,
                    help="census JSON whose CSVs are in --dir; its recorded "
                         "grid shape is checked against the mesh read here")
    a = ap.parse_args()

    census_grid = None
    if a.census is None:
        for cand in ("census_v3.json", "census.json"):
            if (a.dir / cand).exists():
                a.census = a.dir / cand
                break
    if a.census and a.census.exists():
        census_grid = {r["segment"]: r
                       for r in json.loads(a.census.read_text())}
        print(f"# grid shapes checked against {a.census}")
    else:
        print("# WARNING: no census JSON found; grid shapes are UNCHECKED")

    rows = []
    for d in sorted(a.root.iterdir()):
        if not d.is_dir():
            continue
        r = diag(d, a.dir / f"{d.name[:40]}_d0.csv", a.volume, census_grid)
        if r:
            rows.append(r)
            print(f"{r['segment'][:38]:40s} shape={str(r['shape']):>14s} "
                  f"span={r['u_span']:5d} rev={str(r['u_rev']):>8s} "
                  f"s/r={str(r['span_per_rev']):>6s} "
                  f"sep={r.get('sep_max',0):5d} "
                  f"du/rev={str(r.get('du_max_per_rev')):>6s} "
                  f"dth/2pi={str(r.get('ang_sep_max_rev')):>6s} "
                  f"end={r.get('end_frac',0):.2f}", flush=True)
    if a.json:
        a.json.write_text(json.dumps(rows, indent=2))
        print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
