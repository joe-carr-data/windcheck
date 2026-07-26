"""Per-trace self-overlap report, in millimetres rather than grid cells.

Three lessons forced this rewrite of the output:

  * pair counts are not defect counts -- one overlap makes thousands of pairs
  * event counts are sampling-dependent -- they scaled 15x with triangle count
    across a resolution change of the same surface
  * grid cells are not comparable between scrolls or volumes, because the voxel
    size differs (7.91, 9.362, 2.399, 1.129 um)

Physical extent was invariant to 0.3% across a 70x change in triangle count, so
it is the unit that actually measures the surface rather than our sampling of it.

What a mesher can act on:

  overlap_extent_mm     how far apart along the surface the two overlapping
                        parts lie -- for a doubled-back trace this is the span
                        of the duplication
  duplicated_area_mm2   physical area of surface taking part in overlaps
  duplicated_fraction   that area over the trace's total area, so traces of
                        different sizes are comparable
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crossing_analyse import components, events, load       # noqa: E402

from windcheck import tifxyz                                 # noqa: E402

RES_RE = re.compile(r"-(\d+\.?\d*)um\.tifxyz$")


def voxel_um(path: Path) -> float | None:
    m = RES_RE.search(path.name)
    return float(m.group(1)) if m else None


def quad_areas(P: np.ndarray, V: np.ndarray) -> np.ndarray:
    """Area of each quad, in voxels squared, via its two triangles."""
    p00, p10 = P[:-1, :-1], P[1:, :-1]
    p01, p11 = P[:-1, 1:], P[1:, 1:]
    ok = V[:-1, :-1] & V[1:, :-1] & V[:-1, 1:] & V[1:, 1:]
    a = 0.5 * np.linalg.norm(np.cross(p01 - p00, p11 - p00), axis=-1) \
      + 0.5 * np.linalg.norm(np.cross(p11 - p00, p10 - p00), axis=-1)
    a[~ok] = 0.0
    return a


def report(seg_dir: Path, csv: Path, volume: str, cut_mm: float) -> dict | None:
    m = sorted(seg_dir.glob(f"mesh/*{volume}*.tifxyz"))
    if not m:
        return None
    path = m[0]
    vx = voxel_um(path)
    if vx is None:
        return None
    s = tifxyz.read(path)
    P, V = s.points, s.valid
    area = quad_areas(P, V)
    total_mm2 = float(area.sum() * (vx / 1000.0) ** 2)

    # mean physical spacing of one grid step, for converting cell separation
    du = np.linalg.norm(P[:, 1:] - P[:, :-1], axis=-1)
    oku = V[:, 1:] & V[:, :-1]
    pitch_um = float(np.median(du[oku])) * vx if oku.any() else np.nan

    rec = load(csv)
    out = {"segment": seg_dir.name, "voxel_um": vx,
           "pitch_um": pitch_um, "total_area_mm2": round(total_mm2, 1),
           "pairs": int(len(rec))}
    if len(rec) == 0:
        out.update(overlap_extent_mm=0.0, duplicated_area_mm2=0.0,
                   duplicated_fraction=0.0, events_beyond_cut=0)
        return out

    sep_mm = np.maximum(np.abs(rec["v1"] - rec["v2"]),
                        np.abs(rec["u1"] - rec["u2"])) * pitch_um / 1000.0

    uniq = np.unique(np.stack([np.concatenate([rec["v1"], rec["v2"]]),
                               np.concatenate([rec["u1"], rec["u2"]])], 1), axis=0)
    inb = (uniq[:, 0] < area.shape[0]) & (uniq[:, 1] < area.shape[1])
    dup_mm2 = float(area[uniq[inb, 0], uniq[inb, 1]].sum() * (vx / 1000.0) ** 2)

    far = rec[sep_mm > cut_mm]
    n_far = 0
    if len(far):
        lab, _ = components(far)
        n_far = len(events(far, lab))

    out.update(
        overlap_extent_mm=round(float(sep_mm.max()), 1),
        overlap_extent_p99_mm=round(float(np.percentile(sep_mm, 99)), 1),
        duplicated_area_mm2=round(dup_mm2, 2),
        duplicated_fraction=round(dup_mm2 / total_mm2, 5) if total_mm2 else 0.0,
        events_beyond_cut=n_far,
        median_penetration_vx=round(float(np.nanmedian(rec["pen"])), 3),
    )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/scroll5_tifxyz"))
    ap.add_argument("--volume", default="20241024131839")
    ap.add_argument("--dir", type=Path, default=Path("out/crossing"))
    ap.add_argument("--cut-mm", type=float, default=30.0,
                    help="physical separation above which an overlap is "
                         "wrap-scale rather than local")
    ap.add_argument("--json", type=Path,
                    default=Path("out/crossing/physical.json"))
    a = ap.parse_args()

    rows = []
    for d in sorted(a.root.iterdir()):
        if not d.is_dir():
            continue
        r = report(d, a.dir / f"{d.name[:40]}_d0.csv", a.volume, a.cut_mm)
        if r:
            rows.append(r)
    if not rows:
        print("nothing to report")
        return

    rows.sort(key=lambda r: -r["duplicated_fraction"])
    print(f"physical self-overlap report   (wrap-scale cut = {a.cut_mm:.0f} mm)\n")
    print(f"{'trace':38s} {'area mm2':>9s} {'dup mm2':>9s} {'dup %':>7s} "
          f"{'extent mm':>10s} {'ev>cut':>7s}")
    print("-" * 86)
    for r in rows[:22]:
        print(f"{r['segment'][:36]:38s} {r['total_area_mm2']:9,.0f} "
              f"{r['duplicated_area_mm2']:9,.2f} "
              f"{r['duplicated_fraction']*100:6.3f}% "
              f"{r['overlap_extent_mm']:10.1f} {r['events_beyond_cut']:7,d}")

    lab = [r for r in rows if "-w0" in r["segment"]]
    ag = [r for r in rows if "auto_grown" in r["segment"]
          or "auto_trace" in r["segment"]]
    for name, rs in (("labelled", lab), ("multi-wrap", ag)):
        if not rs:
            continue
        d = np.array([r["duplicated_fraction"] for r in rs]) * 100
        e = np.array([r["overlap_extent_mm"] for r in rs])
        n = sum(1 for r in rs if r["events_beyond_cut"] > 0)
        print(f"\n{name} ({len(rs)}): duplicated area median {np.median(d):.3f}%, "
              f"max {d.max():.3f}%")
        print(f"   overlap extent median {np.median(e):.1f} mm, max {e.max():.1f} mm")
        print(f"   traces with wrap-scale overlap (>{a.cut_mm:.0f} mm): {n} of {len(rs)}")

    a.json.write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
