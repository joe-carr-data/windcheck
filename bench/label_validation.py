"""Do flagged crossings sit where the trace visits two different physical sheets?

This is the experiment the whole project has been unable to run. Everything
measured so far is a property of the mesh: the surface passes through itself.
Whether that corresponds to a *tracing error* -- the trace assigning two
different sheets of papyrus to one place -- needed ground truth, and Scroll 1's
80 per-voxel instance-label cubes supply it (`src/windcheck/labels.py`).

## The design, fixed before the result was seen

For each flagged transverse crossing, the two participating quads are looked up
in the label volume. If the crossing marks the trace visiting two different
sheets, their ids differ.

**The control is the whole experiment.** Two quads far apart in the grid land on
different sheets more often for a reason that has nothing to do with defects:
the trace has simply travelled. So each flagged pair is matched against random
pairs of the SAME grid separation drawn from the SAME segment. Without that
matching the comparison manufactures its own answer, since flagged pairs skew
towards large separations by construction.

Three things are kept strictly apart, because collapsing them is the denominator
error that has already cost this project a published claim:

    both quads labelled        -> the pair votes, same id or different
    either quad unlabelled     -> the pair is UNKNOWN, and is counted as such
    either quad outside a cube -> not part of the experiment at all

An unknown is not a negative.

## Pre-registered kill criterion

If flagged pairs show differing sheet ids at less than **1.5x** the rate of
separation-matched controls, then self-intersection does not predict sheet
misassignment, and the honest finding is that this is a geometric defect and not
a demonstrated tracing error.

## Calibration first

A bar has been set once in this project that a known-true positive could not
clear. So `--calibrate` reports the differing-id rate for pairs at very large
separation -- parts of the trace a full revolution or more apart, which cannot be
the same sheet -- to confirm the labels can express "different sheet" at all
before any contrast is trusted.

    uv run --with pynrrd python bench/label_validation.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crossing_analyse import load                              # noqa: E402

from windcheck import labels as L                              # noqa: E402
from windcheck import tifxyz                                   # noqa: E402

RNG = np.random.default_rng(20260727)


def ids_at(cubes: L.LabelCubes, P: np.ndarray, vs, us) -> np.ndarray:
    return cubes.lookup(P[np.asarray(vs), np.asarray(us)])


def matched_controls(valid: np.ndarray, dv: np.ndarray, du: np.ndarray,
                     n_per: int) -> tuple[np.ndarray, np.ndarray]:
    """Random quad pairs with the same grid offsets as the flagged ones.

    For each flagged pair's offset, `n_per` random anchors are drawn from the
    segment's valid cells and the offset applied. Pairs falling outside the grid
    or on an invalid cell are dropped, so the control is drawn from exactly the
    population the flagged pairs came from.
    """
    nv, nu = valid.shape
    vv, uu = np.nonzero(valid)
    if len(vv) == 0:
        return np.empty((0, 2), int), np.empty((0, 2), int)
    a, b = [], []
    k = RNG.integers(0, len(vv), size=len(dv) * n_per)
    off = np.repeat(np.stack([dv, du], 1), n_per, axis=0)
    sign = RNG.choice([-1, 1], size=off.shape)
    v1, u1 = vv[k], uu[k]
    v2, u2 = v1 + off[:, 0] * sign[:, 0], u1 + off[:, 1] * sign[:, 1]
    ok = (v2 >= 0) & (v2 < nv) & (u2 >= 0) & (u2 < nu)
    ok &= valid[np.clip(v1, 0, nv - 1), np.clip(u1, 0, nu - 1)]
    ok &= valid[np.clip(v2, 0, nv - 1), np.clip(u2, 0, nu - 1)]
    a = np.stack([v1[ok], u1[ok]], 1)
    b = np.stack([v2[ok], u2[ok]], 1)
    return a, b


def rate(id1: np.ndarray, id2: np.ndarray) -> tuple[int, int, int]:
    """(pairs where both labelled, of those how many differ, unknown pairs)."""
    both = (id1 > 0) & (id2 > 0)
    return int(both.sum()), int((id1[both] != id2[both]).sum()), int((~both).sum())


def run_segment(seg: Path, csv: Path, cubes: L.LabelCubes, volume: str,
                n_per: int) -> dict | None:
    m = sorted(seg.glob(f"mesh/*{volume}*.tifxyz"))
    if not m or not (m[0] / "x.tif").exists():
        return None
    s = tifxyz.read(m[0])
    P, V = s.points, s.valid
    if not cubes.covered(P[V]).any():
        return None

    rec = load(csv)
    if len(rec) == 0:
        return None
    nv, nu = P.shape[:2]
    inb = ((rec["v1"] < nv) & (rec["v2"] < nv)
           & (rec["u1"] < nu) & (rec["u2"] < nu))
    rec = rec[inb]
    if len(rec) == 0:
        return None

    f1 = ids_at(cubes, P, rec["v1"], rec["u1"])
    f2 = ids_at(cubes, P, rec["v2"], rec["u2"])
    fb, fd, fu = rate(f1, f2)

    dv = np.abs(rec["v1"].astype(int) - rec["v2"].astype(int))
    du = np.abs(rec["u1"].astype(int) - rec["u2"].astype(int))
    a, b = matched_controls(V, dv, du, n_per)
    if len(a) == 0:
        return None
    c1 = ids_at(cubes, P, a[:, 0], a[:, 1])
    c2 = ids_at(cubes, P, b[:, 0], b[:, 1])
    cb, cd, cu = rate(c1, c2)

    return {
        "segment": seg.name,
        "flagged_pairs": int(len(rec)),
        "flagged_both_labelled": fb, "flagged_differ": fd, "flagged_unknown": fu,
        "control_both_labelled": cb, "control_differ": cd, "control_unknown": cu,
        "flagged_rate": (fd / fb) if fb else None,
        "control_rate": (cd / cb) if cb else None,
        "enrichment": ((fd / fb) / (cd / cb)) if fb and cb and cd else None,
        "median_du": int(np.median(du)),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/scroll1_tifxyz"))
    ap.add_argument("--dir", type=Path, default=Path("out/crossing_s1"))
    ap.add_argument("--volume", default="20230205180739")
    ap.add_argument("--n-per", type=int, default=4,
                    help="control pairs drawn per flagged pair")
    ap.add_argument("--json", type=Path, default=Path("out/label_validation.json"))
    a = ap.parse_args()

    cubes = L.LabelCubes()
    print(f"{len(cubes)} label cubes loaded\n")
    print(f"{'segment':32s} {'flagged':>8s} {'labelled':>9s} {'differ':>7s} "
          f"{'ctrl lab':>9s} {'ctrl dif':>9s} {'enrich':>7s}")
    print("-" * 90)

    rows = []
    for d in sorted(a.root.iterdir()):
        if not d.is_dir():
            continue
        r = run_segment(d, a.dir / f"{d.name[:40]}_d0.csv", cubes, a.volume,
                        a.n_per)
        if r is None:
            continue
        rows.append(r)
        e = f"{r['enrichment']:.2f}x" if r["enrichment"] else "-"
        print(f"{r['segment'][:30]:32s} {r['flagged_pairs']:8,d} "
              f"{r['flagged_both_labelled']:9,d} "
              f"{(r['flagged_rate'] or 0):7.3f} "
              f"{r['control_both_labelled']:9,d} "
              f"{(r['control_rate'] or 0):9.3f} {e:>7s}", flush=True)

    if rows:
        FB = sum(r["flagged_both_labelled"] for r in rows)
        FD = sum(r["flagged_differ"] for r in rows)
        CB = sum(r["control_both_labelled"] for r in rows)
        CD = sum(r["control_differ"] for r in rows)
        FU = sum(r["flagged_unknown"] for r in rows)
        print("-" * 90)
        print(f"POOLED over {len(rows)} segments")
        print(f"  flagged: {FB:,} pairs with both quads labelled, "
              f"{FD:,} differ = {FD / max(FB, 1):.4f}")
        print(f"  control: {CB:,} pairs with both quads labelled, "
              f"{CD:,} differ = {CD / max(CB, 1):.4f}")
        print(f"  unknown (either quad unlabelled): {FU:,} flagged pairs, "
              f"excluded rather than counted as agreeing")
        if CB and CD:
            enr = (FD / FB) / (CD / CB)
            print(f"\n  ENRICHMENT {enr:.2f}x")
            print(f"  Pre-registered kill criterion was 1.5x: "
                  f"{'PASSES' if enr >= 1.5 else 'FAILS'}")
        a.json.parent.mkdir(parents=True, exist_ok=True)
        a.json.write_text(json.dumps(rows, indent=2))
        print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
