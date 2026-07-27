"""Positive control: does the period estimator recover a PUBLISHED winding count?

Every number in the revolution-normalised result so far rests on a period
measured from the surface's own geometry, cross-checked against a second
estimator that reads the same geometry. Two estimators agreeing is weaker
evidence than it sounds -- on PHerc1667 they agreed with each other and were
both wrong, because they had been handed the wrong file.

Scroll 1 (`PHercParis4`) supplies the missing external check. Its segments are
named by the range of windings they cover:

    20260701183132-w073-076   ->  windings 73 to 76, so about 4 revolutions
    20260701183146-w118-119   ->  windings 118 to 119, so about 2

That count comes from the publishers' own tracing, not from us and not from the
geometry we measure. So it is a genuine positive control: if the estimator says
`w073-076` covers 4.0 revolutions, the whole normalisation stands on measured
ground rather than on internal consistency.

The comparison is deliberately blunt -- predicted against declared, plus the
ratio -- because a control that needs tuning to pass is not a control.

    uv run python bench/winding_control.py --root data/scroll1_tifxyz
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from revolution_diag import turning_profile                    # noqa: E402

from windcheck import tifxyz                                   # noqa: E402

# e.g. "20260701183132-w073-076" -> (73, 76). Also accepts a single "w073".
NAME_RE = re.compile(r"-w(\d+)(?:-(\d+))?(?:[^0-9]|$)")


def declared_windings(name: str) -> int | None:
    """How many windings the publisher says this segment covers."""
    m = NAME_RE.search(name)
    if not m:
        return None
    lo = int(m.group(1))
    hi = int(m.group(2)) if m.group(2) else lo
    if hi < lo:
        return None
    return hi - lo + 1


def measured_revolutions(path: Path) -> float | None:
    s = tifxyz.read(path)
    theta = turning_profile(s.points, s.valid)
    if theta is None:
        return None
    ucols = np.nonzero(s.valid.any(0))[0]
    if len(ucols) < 8:
        return None
    return float(np.abs(theta[int(ucols[-1])] - theta[int(ucols[0])])) / (2 * np.pi)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/scroll1_tifxyz"))
    # The volume the submission quotes. Defaulting to the coarse 45.532um
    # variant produced 58 rows with three failures while the text described 33
    # rows on this one -- a reader following the default got different numbers
    # from the ones they were reading about.
    ap.add_argument("--volume", default="20230205180739")
    ap.add_argument("--json", type=Path, default=Path("out/winding_control.json"))
    a = ap.parse_args()

    rows = []
    print(f"{'segment':34s} {'declared':>8s} {'measured':>9s} {'ratio':>7s}")
    print("-" * 62)
    for d in sorted(a.root.iterdir()):
        if not d.is_dir():
            continue
        want = declared_windings(d.name)
        if want is None:
            continue
        m = sorted(d.glob(f"mesh/*{a.volume}*.tifxyz"))
        if not m:
            continue
        got = measured_revolutions(m[0])
        if got is None:
            print(f"{d.name[:32]:34s} {want:8d} {'-':>9s} {'-':>7s}")
            rows.append({"segment": d.name, "declared": want, "measured": None})
            continue
        rows.append({"segment": d.name, "declared": want,
                     "measured": round(got, 3),
                     "ratio": round(got / want, 3)})
        print(f"{d.name[:32]:34s} {want:8d} {got:9.3f} {got / want:7.3f}",
              flush=True)

    ok = [r for r in rows if r.get("measured") is not None]
    if ok:
        ratio = np.array([r["ratio"] for r in ok])
        dec = np.array([r["declared"] for r in ok])
        mea = np.array([r["measured"] for r in ok])
        print("-" * 62)
        print(f"{len(ok)} of {len(rows)} segments measurable")
        print(f"  ratio measured/declared: median {np.median(ratio):.3f}, "
              f"IQR {np.percentile(ratio, 25):.3f}-{np.percentile(ratio, 75):.3f}, "
              f"range {ratio.min():.3f}-{ratio.max():.3f}")
        print(f"  within 20% of declared: "
              f"{int((np.abs(ratio - 1) <= 0.2).sum())}/{len(ok)}")
        if len(set(dec.tolist())) > 1:
            r = np.corrcoef(dec, mea)[0, 1]
            print(f"  correlation declared vs measured: r = {r:.3f}")
        print()
        print("A median ratio near 1.0 means the estimator recovers the "
              "publishers' own\nwinding count from geometry alone. Anything "
              "else is the control failing,\nwhich is what it is for.")

    a.json.parent.mkdir(parents=True, exist_ok=True)
    a.json.write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
