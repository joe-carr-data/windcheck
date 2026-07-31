"""The four-corpus comparison, in revolutions instead of millimetres.

Why this file exists. The published headline was a 17.2 mm / 135.7 mm gap
between Scroll 5's labelled windings and its multi-wrap traces. A control on two
further scrolls fired before publication: PHerc0139's labelled windings reach
89 mm and PHerc1667's reach 1237 mm, so an absolute cut in millimetres puts
them on the wrong side. The criterion was a Scroll 5 observation, not a
constant.

The fix is not a bigger threshold. It is a change of unit and a change of
statement.

**Unit.** Separation is expressed as a fraction of the segment's own revolution
period, measured from its own geometry (`revolution_diag.revolution_period`,
cross-checked against `period_cross_check.period_from_neighbour`). Millimetres
are not comparable between scrolls: the voxel size differs, the scrolls differ
in diameter, and a scroll that is crushed has no single circumference.

**Statement.** Not "a trace with separation above X is defective", but a
conditional whose condition is measured rather than assumed:

    among segments that cover at most one revolution, none shows a crossing
    beyond a small fraction of a revolution

That is testable on every corpus, because the covering span is measured from
the same geometry as the separation. It also names the population the earlier
claim silently assumed: "labelled winding" is a filename convention, and only
on Scroll 5 does it happen to coincide with one revolution.

Reads the `revdiag.json` written by `revolution_diag.py` for each corpus.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

# Stratum boundaries, in revolutions of covering span.
#
# 1.005 is a 0.5% tolerance on "covers one revolution", NOT an allowance for
# estimator error -- an earlier comment here claimed it absorbed ~15%, which it
# plainly does not, and a segment measured at 1.032 does carry a 0.944
# separation. Because the period is estimated, a segment near the boundary can
# cross it under a small period error, so `stratum_sensitivity` widens the
# margin on both sides and shows what survives. The headline is quoted from the
# conservative strata, not from these.
SINGLE_MAX = 1.005
MULTI_MIN = 2.0

CORPORA = [
    ("Scroll 1", "out/crossing_s1/revdiag.json", "out/crossing_s1/period.json"),
    ("Scroll 5", "out/crossing/revdiag.json", "out/crossing/period.json"),
    ("PHerc0814", "out/crossing_0814/revdiag.json", "out/crossing_0814/period.json"),
    ("PHerc0139", "out/crossing_0139/revdiag.json", "out/crossing_0139/period.json"),
    ("PHerc1667", "out/crossing_1667/revdiag.json", "out/crossing_1667/period.json"),
]

# Acceptance band for the two period estimators, fixed before looking at what it
# excludes: two geometric estimates of the same quantity should agree within a
# quarter. Segments outside it are dropped and COUNTED -- an earlier version of
# this file described a cross-check it never actually loaded, which is exactly
# the sort of claim that has to be either implemented or deleted.
#
# The cross-check is unavailable on a segment with no neighbouring wrap inside
# itself, which is most single-covering segments: they are the null control, and
# having no second sheet to find is the whole point of them. Coverage is
# reported rather than glossed.
RATIO_LO, RATIO_HI = 0.8, 1.25
MIN_SUPPORT = 200


def cross_check(path: Path) -> dict[str, dict]:
    """Neighbour-period estimates by segment, from period_cross_check.py."""
    p = Path(path)
    if not p.exists():
        return {}
    return {r["segment"]: r for r in json.loads(p.read_text())}


def load(path: Path, xcheck: Path | None = None) -> list[dict]:
    """Segments with a usable period, in column-offset units.

    Separation is a column offset divided by the measured revolution period.
    The obvious-looking alternative -- read the angle straight off the turning
    profile and divide by 2*pi -- was built and TESTED, and it is worse. It is
    noise-dominated at short range: on six Scroll 5 single-covering segments,
    pairs only 14 to 78 columns apart register up to 0.26 revolutions, which is
    94 degrees. Two quads fourteen columns apart are not 84 degrees apart around
    a scroll.

    The cause is geometric rather than incidental. A segment covering less than
    one revolution is a short arc, so its own centroid sits near the arc itself
    and the projected points pass close to the origin, where `arctan2` swings
    freely. The angle is well-conditioned only over long spans.

    So the angular figures are kept in the JSON and used for two things they are
    good at -- confirming that a large separation is not an artifact of
    non-uniform column spacing, and the self-consistency gate in
    `revolution_diag` -- but not as the primary measure. This is not a choice of
    the flattering number: the angular measure makes a claim about 14-column
    pairs that is physically false, which is a falsification, not a preference.

    A span under a fifth of a revolution means the turning integral had almost
    no angle to work with, so its period is not a period. Dropped rather than
    trusted, and counted so the drop is visible.
    """
    rows = json.loads(Path(path).read_text())
    xc = cross_check(xcheck) if xcheck else {}
    out = []
    for r in rows:
        spr, dpr = r.get("span_per_rev"), r.get("du_max_per_rev")
        if spr is None or dpr is None or spr < 0.2:
            continue
        c = xc.get(r["segment"])
        if c and c.get("n_support", 0) >= MIN_SUPPORT:
            ratio = c.get("ratio")
            if ratio is None or ratio != ratio:
                r = {**r, "xcheck": "unavailable"}
            elif not (RATIO_LO <= ratio <= RATIO_HI):
                r = {**r, "xcheck": "rejected", "xcheck_ratio": ratio}
                out.append(r)
                continue
            else:
                r = {**r, "xcheck": "agreed", "xcheck_ratio": ratio}
        else:
            r = {**r, "xcheck": "unavailable"}
        out.append(r)
    return out


def accepted(rows):
    return [r for r in rows if r.get("xcheck") != "rejected"]


def band(rows, lo, hi):
    return [r for r in rows if lo <= r["span_per_rev"] < hi]


def fmt(rows):
    if not rows:
        return f"{'-':>5s} {'':>25s}"
    d = np.array([r["du_max_per_rev"] for r in rows])
    n_any = int((d > 0.15).sum())
    return (f"{len(rows):5d}  max {d.max():6.3f}  median {np.median(d):6.3f}  "
            f">0.15 rev: {n_any:3d}/{len(rows)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", type=Path,
                    default=Path("out/revolution_summary.json"))
    a = ap.parse_args()

    print("Self-overlap separation as a fraction of the segment's own "
          "revolution period\n")
    print(f"{'corpus':12s} {'covering span':22s} separation du_max / revolution")
    print("-" * 92)

    summary = []
    for name, path, xc in CORPORA:
        p = Path(path)
        if not p.exists():
            print(f"{name:12s} (no revdiag.json -- run revolution_diag.py)")
            continue
        allrows = load(p, Path(xc) if xc else None)
        rows = accepted(allrows)
        raw = len(json.loads(p.read_text()))
        bands = [(f"<= 1 revolution", band(rows, 0.0, SINGLE_MAX)),
                 (f"1 - 2 revolutions", band(rows, SINGLE_MAX, MULTI_MIN)),
                 (f"> 2 revolutions", band(rows, MULTI_MIN, 1e9))]
        for i, (label, rs) in enumerate(bands):
            head = name if i == 0 else ""
            print(f"{head:12s} {label:22s} {fmt(rs)}")
            summary.append({"corpus": name, "band": label, "n": len(rs),
                            "max_du_per_rev": (max(r["du_max_per_rev"] for r in rs)
                                               if rs else None),
                            "segments_above_0.15": [r["segment"] for r in rs
                                                    if r["du_max_per_rev"] > 0.15]})
        nrej = sum(1 for r in allrows if r.get("xcheck") == "rejected")
        nagree = sum(1 for r in rows if r.get("xcheck") == "agreed")
        drop_span = raw - len(allrows)
        if drop_span:
            print(f"{'':12s} ({drop_span} dropped: span under 0.2 revolution, "
                  f"period not measurable)")
        print(f"{'':12s} cross-check: {nagree} agreed, {nrej} rejected "
              f"(ratio outside [{RATIO_LO}, {RATIO_HI}]), "
              f"{len(allrows) - nagree - nrej} unavailable")
        print()

    pooled = [r for _, p, x in CORPORA if Path(p).exists()
              for r in accepted(load(Path(p), Path(x)))]
    single = band(pooled, 0.0, SINGLE_MAX)
    multi = band(pooled, MULTI_MIN, 1e9)
    unfiltered = [r for _, p, x in CORPORA if Path(p).exists()
                  for r in load(Path(p), Path(x))]
    print()
    strata("POOLED, no cross-check exclusion at all (the conservative headline)",
           unfiltered)
    strata(f"POOLED, cross-check enforced at [{RATIO_LO}, {RATIO_HI}]", pooled)
    sensitivity(unfiltered)
    stratum_sensitivity(unfiltered)
    bands = band_structure(pooled)
    summary.append({"band_structure": bands})

    a.json.parent.mkdir(parents=True, exist_ok=True)
    a.json.write_text(json.dumps(summary, indent=2))
    print(f"\nwrote {a.json}")


def rank_auc(lo_group, hi_group) -> float:
    """Probability a randomly drawn long-covering segment separates further.

    Reported alongside the extreme-value gap because a gap is two numbers and
    an AUC is all of them. A gap can be produced by one lucky pair of extremes;
    an AUC near 1 cannot.
    """
    a = np.array([r["du_max_per_rev"] for r in lo_group])
    b = np.array([r["du_max_per_rev"] for r in hi_group])
    if not len(a) or not len(b):
        return float("nan")
    wins = (b[:, None] > a[None, :]).sum() + 0.5 * (b[:, None] == a[None, :]).sum()
    return float(wins) / (len(a) * len(b))


def strata(title: str, rows) -> None:
    """The two exposure-conditioned strata, with the gap and the AUC."""
    lo, hi = band(rows, 0.0, SINGLE_MAX), band(rows, MULTI_MIN, 1e9)
    if not lo or not hi:
        return
    s = np.array([r["du_max_per_rev"] for r in lo])
    m = np.array([r["du_max_per_rev"] for r in hi])
    m_nz = m[m > 0]
    print("=" * 92)
    print(title)
    print(f"  estimated <=1-period stratum   n={len(s):3d}   "
          f"median {np.median(s):.3f}  max {s.max():.3f} rev")
    print(f"  estimated >2-period stratum    n={len(m):3d}   "
          f"median {np.median(m):.3f}  min nonzero {m_nz.min():.3f} rev "
          f"({int((m == 0).sum())} with no crossing)")
    if s.max() < m_nz.min():
        print(f"  GAP {s.max():.3f} <-> {m_nz.min():.3f} rev = "
              f"{m_nz.min() / max(s.max(), 1e-9):.1f}x, nothing between")
    else:
        print("  NO GAP: the strata overlap in these units")
    print(f"  rank AUC {rank_auc(lo, hi):.3f}")


def sensitivity(unfiltered) -> None:
    """Does the headline depend on where the acceptance band is drawn?

    It is the question to ask of any threshold that was not fixed in advance,
    and the honest way to answer it is to vary the threshold and show the table
    rather than to defend one value.
    """
    global RATIO_LO, RATIO_HI
    keep = (RATIO_LO, RATIO_HI)
    print("=" * 92)
    print("SENSITIVITY to the cross-check acceptance band")
    print(f"  {'band':>12s} {'rejected':>9s} {'<=1 n':>6s} {'max':>7s} "
          f"{'>2 n':>5s} {'min nz':>7s} {'gap':>7s}")
    for lo, hi in ((0.0, 1e9), (0.9, 1.11), (0.85, 1.18),
                   (0.8, 1.25), (0.7, 1.43), (0.5, 2.0)):
        RATIO_LO, RATIO_HI = lo, hi
        rows = [r for _, p, x in CORPORA if Path(p).exists()
                for r in accepted(load(Path(p), Path(x)))]
        nrej = sum(1 for r in unfiltered
                   if r.get("xcheck_ratio") is not None
                   and not (lo <= r["xcheck_ratio"] <= hi))
        a = np.array([r["du_max_per_rev"] for r in band(rows, 0.0, SINGLE_MAX)])
        b = np.array([r["du_max_per_rev"] for r in band(rows, MULTI_MIN, 1e9)])
        bnz = b[b > 0]
        lab = "none" if hi > 100 else f"[{lo},{hi}]"
        print(f"  {lab:>12s} {nrej:9d} {len(a):6d} {a.max():7.3f} "
              f"{len(b):5d} {bnz.min():7.3f} "
              f"{bnz.min() / max(a.max(), 1e-9):6.1f}x")
    RATIO_LO, RATIO_HI = keep


def stratum_sensitivity(rows) -> None:
    """Widen the margin around the stratum boundaries and see what survives.

    The period is estimated, so a segment measured at 1.03 revolutions might
    truly cover 0.98. Pulling the boundaries apart discards the segments whose
    membership a plausible period error could flip, and keeps only those no
    error of that size could move. If the gap holds under the widest margin,
    the boundary was not doing the work.
    """
    global SINGLE_MAX, MULTI_MIN
    keep = (SINGLE_MAX, MULTI_MIN)
    print("=" * 92)
    print("SENSITIVITY to the stratum boundaries (no cross-check exclusion)")
    print(f"  {'<= x rev':>10s} {'>= y rev':>9s} {'<=1 n':>6s} {'max':>7s} "
          f"{'>2 n':>5s} {'min nz':>7s} {'gap':>7s} {'AUC':>6s}")
    for lo_b, hi_b in ((1.005, 2.0), (0.95, 2.2), (0.90, 2.35), (0.85, 2.48)):
        SINGLE_MAX, MULTI_MIN = lo_b, hi_b
        a_rows, b_rows = band(rows, 0.0, lo_b), band(rows, hi_b, 1e9)
        if not a_rows or not b_rows:
            continue
        a = np.array([r["du_max_per_rev"] for r in a_rows])
        b = np.array([r["du_max_per_rev"] for r in b_rows])
        bnz = b[b > 0]
        print(f"  {lo_b:10.3f} {hi_b:9.2f} {len(a):6d} {a.max():7.3f} "
              f"{len(b):5d} {bnz.min():7.3f} "
              f"{bnz.min() / max(a.max(), 1e-9):6.1f}x "
              f"{rank_auc(a_rows, b_rows):6.3f}")
    SINGLE_MAX, MULTI_MIN = keep


def band_structure(rows, min_ratio: float = 1.4,
                   min_gap_rev: float = 0.05) -> dict:
    """Find the bands in the pooled distribution rather than imposing them.

    The earlier version of this result chose a cut and defended it. This asks
    the data where its gaps are: sort every segment's maximum separation, take
    the ratio between consecutive values, and split wherever that ratio exceeds
    `min_ratio`. A threshold that has to be argued for is a weakness; a gap that
    is simply the largest feature in the distribution is not.

    A split also has to clear `min_gap_rev` in absolute terms. Ratios alone
    split noise: 0.002 against 0.003 revolutions is a 1.5x ratio and a
    difference of a fifth of a degree. 0.05 revolutions is 18 degrees, which is
    the smallest gap that means anything on a scroll.
    """
    vals = sorted((r["du_max_per_rev"], r["span_per_rev"], r["segment"])
                  for r in rows if r["du_max_per_rev"] > 0)
    if len(vals) < 3:
        return {}
    cuts = [i for i in range(len(vals) - 1)
            if vals[i + 1][0] / max(vals[i][0], 1e-9) >= min_ratio
            and vals[i + 1][0] - vals[i][0] >= min_gap_rev]

    print("=" * 92)
    print(f"BAND STRUCTURE, pooled -- splits where consecutive maxima differ by "
          f"{min_ratio}x AND {min_gap_rev} rev or more")
    print(f"  {len(vals)} segments show a crossing; "
          f"{len(rows) - len(vals)} show none at all")
    out, lo = [], 0
    for hi in cuts + [len(vals) - 1]:
        seg = vals[lo:hi + 1]
        d = [s[0] for s in seg]
        sp = [s[1] for s in seg]
        entry = {"n": len(seg), "sep_min": d[0], "sep_max": d[-1],
                 "covering_span_min": min(sp), "covering_span_max": max(sp)}
        out.append(entry)
        ratio = (f"  <-- gap {vals[hi + 1][0] / max(d[-1], 1e-9):5.2f}x"
                 if hi + 1 < len(vals) else "")
        print(f"  n={len(seg):3d}   separation {d[0]:6.3f} - {d[-1]:6.3f} rev   "
              f"covering span {min(sp):5.2f} - {max(sp):5.2f} rev{ratio}")
        lo = hi + 1
    ratios = sorted(((vals[i + 1][0] / max(vals[i][0], 1e-9), vals[i][0],
                      vals[i + 1][0])
                     for i in range(len(vals) - 1)
                     if vals[i + 1][0] - vals[i][0] >= min_gap_rev),
                    reverse=True)
    if ratios:
        print("  gaps by size:")
        for r, lo_v, hi_v in ratios[:4]:
            print(f"    {lo_v:6.3f} -> {hi_v:6.3f} rev   {r:5.2f}x")
    return {"min_ratio": min_ratio, "min_gap_rev": min_gap_rev, "bands": out,
            "gaps": [{"ratio": round(r, 2), "below": lo_v, "above": hi_v}
                     for r, lo_v, hi_v in ratios[:4]]}


if __name__ == "__main__":
    main()
