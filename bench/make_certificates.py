"""Emit a certificate and a VC3D-loadable overlay for a few example traces.

Two things this file gets right that an earlier version did not.

**One point per far-field event, grouped by the far-field components.** The
representatives used to be keyed off `lab`, the component map built from *all*
crossings, while the event count came from `l2`, built from the far ones only.
Two far regions joined by a chain of local crossings collapse into one key under
`lab`, so the overlay carried 357 points against 379 reported events. The counts
now come from the same map and an assert holds them together -- an overlay that
disagrees with the certificate beside it is worse than no overlay.

**A cut in revolutions, not millimetres.** The 30 mm cut held on Scroll 5 and
inverted on PHerc0139 and PHerc1667. See `docs/submission.md` section 4.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crossing_analyse import components, events, load           # noqa: E402
from revolution_diag import revolution_period                   # noqa: E402

from windcheck import tifxyz                                    # noqa: E402
from windcheck.certificate import (WRAP_SCALE_CUT_REV,          # noqa: E402
                                   certificate, write_collection)

SEGMENTS = [
    "20251115002745-auto_grown_20251115002740308_5_flatboi",
    "20251109232817-w065_20251109232817724_flatboi",
    "20250917143559-w062_20250917143559205_flatboi",
]
ROOT = Path("data/scroll5_tifxyz")
VOLUME = "20241024131839"
OUT = Path("out/certificates")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    phys = {r["segment"]: r
            for r in json.loads(Path("out/crossing/physical.json").read_text())}
    made = []

    for seg in SEGMENTS:
        mp = sorted((ROOT / seg).glob(f"mesh/*{VOLUME}*.tifxyz"))[0]
        s = tifxyz.read(mp)
        P = s.points
        rec = load(Path("out/crossing") / f"{seg[:40]}_d0.csv")
        p = phys.get(seg, {})
        pitch = p.get("pitch_um", 158.0)

        period = revolution_period(P, s.valid)
        ucols = np.nonzero(s.valid.any(0))[0]
        span_rev = (float(len(ucols)) / period) if period == period else None

        pts: list[list[float]] = []
        n_events = n_far = 0
        sep_rev = None
        if len(rec) and period == period:
            du = np.abs(rec["u1"].astype(np.int64) - rec["u2"].astype(np.int64))
            sep_rev = float(du.max()) / period
            lab, _ = components(rec)
            n_events = len(events(rec, lab))

            far = rec[du > WRAP_SCALE_CUT_REV * period]
            if len(far):
                l2, _ = components(far)
                n_far = len(events(far, l2))
                # Representative per event, from the SAME component map the
                # count came from: the pair with the deepest penetration.
                by_event: dict[tuple, np.void] = {}
                for r in far:
                    k = tuple(sorted((l2[(int(r["v1"]), int(r["u1"]))],
                                      l2[(int(r["v2"]), int(r["u2"]))])))
                    if k not in by_event or r["pen"] > by_event[k]["pen"]:
                        by_event[k] = r
                pts = [P[int(r["v1"]), int(r["u1"])].tolist()
                       for r in by_event.values()]
                assert len(pts) == n_far, (
                    f"{seg}: overlay has {len(pts)} points against {n_far} "
                    "events -- the two must be built from one component map")

        if pts:
            write_collection(
                OUT / f"{seg[:30]}_points.json",
                f"windcheck-overlaps-{seg[:22]}", pts,
                metadata={"tool": "windcheck",
                          "cut_revolutions": WRAP_SCALE_CUT_REV})

        cert = certificate(
            segment=seg, mesh_path=mp, voxel_um=p.get("voxel_um", 7.91),
            params={"exclude": 1, "cell": 40.0, "maxedge": 60.0,
                    "touch_tol": 1e-3, "diagonal": 0, "stride": 1},
            total_area_mm2=p.get("total_area_mm2", 0.0),
            participating_quad_area_mm2=p.get("duplicated_area_mm2", 0.0),
            parameter_span_mm_estimate=p.get("overlap_extent_mm", 0.0),
            separation_revolutions=sep_rev,
            covering_span_revolutions=span_rev,
            revolution_period_columns={
                "turning_estimate": (round(period, 1) if period == period
                                     else None),
                "method": "cumulative turning of the surface centreline; "
                          "cross-checked by bench/period_cross_check.py",
                "pitch_um": pitch,
            },
            n_pairs=len(rec), n_events=n_events, events_beyond_cut=n_far,
            median_penetration_vx=p.get("median_penetration_vx", 0.0))
        (OUT / f"{seg[:30]}_certificate.json").write_text(
            json.dumps(cert, indent=2))
        made.append((seg[-20:], len(pts), n_far, sep_rev, span_rev,
                     cert["verdict"]))

    for name, npts, nfar, sr, spr, verdict in made:
        srs = f"{sr:.3f}" if sr is not None else "  -  "
        sps = f"{spr:.2f}" if spr is not None else " - "
        print(f"  {name:22s} points {npts:4d}  events>cut {nfar:4d}  "
              f"sep {srs} rev  span {sps} rev  {verdict}")


if __name__ == "__main__":
    main()
