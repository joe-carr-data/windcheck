"""Command-line entry point for windcheck."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .catalog import fetch_catalog, load_segments, winding_run

DATA = Path(__file__).resolve().parents[2] / "data"


def cmd_status(args: argparse.Namespace) -> int:
    """Report the structure of a sample's published segmentation corpus."""
    catalog = fetch_catalog(DATA / "catalog.json.gz")
    segments = load_segments(catalog, args.sample)

    labelled = [s for s in segments if s.winding is not None]
    auto = [s for s in segments if s.winding is None]
    lo, hi, gaps = winding_run(labelled)

    print(f"sample {args.sample}: {len(segments)} segments")
    print(f"  labelled windings : {len(labelled)}  (w{lo:03d}-w{hi:03d})")
    print(f"  gaps in the run   : {gaps if gaps else 'none (contiguous)'}")
    print(f"  auto_grown traces : {len(auto)}")

    if labelled and auto:
        lu = sorted(s.grid[0] for s in labelled)
        au = sorted(s.grid[0] for s in auto)
        med_l, med_a = lu[len(lu) // 2], au[len(au) // 2]
        print(f"\n  median u-extent, labelled : {med_l}")
        print(f"  median u-extent, auto     : {med_a}  ({med_a / med_l:.1f}x)")
        print("  -> auto_grown traces span multiple wraps; single windings cannot")
        print("     sheet-switch, which is what makes them the null control.")

    if args.verbose:
        print(f"\n  {'wind':>5} {'grid(u,v)':>14}  long_id")
        for s in segments:
            w = f"w{s.winding:03d}" if s.winding is not None else "auto"
            print(f"  {w:>5} {list(s.grid)!s:>14}  {s.long_id}")
    return 0


def _discover(root: Path) -> dict[str, Path]:
    """Map segment name -> its tifxyz directory, for any published layout.

    Two layouts exist in the open-data bucket: `mesh/<id>-on-<volume>.tifxyz/`
    and `mesh/intermediate/tifxyz_original/`. Both are handled.
    """
    found: dict[str, Path] = {}
    for z in sorted(root.rglob("z.tif")):
        d = z.parent
        if not (d / "x.tif").exists():
            continue
        found.setdefault(d.relative_to(root).parts[0], d)
    # OBJ meshes are accepted too, but only where no tifxyz was found for that
    # segment: the grid form carries concentration statistics the mesh cannot.
    for o in sorted(root.rglob("*.obj")):
        found.setdefault(o.relative_to(root).parts[0] + ":" + o.stem, o)
    return found


def cmd_selfgap(args: argparse.Namespace) -> int:
    """Self-gap analysis over every trace under a directory."""
    from .selfgap import analyse

    root = Path(args.path)
    if not root.exists():
        print(f"no such path: {root}", file=sys.stderr)
        return 2

    if root.suffix.lower() == ".obj":
        targets = {root.stem: root}
    elif (root / "z.tif").exists():
        targets = {root.name: root}
    else:
        targets = _discover(root)
    if not targets:
        print(f"no tifxyz surfaces found under {root}", file=sys.stderr)
        return 2

    rows, skipped = [], 0
    # All percentages are of valid queries SUBMITTED, not of the ones that came
    # back finite. `cover` is how many came back finite, shown so the reader can
    # see how much of the trace the search radius actually reached.
    hdr = (f"{'trace':>34} {'u':>6} {'cover':>6} {'<2.0vx':>7} {'flag%':>7} "
           f"{'blob':>8} {'blob%':>7} {'top5':>6} {'du p10':>7} {'valid':>7}")
    print(hdr)
    print("-" * len(hdr))
    for name, path in targets.items():
        r = analyse(path, name=name, stride=args.stride, exclude_u=args.exclude_u)
        if r is None:
            skipped += 1
            continue
        rows.append(r.as_dict())
        print(f"{r.name[:34]:>34} {r.u_extent:>6} {r.coverage:>6.1%} "
              f"{r.frac_below_grower_th * 100:>7.3f} "
              f"{r.frac_flagged * 100:>7.3f} {r.largest_blob:>8,} "
              f"{r.blob_fraction * 100:>7.3f} {r.top5_share:>5.0%} "
              f"{r.du_p10:>7.0f} {'OK' if r.valid else 'REJECT':>7}")

    print(f"\n  {len(rows)} analysed, {skipped} skipped (single-wrap: no previous "
          f"wrap to compare against)")
    for r in rows:
        if not r["valid"]:
            print(f"  REJECT {r['name']}: {r['reason']}")

    if args.json:
        Path(args.json).write_text(json.dumps(rows, indent=1) + "\n")
        print(f"  wrote {args.json}")
    return 0


def cmd_calibrate(args: argparse.Namespace) -> int:
    """Re-derive the sheet-separation table quoted in the README.

    Measures point-to-surface distance from each labelled winding to the sheet
    `k` windings away. If the labels are a true radial index, the medians must
    come out proportional to `k` -- and they do, which is what licenses reading
    a ~2x gap as "one wrap was skipped".
    """
    import contextlib
    import io

    import numpy as np

    from . import atlas

    entries = atlas.discover(Path(args.path))
    refs = {e.winding: e for e in entries if e.winding is not None}
    if len(refs) < 4:
        print(f"need labelled wNNN segments; found {len(refs)}", file=sys.stderr)
        return 2

    work = Path("out")
    work.mkdir(exist_ok=True)
    per_k: dict[int, list] = {k: [] for k in range(1, args.max_k + 1)}
    for w in sorted(refs):
        pts, _ = atlas.sample_points(refs[w], stride=args.stride)
        atlas.write_queries(pts, work / "cal_q.bin")
        for k, bucket in per_k.items():
            if w + k not in refs:
                continue
            atlas.write_atlas([refs[w + k]], work / "cal_a.bin")
            with contextlib.redirect_stderr(io.StringIO()):
                r = atlas.run_engine(work / "cal_a.bin", work / "cal_q.bin",
                                     work / "cal_r.bin", threads=args.threads)
            d = r["d1"]
            bucket.append(d[np.isfinite(d)])

    print(f"sheet separation, {sum(len(v) for v in per_k.values())} surface pairs\n")
    print(f"{'gap':>4} {'pairs':>7} {'median vx':>10} {'ratio':>7}")
    base = None
    for k, chunks in per_k.items():
        if not chunks:
            continue
        a = np.concatenate(chunks)
        a = a[a < args.max_dist]        # drop non-overlapping tails
        med = float(np.median(a))
        base = base or med
        print(f"{k:>4} {len(chunks):>7} {med:>10.2f} {med / base:>7.3f}")
    print("\n  ratios must track the gap (1, 2, 3, ...) for the labels to be a")
    print("  radial index; that is what makes a doubled gap mean a skipped wrap.")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="windcheck",
        description="Cross-wrap consistency checking for scroll segmentations.",
    )
    p.add_argument("--version", action="version", version=f"windcheck {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("status", help="summarise a sample's segmentation corpus")
    s.add_argument("--sample", default="PHerc0172")
    s.add_argument("-v", "--verbose", action="store_true")
    s.set_defaults(func=cmd_status)

    g = sub.add_parser("selfgap", help="self-gap analysis of traced surfaces")
    g.add_argument("path", help="a tifxyz directory, or a directory of segments")
    g.add_argument("--stride", type=int, default=3, help="grid subsampling (default 3)")
    g.add_argument("--exclude-u", type=int, default=60,
                   help="columns excluded around each query's own u (default 60)")
    g.add_argument("--json", help="also write results as JSON to this path")
    g.set_defaults(func=cmd_selfgap)

    c = sub.add_parser("calibrate", help="re-derive the sheet-separation table")
    c.add_argument("path", help="directory of labelled wNNN segments")
    c.add_argument("--max-k", type=int, default=3, help="largest winding gap (default 3)")
    c.add_argument("--stride", type=int, default=8)
    c.add_argument("--max-dist", type=float, default=200.0,
                   help="ignore pairs beyond this, where the sheets do not overlap")
    c.add_argument("--threads", type=int, default=12)
    c.set_defaults(func=cmd_calibrate)

    args = p.parse_args(argv)
    return int(args.func(args))
