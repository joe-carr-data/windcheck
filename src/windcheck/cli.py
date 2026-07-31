"""Command-line entry point for windcheck.

Two commands do the work, and the split between them is the whole point:

    windcheck check <tifxyz>              report only; never writes a mesh
    windcheck transform <tifxyz> --out D  the explicit change to geometry

`check` is the default posture. It reads a segment, runs the both-diagonal
non-adjacent transverse-contact census on it, prints a verdict, and leaves
behind a certificate and a point collection the viewer can open. It does not
modify the input and it does not produce a transformed mesh -- there is no flag
on `check` that makes it do either. Transformation happens only when the user
types `transform`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .catalog import fetch_catalog, load_segments, winding_run

DATA = Path(__file__).resolve().parents[2] / "data"


def cmd_check(args: argparse.Namespace) -> int:
    """Report-only census of one segment."""
    from . import pipeline

    try:
        r = pipeline.check_segment(
            Path(args.path), Path(args.out), volume=args.volume,
            threads=args.threads, cell=args.cell, maxedge=args.maxedge)
    except (FileNotFoundError, pipeline.CensusRefusal) as e:
        print(f"windcheck: {e}", file=sys.stderr)
        return 2

    print(f"\n{r['segment']}")
    print(f"  mesh                {r['mesh']}")
    print(f"  grid                {r['grid'][0]} x {r['grid'][1]}"
          f"{'':6s}triangles  {r['triangles']:,}")
    print("  transverse contacts "
          f"diagonal 0 = {r['transverse_d0']:,}   "
          f"diagonal 1 = {r['transverse_d1']:,}")
    print("  coplanar contacts   "
          f"diagonal 0 = {r['coplanar_d0']:,}   "
          f"diagonal 1 = {r['coplanar_d1']:,}")
    print(f"  crossing events     {r['events']:,}")
    print(f"  VERDICT             {r['verdict']}")
    print(f"  elapsed             {r['wall_seconds']:.1f} s")
    print(f"\n  certificate  {r['certificate']}")
    print(f"  overlay      {r['points']}   <- open in VC3D "
          f"({r['n_points']} points)")
    print("\n  Report only: the input was not modified and no mesh was "
          "written.\n  Run `windcheck transform` to change geometry.")
    return 0


def cmd_transform(args: argparse.Namespace) -> int:
    """Explicitly transform one segment under the frozen policy."""
    from . import pipeline

    try:
        r = pipeline.transform_segment(
            Path(args.path), Path(args.out), volume=args.volume,
            threads=args.threads, cell=args.cell, maxedge=args.maxedge,
            displaced=(Path(args.displaced) if args.displaced else None),
            use_displaced=not args.no_displacement)
    except (FileNotFoundError, ValueError, pipeline.CensusRefusal) as e:
        print(f"windcheck: {e}", file=sys.stderr)
        return 2

    print(f"\n{r['segment']}")
    if r.get("before"):
        b, a = r["before"], r.get("after") or r["before"]
        print(f"  before              transverse d0 = {b['d0']['transverse']:,}"
              f"   d1 = {b['d1']['transverse']:,}")
        print(f"  after               transverse d0 = {a['d0']['transverse']:,}"
              f"   d1 = {a['d1']['transverse']:,}")
    if r.get("displacement"):
        print("  displacement        applied (bounded displacement base)")
    if r.get("removed_quads") is not None:
        print(f"  removed quads       {r['removed_quads']:,}")
    if r.get("retained_fraction") is not None:
        print(f"  retained area       {r['retained_fraction'] * 100:.4f}%")
    print(f"  STATUS              {r['status']}")
    print(f"  elapsed             {r['wall_seconds']:.1f} s")
    if r.get("mesh_out"):
        print(f"\n  mesh         {r['mesh_out']}")
    print(f"  certificate  {r['certificate']}")
    if not r.get("clean"):
        print("\n  No clean claim is made for this output.")
    return 0 if r.get("clean") else 1


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

    if args.verbose:
        print(f"\n  {'wind':>5} {'grid(u,v)':>14}  long_id")
        for s in segments:
            w = f"w{s.winding:03d}" if s.winding is not None else "auto"
            print(f"  {w:>5} {list(s.grid)!s:>14}  {s.long_id}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="windcheck",
        description="Find where a traced scroll surface passes through "
                    "itself, and -- only when asked -- repair it.")
    p.add_argument("--version", action="version",
                   version=f"windcheck {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    k = sub.add_parser(
        "check", help="report-only: census one surface for self-intersection",
        description="Run the both-diagonal non-adjacent transverse-contact "
                    "census on one segment and report. Never modifies the "
                    "input; never writes a mesh.")
    k.add_argument("path", help="a .tifxyz directory, or a segment directory")
    k.add_argument("--out", default="out/check",
                   help="where the certificate and overlay go "
                        "(default out/check)")
    k.add_argument("--volume", default="",
                   help="substring selecting which published volume to read "
                        "when the path is a segment directory")
    k.add_argument("--threads", type=int, default=0, help="0 = all cores")
    k.add_argument("--cell", type=float, default=40.0,
                   help="broad-phase grid cell in voxels (default 40)")
    k.add_argument("--maxedge", type=float, default=60.0,
                   help="drop quads with any edge longer than this (voxels)")
    k.set_defaults(func=cmd_check)

    t = sub.add_parser(
        "transform", help="explicitly repair one surface and emit a new mesh",
        description="Bounded displacement where applicable, then one "
                    "certified excision of every residual transverse contact, "
                    "under the frozen scheduling policy. Emits one aggregate "
                    ".tifxyz plus its certificate.")
    t.add_argument("path", help="a .tifxyz directory, or a segment directory")
    t.add_argument("--out", required=True,
                   help="output directory for the mesh and its certificate")
    t.add_argument("--volume", default="",
                   help="substring selecting which published volume to read "
                        "when the path is a segment directory")
    t.add_argument("--threads", type=int, default=0, help="0 = all cores")
    t.add_argument("--cell", type=float, default=40.0)
    t.add_argument("--maxedge", type=float, default=60.0,
                   help="drop quads with any edge longer than this (voxels)")
    t.add_argument("--displaced", default=None,
                   help="the bounded-displacement repair of this segment to "
                        "cut instead of the published mesh")
    t.add_argument("--no-displacement", action="store_true",
                   help="cut the published mesh even when a displacement "
                        "repair for this segment is on disk")
    t.set_defaults(func=cmd_transform)

    s = sub.add_parser("status", help="summarise a sample's segmentation corpus")
    s.add_argument("--sample", default="PHerc0172")
    s.add_argument("-v", "--verbose", action="store_true")
    s.set_defaults(func=cmd_status)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))
