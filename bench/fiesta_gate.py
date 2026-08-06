"""Fiesta gate, day 2-3: GGEdge group-submesh census (prereg Addendum 1).

The pre-registered primary analysis for the Scroll Fiesta experiment
(notes/PREREG-FIESTA-GATE.md + Addendum 1): the solver's REAL merge
constraints are GGEdges over (cube, winding-group) nodes, dumped by the
instrumented scroll_whole (SF_GG_DUMP). This script

  attribute  maps every placed-mesh face to its winding group and reports
             the attribution accounting per cube;
  gate       builds ONE triangle-soup atlas over every (cube, gid) that
             appears as a GGEdge endpoint, runs the C++ census ONCE (one
             diagonal pass -- soups carry their triangulation; all pair
             work on the engine's thread pool), and joins cross-surface
             verdicts back onto the edge list;
  spotcheck  re-reads sampled triangles from the raw OBJ text through an
             independent parser and compares them to what the atlas was
             fed (the prereg's geometry-preserving conversion check).

FROZEN ATTRIBUTION RULE (recorded before any gate measurement was made):
a face belongs to group g iff ALL THREE of its vertices carry gid == g
with g >= 0. Faces with mixed or unmapped (-1) vertex gids belong to no
group; they are counted and reported, never silently dropped. Vertex i of
<id>_placed.obj corresponds to group.i32[i]: PlacedCube_finalize writes
placed.obj from mesh.obj's full vertex array in order (ObjIO_write_uv_
masked writes all nv vertices; facekeep only filters faces), and
group.i32 is written per mesh vertex in Pass B. OBJ vertex lines are
"v z y x" (their CLI convention) and are swapped to (x, y, z) on load.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
for _p in (REPO_ROOT / "src",):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from windcheck import atlas as wc_atlas                          # noqa: E402
from windcheck import objmesh                                    # noqa: E402
from windcheck.pipeline import ENGINE                            # noqa: E402

# Engine invocation, identical shape to the soup parity test:
# <atlas> <out.csv> <threads> <cell> <exclude> <diagonal> <maxedge>
CELL = "40.0"
EXCLUDE = "1"
DIAGONAL = "0"
MAXEDGE = "60.0"


# ------------------------------------------------------------------ dump

def read_gg_dump(path: Path) -> tuple[list[dict], list[dict]]:
    """Parse the SF_GG_DUMP TSV: node and edge rows."""
    nodes, edges = [], []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if f[0] == "node":
            nodes.append({"id": int(f[1]), "cube_idx": int(f[2]),
                          "cube": f[3], "gid": int(f[4]), "k": int(f[5]),
                          "du": float(f[6]), "prior": float(f[7]),
                          "n_prior": int(f[8]), "comp": int(f[9])})
        elif f[0] == "edge":
            edges.append({"a": int(f[1]), "b": int(f[2]),
                          "cube_a": f[3], "gid_a": int(f[4]),
                          "cube_b": f[5], "gid_b": int(f[6]),
                          "obs": int(f[7]), "du": float(f[8]),
                          "conf": float(f[9]), "mad": float(f[10]),
                          "n": int(f[11]), "in_tree": int(f[12])})
    if not nodes or not edges:
        raise SystemExit(f"{path}: no nodes/edges parsed")
    return nodes, edges


# ------------------------------------------------------------ attribution

def load_cube(placed: Path, cube: str) -> tuple[objmesh.Mesh, np.ndarray]:
    """A cube's placed mesh (verts swapped zyx -> xyz) and per-vertex gid."""
    mesh = objmesh.read(placed / f"{cube}_placed.obj")
    gid = np.fromfile(placed / f"{cube}_group.i32", dtype="<i4")
    if len(gid) != len(mesh.verts):
        raise SystemExit(f"{cube}: group.i32 has {len(gid)} entries but "
                         f"placed.obj has {len(mesh.verts)} vertices -- the "
                         "per-vertex correspondence does not hold")
    mesh.verts = mesh.verts[:, ::-1].copy()      # file order is z y x
    return mesh, gid

def face_groups(mesh: objmesh.Mesh, gid: np.ndarray) -> tuple[np.ndarray, dict]:
    """Per-face group id under the frozen rule; -2 marks mixed/unmapped."""
    g3 = gid[mesh.tris]                          # (T, 3)
    same = (g3[:, 0] == g3[:, 1]) & (g3[:, 0] == g3[:, 2]) & (g3[:, 0] >= 0)
    fg = np.where(same, g3[:, 0], -2).astype(np.int32)
    stats = {"n_faces": int(len(fg)),
             "n_attributed": int(same.sum()),
             "n_mixed_or_unmapped": int((~same).sum()),
             "groups": {int(g): int(c) for g, c in
                        zip(*np.unique(fg[same], return_counts=True))}}
    return fg, stats


def cmd_attribute(args) -> int:
    placed = Path(args.placed)
    out = {}
    for objp in sorted(placed.glob("*_placed.obj")):
        cube = objp.name[:-len("_placed.obj")]
        mesh, gid = load_cube(placed, cube)
        _, stats = face_groups(mesh, gid)
        stats["gid_values_in_group_i32"] = sorted(
            int(g) for g in np.unique(gid))
        out[cube] = stats
    tot = Counter()
    for s in out.values():
        tot["faces"] += s["n_faces"]
        tot["attributed"] += s["n_attributed"]
        tot["mixed_or_unmapped"] += s["n_mixed_or_unmapped"]
    doc = {"rule": ("face -> gid iff all three vertex gids equal and >= 0; "
                    "mixed/unmapped faces excluded from every submesh and "
                    "counted"),
           "totals": dict(tot), "cubes": out}
    Path(args.out).write_text(json.dumps(doc, indent=1))
    print(f"{len(out)} cubes; faces {tot['faces']}, attributed "
          f"{tot['attributed']} ({tot['attributed']/max(tot['faces'],1):.4f})"
          f", mixed/unmapped {tot['mixed_or_unmapped']}")
    return 0


# ------------------------------------------------------------------ gate

def build_submeshes(placed: Path, endpoints: list[tuple[str, int]]):
    """Triangle arrays for every (cube, gid) endpoint, one cube load each."""
    by_cube = defaultdict(list)
    for cube, g in endpoints:
        by_cube[cube].append(g)
    tris = {}
    counts = {}
    for cube, gids in sorted(by_cube.items()):
        mesh, gid = load_cube(placed, cube)
        fg, _ = face_groups(mesh, gid)
        xyz = mesh.triangle_xyz()                # (T, 3, 3) xyz after swap
        for g in gids:
            sel = fg == g
            tris[(cube, g)] = xyz[sel]
            counts[(cube, g)] = int(sel.sum())
    return tris, counts


def cmd_gate(args) -> int:
    placed, work = Path(args.placed), Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    nodes, edges = read_gg_dump(Path(args.dump))

    endpoints = sorted({(e["cube_a"], e["gid_a"]) for e in edges} |
                       {(e["cube_b"], e["gid_b"]) for e in edges})
    tris, counts = build_submeshes(placed, endpoints)
    order = [ep for ep in endpoints if counts[ep] > 0]
    empty = [ep for ep in endpoints if counts[ep] == 0]
    index = {ep: i for i, ep in enumerate(order)}

    abin = work / "gg_groups.bin"
    wc_atlas.write_atlas_soups([tris[ep] for ep in order], abin)
    csvp = work / "gg_groups.csv"
    cmd = [str(ENGINE), str(abin), str(csvp), str(args.threads),
           CELL, EXCLUDE, DIAGONAL, MAXEDGE]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"engine failed rc={r.returncode}: {r.stderr[-800:]}")

    # cross-surface verdict counts per unordered surface pair
    pair_counts: dict[tuple[int, int], Counter] = defaultdict(Counter)
    with open(csvp) as fh:
        header = fh.readline().strip().split(",")
        i1, i2 = header.index("surf1"), header.index("surf2")
        iv = header.index("verdict")
        for line in fh:
            f = line.rstrip("\n").split(",")
            if f[i1] == f[i2]:
                continue
            key = tuple(sorted((int(f[i1]), int(f[i2]))))
            pair_counts[key][f[iv]] += 1

    out_edges = []
    verdicts = Counter()
    for e in edges:
        pa, pb = (e["cube_a"], e["gid_a"]), (e["cube_b"], e["gid_b"])
        rec = dict(e)
        rec["tris_a"], rec["tris_b"] = counts[pa], counts[pb]
        if counts[pa] == 0 or counts[pb] == 0:
            rec["verdict"] = "not_testable_empty_submesh"
        else:
            key = tuple(sorted((index[pa], index[pb])))
            c = pair_counts.get(key, Counter())
            rec["transverse"] = int(c.get("transverse", 0))
            rec["grazing"] = int(c.get("grazing", 0))
            rec["coplanar"] = int(c.get("coplanar", 0))
            rec["verdict"] = ("transverse_conflict" if rec["transverse"] > 0
                              else "no_transverse_conflict")
        verdicts[rec["verdict"]] += 1
        out_edges.append(rec)

    # every cross-surface contact in the atlas, edge or not, for the record
    n_pairs_with_transverse = sum(
        1 for c in pair_counts.values() if c.get("transverse", 0) > 0)

    doc = {"what": ("prereg Addendum 1 primary analysis: windcheck census "
                    "over the group submeshes of every admitted GGEdge, "
                    "one soup atlas, one engine call, one diagonal pass"),
           "engine_args": {"cell": CELL, "exclude": EXCLUDE,
                           "diagonal": DIAGONAL, "maxedge": MAXEDGE,
                           "threads": args.threads},
           "n_nodes": len(nodes), "n_edges": len(edges),
           "n_endpoint_submeshes": len(endpoints),
           "n_empty_submeshes": len(empty),
           "empty_submeshes": [f"{c}:{g}" for c, g in empty],
           "verdicts": dict(verdicts),
           "n_atlas_pairs_with_transverse_any": n_pairs_with_transverse,
           "edges": out_edges}
    Path(args.out).write_text(json.dumps(doc, indent=1))
    print(f"edges {len(edges)}: {dict(verdicts)}; "
          f"atlas cross-surface pairs with transverse (incl. non-edges): "
          f"{n_pairs_with_transverse}")
    return 0


# --------------------------------------------------------------- assembly

def cmd_assembly(args) -> int:
    """Census the assembled output: every (cube, gid) submesh of every
    placed cube as its own soup surface, one atlas, one engine call.

    Cross-surface rows are the assembly's transverse contacts -- between
    groups of one cube (different wraps of the same cube mesh) and across
    cubes. WITHIN one (cube, gid) submesh nothing is tested (a soup has no
    adjacency); that limitation is disclosed in the record, not hidden.
    Contact sites (sx, sy, sz when present in the CSV) are kept for the
    pre-registered co-location test.
    """
    placed, work = Path(args.placed), Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    cubes = sorted(p.name[:-len("_placed.obj")]
                   for p in placed.glob("*_placed.obj"))
    order = []
    tri_arrays = []
    for cube in cubes:
        mesh, gid = load_cube(placed, cube)
        fg, _ = face_groups(mesh, gid)
        xyz = mesh.triangle_xyz()
        for g in sorted(int(x) for x in np.unique(fg) if x >= 0):
            sel = fg == g
            if not sel.any():
                continue
            order.append((cube, g))
            tri_arrays.append(xyz[sel])
    abin = work / "assembly.bin"
    wc_atlas.write_atlas_soups(tri_arrays, abin)
    csvp = work / "assembly.csv"
    cmd = [str(ENGINE), str(abin), str(csvp), str(args.threads),
           CELL, EXCLUDE, DIAGONAL, MAXEDGE]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"engine failed rc={r.returncode}: {r.stderr[-800:]}")

    rows = []
    with open(csvp) as fh:
        header = fh.readline().strip().split(",")
        i1, i2 = header.index("surf1"), header.index("surf2")
        iv = header.index("verdict")
        for line in fh:
            f = line.rstrip("\n").split(",")
            if f[i1] == f[i2] or f[iv] != "transverse":
                continue
            a, b = int(f[i1]), int(f[i2])
            rows.append({"a": f"{order[a][0]}:{order[a][1]}",
                         "b": f"{order[b][0]}:{order[b][1]}",
                         "same_cube": order[a][0] == order[b][0],
                         "fields": dict(zip(header, f))})
    pair_c = Counter((r["a"], r["b"]) for r in rows)
    doc = {"what": ("assembly census: every (cube,gid) submesh a soup "
                    "surface; cross-surface transverse only. Within-"
                    "(cube,gid) contacts are UNTESTED by construction "
                    "(soups carry no adjacency) -- disclosed limitation"),
           "engine_args": {"cell": CELL, "exclude": EXCLUDE,
                           "diagonal": DIAGONAL, "maxedge": MAXEDGE},
           "n_surfaces": len(order),
           "n_transverse_rows": len(rows),
           "n_pairs_with_transverse": len(pair_c),
           "n_same_cube_rows": sum(1 for r in rows if r["same_cube"]),
           "n_cross_cube_rows": sum(1 for r in rows if not r["same_cube"]),
           "rows": rows}
    Path(args.out).write_text(json.dumps(doc, indent=1))
    print(f"assembly: {len(order)} group submeshes; transverse rows "
          f"{len(rows)} over {len(pair_c)} pairs "
          f"(same-cube {doc['n_same_cube_rows']}, cross-cube "
          f"{doc['n_cross_cube_rows']}); csv {csvp}")
    return 0


# ----------------------------------------------------------------- unroll

def read_faces_tsv(path: Path):
    """faces.tsv -> (W, H, cube_id[], gid[]) with face index as position."""
    W = H = nf = None
    cubes, gids = [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith("# facemap"):
            parts = dict(p.split("=") for p in line.split()[3:])
            W, H, nf = int(parts["W"]), int(parts["H"]), int(parts["nf"])
            continue
        if line.startswith("#") or not line:
            continue
        f, cube, gid = line.split("\t")
        assert int(f) == len(cubes)
        cubes.append(cube)
        gids.append(int(gid))
    assert nf == len(cubes), (nf, len(cubes))
    return W, H, cubes, np.asarray(gids, dtype=np.int32)


def cmd_unroll(args) -> int:
    """Census one unroll export and classify events per the FROZEN prereg
    rule (PREREG-UNROLL-GATE.md): cell source set = {(cube,gid) of facemap
    at its four corner pixels, omitting empty}; an event is cross_source
    iff both sets are non-empty and disjoint, same_source iff they
    intersect, unattributed iff either is empty."""
    import tifffile
    sys.path.insert(0, str(REPO_ROOT / "bench"))
    from crossing_census import census_one
    from excise_segment import CENSUS, parse_census_csv

    export, work = Path(args.export), Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    tag = args.tag

    W, H, cubes, face_gid = read_faces_tsv(
        Path(args.facemap_prefix + "_faces.tsv"))
    fm = np.fromfile(args.facemap_prefix + "_facemap.i32", dtype="<i4")
    assert len(fm) == W * H, (len(fm), W, H)
    fm = fm.reshape(H, W)

    # instrumentation consistency: facemap empty <=> provenance state 0
    prov_path = export / "provenance.tif"
    prov_ok = None
    if prov_path.exists():
        prov = np.asarray(tifffile.imread(prov_path))
        assert prov.shape == (H, W), (prov.shape, H, W)
        prov_ok = bool(((fm == -1) == (prov == 0)).all())
        if not prov_ok:
            n_bad = int(((fm == -1) != (prov == 0)).sum())
            raise SystemExit(f"facemap/provenance DISAGREE on {n_bad} px -- "
                             "instrumentation bug, stopping per prereg")

    row = census_one(export, tag, CENSUS["exclude"], CENSUS["cell"],
                     CENSUS["threads"], CENSUS["maxedge"], work)
    if row is None:
        raise SystemExit(f"{export}: below census validity floor")

    def cell_sources(v, u):
        out = set()
        for dv in (0, 1):
            for du in (0, 1):
                vv, uu = v + dv, u + du
                if 0 <= vv < H and 0 <= uu < W:
                    f = int(fm[vv, uu])
                    if f >= 0:
                        out.add((cubes[f], int(face_gid[f])))
        return out

    events = []
    counts = Counter()
    for diag in (0, 1):
        rec = parse_census_csv(Path(row[f"csv_d{diag}"]), diag,
                               row["grid"][0], row["grid"][1])
        for r in rec["rows"]:
            if r["verdict"] != "transverse":
                continue
            (v1, u1), (v2, u2) = r["q1"], r["q2"]
            sa = cell_sources(v1, u1)
            sb = cell_sources(v2, u2)
            if not sa or not sb:
                cls = "unattributed"
            elif sa.isdisjoint(sb):
                cls = "cross_source"
            else:
                cls = "same_source"
            counts[f"d{diag}_{cls}"] += 1
            counts[cls] += 1
            events.append({"diag": diag, "class": cls,
                           "v1": v1, "u1": u1, "v2": v2, "u2": u2,
                           "penetration": r["penetration"],
                           "angle_deg": r["angle_deg"],
                           "sources_a": sorted(f"{c0}:{g0}"
                                               for c0, g0 in sa),
                           "sources_b": sorted(f"{c0}:{g0}"
                                               for c0, g0 in sb)})

    doc = {"export": str(export), "tag": tag,
           "prereg": "notes/PREREG-UNROLL-GATE.md",
           "census": {k: row[k] for k in
                      ("grid", "valid_cells", "params", "d0", "d1")},
           "facemap_provenance_consistent": prov_ok,
           "n_transverse_rows": len(events),
           "classes": dict(counts),
           "events": events}
    Path(args.out).write_text(json.dumps(doc, indent=1))
    print(f"{tag}: grid {row['grid']} valid {row['valid_cells']}; "
          f"transverse d0={row['d0']['transverse']} "
          f"d1={row['d1']['transverse']}; classes {dict(counts)}; "
          f"provenance_consistent={prov_ok}")
    return 0


def cmd_plant(args) -> int:
    """Positive control: copy the export, shift a band of columns by
    +0.5 vx in x, census the copy. Must produce transverse events."""
    import shutil
    import tifffile
    sys.path.insert(0, str(REPO_ROOT / "bench"))
    from crossing_census import census_one
    from excise_segment import CENSUS

    export, work = Path(args.export), Path(args.workdir)
    dst = work / "planted.tifxyz"
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(export, dst)
    x = np.asarray(tifffile.imread(dst / "x.tif")).copy()
    H, W = x.shape
    lo, hi = W // 3, W // 3 + max(W // 20, 8)
    valid = x > 0
    x[:, lo:hi][valid[:, lo:hi]] += 0.5
    tifffile.imwrite(dst / "x.tif", x)
    row = census_one(dst, "planted", CENSUS["exclude"], CENSUS["cell"],
                     CENSUS["threads"], CENSUS["maxedge"], work)
    t = (row["d0"]["transverse"], row["d1"]["transverse"]) if row else None
    print(f"planted control: transverse d0/d1 = {t} -> "
          f"{'DETECTS' if t and max(t) > 0 else 'FAILS'}")
    return 0 if t and max(t) > 0 else 1


# -------------------------------------------------------------- spotcheck

def cmd_spotcheck(args) -> int:
    """Independent re-parse of sampled triangles from the raw OBJ text."""
    placed = Path(args.placed)
    rng = np.random.default_rng(20260805)
    cubes = sorted(placed.glob("*_placed.obj"))
    picks = rng.choice(len(cubes), size=min(args.cubes, len(cubes)),
                       replace=False)
    checked = 0
    for ci in picks:
        objp = cubes[int(ci)]
        cube = objp.name[:-len("_placed.obj")]
        mesh, gid = load_cube(placed, cube)
        xyz = mesh.triangle_xyz()
        # independent parse: raw text, no shared code path
        vlines = [l for l in objp.read_text().splitlines()
                  if l.startswith("v ")]
        flines = [l for l in objp.read_text().splitlines()
                  if l.startswith("f ")]
        assert len(vlines) == len(mesh.verts)
        assert len(flines) == len(mesh.tris)
        for ti in rng.choice(len(flines), size=min(args.per_cube,
                                                   len(flines)),
                             replace=False):
            toks = flines[int(ti)].split()[1:4]
            vidx = [int(t.split("/")[0]) - 1 for t in toks]
            got = xyz[int(ti)]
            for corner, vi in enumerate(vidx):
                z, y, x = (np.float32(v) for v in
                           vlines[vi].split()[1:4])
                want = np.array([x, y, z], dtype=np.float32)
                if not np.array_equal(got[corner], want):
                    raise SystemExit(
                        f"MISMATCH {cube} face {ti} corner {corner}: "
                        f"atlas {got[corner]} vs obj {want}")
            checked += 1
    print(f"spotcheck OK: {checked} triangles across {len(picks)} cubes, "
          "atlas geometry == raw OBJ text (zyx->xyz)")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("attribute")
    a.add_argument("--placed", required=True)
    a.add_argument("--out", required=True)
    g = sub.add_parser("gate")
    g.add_argument("--placed", required=True)
    g.add_argument("--dump", required=True)
    g.add_argument("--workdir", required=True)
    g.add_argument("--out", required=True)
    g.add_argument("--threads", type=int, default=0)
    s = sub.add_parser("spotcheck")
    s.add_argument("--placed", required=True)
    s.add_argument("--cubes", type=int, default=12)
    s.add_argument("--per-cube", type=int, default=40)
    m = sub.add_parser("assembly")
    m.add_argument("--placed", required=True)
    m.add_argument("--workdir", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--threads", type=int, default=0)
    u = sub.add_parser("unroll")
    u.add_argument("--export", required=True)
    u.add_argument("--facemap-prefix", required=True)
    u.add_argument("--workdir", required=True)
    u.add_argument("--out", required=True)
    u.add_argument("--tag", required=True)
    p = sub.add_parser("plant")
    p.add_argument("--export", required=True)
    p.add_argument("--workdir", required=True)
    args = ap.parse_args(argv)
    return {"attribute": cmd_attribute, "gate": cmd_gate,
            "spotcheck": cmd_spotcheck, "assembly": cmd_assembly,
            "unroll": cmd_unroll, "plant": cmd_plant}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
