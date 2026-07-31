"""The two things a user actually does to one segment: check it, or transform it.

`check` is the default posture of the tool and is strictly report-only. It runs
the both-diagonal, non-adjacent transverse-contact census on a `.tifxyz`
segment, prints a verdict, writes a machine-readable certificate, and writes the
crossing sites as a point collection the viewer can load. It never writes a mesh
and never touches the input directory.

`transform` is the explicit path. It applies the frozen scheduling policy
implemented in `windcheck.excise` -- bounded displacement first where a
displacement-repaired mesh exists for the segment, then one certified excision
of every residual transverse contact -- and emits one aggregate `.tifxyz` plus
its certificate. Nothing here runs unless the user asked for it by name.

Both paths share one authority for the word "clean": the C++ census in
`engines/selfcross`, run under BOTH quad triangulations (diagonal 0 and
diagonal 1), with Chebyshev adjacency exclusion so that quads sharing a corner
are never reported against each other. A mesh is clean when both diagonals
report zero transverse contacts.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import time
from collections import Counter
from math import fsum
from pathlib import Path

import numpy as np
import tifffile

from . import atlas, tifxyz
from .certificate import write_collection
from .excise import (FROZEN_POLICY_VERSION, GEOMETRY_STATUS_CLEAN,
                     HYBRID_INVALIDATION, MISSING, REDUCTION_RULE,
                     RETAINED_BIT_IDENTITY, SELECTION_STATUS_RULE,
                     STALENESS_WARNING, frozen_policy_hash, quad_corners,
                     quad_triangle_corners, select_global_frozen)
from .intrinsic import retained_quads
from .manifest import (MANIFEST_RULE, BaseKind, mesh_manifest,
                       verify_base_kind)
from .provenance import release_provenance

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / "engines" / "selfcross"
ENGINE_SRC = ROOT / "engines" / "selfcross.cpp"

AXES = ("x", "y", "z")

#: Census parameters. These are not tuning knobs: every certificate in this
#: project was produced under them, and changing one makes two results
#: incomparable, so they travel into the certificate verbatim.
CENSUS = {"cell": 40.0, "exclude": 1, "maxedge": 60.0, "threads": 0,
          "touch_tol": 1e-3}

SCHEMA_V2_HEADER = "v1,u1,v2,u2,verdict,penetration,angle_deg,tri1,tri2"
VERDICTS = ("transverse", "coplanar", "grazing")

#: Where a bounded-displacement repair of a segment is looked for, relative to
#: the working directory, when the user does not name one with `--displaced`.
DISPLACED_ROOT = Path("out/repaired/multi/meshes")

CLEAN_DEFINITION = (
    "CLEAN means: the mesh, censused by engines/selfcross under BOTH quad "
    "triangulations with Chebyshev adjacency exclusion, reports zero "
    "transverse contacts. Coplanar and grazing contacts are recorded but are "
    "not crossings and do not affect the verdict.")

CHECK_CAVEATS = [
    "Deterministic floating-point validator, not exact predicates.",
    "Reports that the surface meets itself and where; does not identify which "
    "branch is wrong, nor establish a cause.",
    "Adjacent quads (sharing a corner) are excluded by construction, so a "
    "reported contact is always non-local in the grid.",
    "Event counts depend on the grid-connectivity clustering rule; the "
    "per-diagonal contact counts do not.",
    "The mesh fingerprint is a partial one: file names, sizes and a bounded "
    "head of each file, not full content.",
    "This command is report-only. It does not modify the input and does not "
    "write a mesh.",
]

TRANSFORM_CAVEATS = [
    "The clean claim is made by the C++ census of the EMITTED arrays reloaded "
    "from disk, under both diagonals -- never by the solver that chose the cut.",
    "Retained coordinates are bit-identical to the input; the only cells whose "
    "coordinates differ are the excised ones.",
    "Excised cells are marked both ways (mask sidecar and the -1 coordinate "
    "marker), so a consumer honouring either convention sees the same surface.",
    STALENESS_WARNING,
]


class CensusRefusal(RuntimeError):
    """The mesh is below the census validity threshold; no verdict is defined."""


# --------------------------------------------------------------- discovery

def find_mesh(target: Path, volume: str = "") -> Path | None:
    """Accept a `.tifxyz` directly, or a segment directory to search inside."""
    target = Path(target)
    if (target / "x.tif").exists() and (target / "z.tif").exists():
        return target
    for pat in (f"mesh/*{volume}*.tifxyz", f"*{volume}*.tifxyz"):
        hits = sorted(target.glob(pat))
        if hits:
            return hits[0]
    return None


def segment_name(mesh: Path) -> str:
    """The segment a mesh belongs to: its grandparent when it sits in `mesh/`."""
    mesh = Path(mesh)
    if mesh.parent.name == "mesh" and mesh.parent.parent.name:
        return mesh.parent.parent.name
    return mesh.name.removesuffix(".tifxyz")


def find_displaced(name: str, root: Path = DISPLACED_ROOT) -> Path | None:
    """A bounded-displacement repair of this segment, if one is on disk."""
    cand = Path(root) / f"{name}_repaired.tifxyz"
    if cand.is_dir() and (cand / "x.tif").exists():
        return cand
    return None


def sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mesh_hashes(mesh: Path) -> dict:
    """Mesh identity: the canonical manifest of every file a reader consumes.

    Not the coordinate planes alone. `mask.tif` decides which triangles
    exist -- it is the carrier of an excision -- and `meta.json` is
    provenance-relevant, so both are inside the digest. The returned dict
    carries `files` (one explicit row per declared file, present or absent),
    `digest` (sha256 over the frozen serialisation of those rows) and, as a
    view derived from the same rows, the `x`/`y`/`z`/`mask` spellings older
    records used.
    """
    return mesh_manifest(mesh)


def _voxel_um(mesh: Path) -> float | None:
    import re
    m = re.search(r"-(\d+\.?\d*)um\.tifxyz$", Path(mesh).name)
    return float(m.group(1)) if m else None


# ------------------------------------------------------------ the census

def run_engine(mesh: Path, tag: str, work: Path, diag: int,
               params: dict | None = None) -> tuple[Path, dict]:
    """One engine pass over one triangulation. Returns (csv path, counts)."""
    p = dict(CENSUS, **(params or {}))
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    abin = work / f"{tag[:40]}_atlas.bin"
    atlas.write_atlas([_AtlasEntry(mesh)], abin)
    csv = work / f"{tag[:40]}_d{diag}.csv"
    r = subprocess.run(
        [str(ENGINE), str(abin), str(csv), str(int(p["threads"])),
         str(p["cell"]), str(int(p["exclude"])), str(diag), str(p["maxedge"])],
        capture_output=True, text=True, check=True)
    counts = json.loads(r.stdout.strip().splitlines()[-1])
    return csv, {k: counts[k] for k in
                 ("triangles", "quads_dropped", "pairs_tested",
                  "transverse", "coplanar", "grazing")}


class _AtlasEntry:
    def __init__(self, path: Path):
        self.path, self.winding = Path(path), None


def parse_census_csv(csv: Path, diag: int, nv: int, nu: int) -> dict:
    """Parse one schema-v2 census CSV, keeping the engine's own identities.

    Rows are kept verbatim as `(diag, v1, u1, tri1, v2, u2, tri2, verdict)`, so
    a contact multiset from one mesh can be compared against another's. Rows
    outside the grid's quad range are counted rather than dropped: the engine
    never emits them, so a nonzero count is evidence of a grid mismatch.
    """
    text = Path(csv).read_text().splitlines()
    if not text:
        raise ValueError(f"empty census CSV: {csv}")
    header = text[0].strip()
    if header != SCHEMA_V2_HEADER:
        raise ValueError(f"unexpected census schema in {csv}: {header!r}")
    rows: list[dict] = []
    multiset: Counter = Counter()
    out_of_range = 0
    for line in text[1:]:
        if not line.strip():
            continue
        p = line.split(",")
        verdict = p[4]
        if verdict not in VERDICTS:
            raise ValueError(f"unknown verdict {verdict!r} in {csv}")
        v1, u1, v2, u2 = int(p[0]), int(p[1]), int(p[2]), int(p[3])
        t1, t2 = int(p[7]), int(p[8])
        if not (0 <= v1 < nv - 1 and 0 <= u1 < nu - 1
                and 0 <= v2 < nv - 1 and 0 <= u2 < nu - 1):
            out_of_range += 1
            continue
        key = (diag, v1, u1, t1, v2, u2, t2, verdict)
        multiset[key] += 1
        rows.append({"key": key, "diag": diag, "q1": (v1, u1), "q2": (v2, u2),
                     "t1": t1, "t2": t2, "verdict": verdict,
                     "penetration": float(p[5]), "angle_deg": float(p[6])})
    return {"rows": rows, "multiset": multiset, "out_of_range": out_of_range,
            "header": header, "n_lines": len(text) - 1,
            "path": str(csv), "sha256": sha(csv)}


def census(mesh: Path, tag: str, work: Path, nv: int, nu: int,
           params: dict | None = None) -> dict:
    """The authoritative both-diagonal census of a mesh on disk."""
    surf_valid = int(np.asarray(tifxyz.read(mesh).valid, bool).sum())
    if surf_valid < 5000:
        raise CensusRefusal(
            f"{mesh}: {surf_valid} valid cells is below the census validity "
            "threshold of 5000; no verdict is defined on it")
    per_diag, engine = {}, {}
    for d in (0, 1):
        csv, counts = run_engine(mesh, tag, work, d, params)
        engine[f"d{d}"] = counts
        per_diag[d] = parse_census_csv(csv, d, nv, nu)
    multiset: Counter = Counter()
    for d in (0, 1):
        multiset.update(per_diag[d]["multiset"])
    rows = [r for d in (0, 1) for r in per_diag[d]["rows"]]
    return {"tag": tag, "engine": engine, "diag": per_diag,
            "multiset": multiset, "rows": rows,
            "transverse_total": sum(engine[f"d{d}"]["transverse"]
                                    for d in (0, 1)),
            "clean": all(engine[f"d{d}"]["transverse"] == 0 for d in (0, 1))}


# ------------------------------------------------------------- geometry

def quad_area_grids(P64: np.ndarray) -> dict:
    """Per-quad areas for the whole grid, under each diagonal and canonical.

    The area of a quad under diagonal `d` is the sum of its two triangle areas
    under that tessellation; the canonical value is the mean of the two, so it
    does not privilege either triangulation.
    """
    p00, p10 = P64[:-1, :-1], P64[1:, :-1]
    p01, p11 = P64[:-1, 1:], P64[1:, 1:]

    def tri(a, b, c):
        return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=-1)

    a0 = tri(p00, p01, p11) + tri(p00, p11, p10)
    a1 = tri(p00, p01, p10) + tri(p01, p11, p10)
    return {0: a0, 1: a1, "canonical": 0.5 * (a0 + a1)}


def coverage_from_rows(rows: list[dict], Q: np.ndarray) -> list[dict]:
    """Map transverse census rows to the coverage sets the selector expects.

    Only TRANSVERSE rows become constraints: the acceptance criterion is zero
    transverse contacts, and coplanar or grazing contacts are recorded in the
    multiset but are not crossings. Coverage is the eight-corner set
    `corners(quad1) | corners(quad2)`; invalidating any one of those eight
    vertices destroys at least one of the two crossing triangles.
    """
    out = []
    for r in rows:
        if r["verdict"] != "transverse":
            continue
        (v1, u1), (v2, u2) = r["q1"], r["q2"]
        for (v, u) in (r["q1"], r["q2"]):
            if not Q[v, u]:
                raise ValueError(
                    f"censused quad ({v},{u}) is not a retained quad: "
                    "engine and library retention semantics diverge")
        c1 = quad_triangle_corners(v1, u1, r["diag"])[r["t1"]]
        c2 = quad_triangle_corners(v2, u2, r["diag"])[r["t2"]]
        cov1, cov2 = set(quad_corners(v1, u1)), set(quad_corners(v2, u2))
        if cov1 & cov2:
            raise ValueError(
                f"crossing quads ({v1},{u1}) and ({v2},{u2}) share a corner: "
                "adjacency exclusion violated")
        out.append({"key": r["key"], "q1": r["q1"], "q2": r["q2"],
                    "t1": r["t1"], "t2": r["t2"], "diag": r["diag"],
                    "corners1": c1, "corners2": c2,
                    "participants": sorted(set(c1) | set(c2)),
                    "coverage": sorted(cov1 | cov2)})
    return out


def crossing_sites(rows: list[dict], P: np.ndarray) -> list[dict]:
    """One representative volume location per grid-connected crossing region.

    A crossing is a curve over many quads; listing every contact would drown
    the viewer. Participating quads are grouped by grid adjacency, each pair of
    groups is one event, and the event's site is the midpoint of the deepest
    contact in it -- which is inside the overlap, so clicking it lands on the
    problem rather than near it.
    """
    trans = [r for r in rows if r["verdict"] == "transverse"]
    if not trans:
        return []
    quads = set()
    for r in trans:
        quads.add(r["q1"])
        quads.add(r["q2"])
    parent = {q: q for q in quads}

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

    best: dict = {}
    for r in trans:
        k = tuple(sorted((find(r["q1"]), find(r["q2"]))))
        if k not in best or r["penetration"] > best[k]["penetration"]:
            best[k] = r
    sites = []
    for k in sorted(best):
        r = best[k]
        (v1, u1), (v2, u2) = r["q1"], r["q2"]
        mid = 0.5 * (np.asarray(P[v1, u1], np.float64)
                     + np.asarray(P[v2, u2], np.float64))
        sites.append({"p": [float(x) for x in mid],
                      "quad_a": [int(v1), int(u1)],
                      "quad_b": [int(v2), int(u2)],
                      "diagonal": int(r["diag"]),
                      "penetration_vx": float(r["penetration"]),
                      "angle_deg": float(r["angle_deg"])})
    return sites


# ------------------------------------------------------------- emission

def emit_tifxyz(src: Path, dst: Path, valid_out: np.ndarray,
                excised: np.ndarray, *, with_mask: bool = True) -> dict:
    """Emit one aggregate `.tifxyz` under hybrid invalidation.

    tifxyz marks a missing cell either by setting x=y=z=-1 or through a mask
    sidecar. Writing only one of them is an interoperability defect: a consumer
    implementing the other convention reads the excised cells as present and
    reconstructs exactly the crossing that was removed. So both are written,
    and every retained coordinate is copied through bit-identically.
    """
    src, dst = Path(src), Path(dst)
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    for stale in ("mask.tif", "mask.png"):
        if (src / stale).exists():
            raise ValueError(
                f"input carries {stale}: emitting here would overwrite input "
                "validity semantics; refusing")
    for ax in AXES:
        plane = np.asarray(tifffile.imread(src / f"{ax}.tif"))
        if plane.shape != excised.shape:
            raise ValueError(f"{ax}.tif shape {plane.shape} != grid "
                             f"{excised.shape}")
        out = plane.copy()
        out[excised] = MISSING
        tifffile.imwrite(dst / f"{ax}.tif", out)
    if (src / "meta.json").exists():
        shutil.copy(src / "meta.json", dst / "meta.json")
    if with_mask:
        tifffile.imwrite(dst / "mask.tif", valid_out.astype(np.uint8))
    return {"invalidation_carrier": HYBRID_INVALIDATION,
            "mask_written": with_mask,
            "coordinate_guarantee": RETAINED_BIT_IDENTITY,
            "files": {f.name: sha(f) for f in sorted(dst.iterdir())}}


# ---------------------------------------------------------------- check

def check_segment(target: Path, out: Path, *, volume: str = "",
                  threads: int = 0, cell: float = CENSUS["cell"],
                  maxedge: float = CENSUS["maxedge"],
                  name: str | None = None) -> dict:
    """Report-only. Census, verdict, certificate, point collection. No mesh."""
    t0 = time.time()
    mesh = find_mesh(target, volume)
    if mesh is None:
        raise FileNotFoundError(f"no .tifxyz surface found under {target}")
    if not ENGINE.exists():
        raise FileNotFoundError(
            f"census engine not built: {ENGINE}\n"
            "  clang++ -O3 -std=c++17 -pthread -o engines/selfcross "
            "engines/selfcross.cpp")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    seg = name or segment_name(mesh)
    params = dict(CENSUS, threads=threads, cell=cell, maxedge=maxedge)

    surf = tifxyz.read(mesh)
    P = np.asarray(surf.points, np.float64)
    V = np.asarray(surf.valid, bool)
    nv, nu = surf.shape
    Q = retained_quads(P, V, maxedge)

    work = out / f"work_{hashlib.sha256(seg.encode()).hexdigest()[:12]}"
    cen = census(mesh, f"{seg[:40]}_check", work, nv, nu, params)
    sites = crossing_sites(cen["rows"], P)

    points_path = out / f"{seg[:60]}_points.json"
    write_collection(points_path, f"windcheck-crossings-{seg[:32]}",
                     [s["p"] for s in sites],
                     metadata={"tool": "windcheck check",
                               "n_crossing_events": len(sites)})

    wall = time.time() - t0
    cert = {
        "tool": "windcheck check",
        "record_kind": "census certificate",
        "schema": "windcheck_check/v1",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "segment": seg,
        "report_only": True,
        "report_only_note": ("This command read the input and wrote nothing "
                             "into it. No mesh was produced. Use `windcheck "
                             "transform` to change geometry."),
        "provenance": release_provenance(),
        "mesh_identity_rule": MANIFEST_RULE,
        "mesh": {"path": str(mesh), "hashes": mesh_hashes(mesh),
                 "grid_shape": [int(nv), int(nu)],
                 "n_valid_vertices": int(V.sum()),
                 "n_retained_quads": int(Q.sum()),
                 "voxel_um": _voxel_um(mesh)},
        "census": {
            "engine": str(ENGINE),
            "engine_sha256": sha(ENGINE) if ENGINE.exists() else None,
            "engine_source_sha256": (sha(ENGINE_SRC) if ENGINE_SRC.exists()
                                     else None),
            "parameters": dict(params, diagonals=[0, 1]),
            "csv": {f"d{d}": {"path": cen["diag"][d]["path"],
                              "sha256": cen["diag"][d]["sha256"],
                              "n_rows": cen["diag"][d]["n_lines"],
                              "out_of_range_rows": cen["diag"][d]["out_of_range"]}
                    for d in (0, 1)}},
        "measurements": {
            "diagonal_0": cen["engine"]["d0"],
            "diagonal_1": cen["engine"]["d1"],
            "transverse_d0": cen["engine"]["d0"]["transverse"],
            "transverse_d1": cen["engine"]["d1"]["transverse"],
            "transverse_both_diagonals": cen["transverse_total"],
            "transverse_min_of_diagonals": min(cen["engine"]["d0"]["transverse"],
                                               cen["engine"]["d1"]["transverse"]),
            "crossing_events": len(sites),
            "max_penetration_vx": (max((s["penetration_vx"] for s in sites),
                                       default=0.0)),
        },
        "clean": bool(cen["clean"]),
        "verdict": ("clean: no transverse self-intersection under either "
                    "triangulation" if cen["clean"] else
                    f"NOT clean: {cen['transverse_total']} transverse contacts "
                    f"over {len(sites)} crossing events"),
        "clean_definition": CLEAN_DEFINITION,
        "crossing_sites": sites[:2000],
        "crossing_sites_truncated_to": 2000,
        "outputs": {"point_collection": str(points_path),
                    "point_collection_format": (
                        "volume-cartographer PointCollections JSON; one "
                        "collection, one point per crossing event, `p` is "
                        "[x, y, z] in volume voxels")},
        "wall_seconds": round(wall, 2),
        "caveats": CHECK_CAVEATS,
    }
    cert_path = out / f"{seg[:60]}_check_certificate.json"
    cert_path.write_text(json.dumps(cert, indent=1, default=str) + "\n")

    return {"segment": seg, "mesh": mesh, "grid": [int(nv), int(nu)],
            "triangles": sum(cen["engine"][f"d{d}"]["triangles"]
                             for d in (0, 1)) // 2,
            "clean": bool(cen["clean"]),
            "transverse_d0": cen["engine"]["d0"]["transverse"],
            "transverse_d1": cen["engine"]["d1"]["transverse"],
            "coplanar_d0": cen["engine"]["d0"]["coplanar"],
            "coplanar_d1": cen["engine"]["d1"]["coplanar"],
            "events": len(sites), "verdict": cert["verdict"],
            "certificate": cert_path, "points": points_path,
            "n_points": len(sites), "wall_seconds": wall, "out": out}


# ------------------------------------------------------------ transform

def transform_segment(target: Path, out: Path, *, volume: str = "",
                      threads: int = 0, cell: float = CENSUS["cell"],
                      maxedge: float = CENSUS["maxedge"],
                      displaced: Path | None = None,
                      use_displaced: bool = True,
                      name: str | None = None) -> dict:
    """The explicit transformation: frozen policy, one aggregate mesh out.

    Bounded displacement where applicable -- a displacement-repaired mesh for
    this segment becomes the base that is cut -- then ONE certified excision of
    every residual transverse contact, chosen under the frozen policy in
    `windcheck.excise`. The clean claim comes from recensusing the emitted
    arrays reloaded from disk, never from the selector.
    """
    t0 = time.time()
    original = find_mesh(target, volume)
    if original is None:
        raise FileNotFoundError(f"no .tifxyz surface found under {target}")
    if not ENGINE.exists():
        raise FileNotFoundError(f"census engine not built: {ENGINE}")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    seg = name or segment_name(original)
    params = dict(CENSUS, threads=threads, cell=cell, maxedge=maxedge)

    base = Path(displaced) if displaced else (
        find_displaced(seg) if use_displaced else None)
    if base is not None and not (base / "x.tif").exists():
        raise FileNotFoundError(f"not a .tifxyz surface: {base}")
    mesh = base or original
    # Declared from the selector, then confirmed against content. The
    # decision is never a path comparison: paths move between downloaded
    # archives and fresh workdirs, manifests do not.
    declared_kind = (BaseKind.DISPLACEMENT_REPAIRED if base is not None
                     else BaseKind.ORIGINAL)
    input_manifest = mesh_hashes(mesh)
    original_manifest = mesh_hashes(original)
    base_kind_check = verify_base_kind(declared_kind, input_manifest,
                                       original_manifest, label=seg)
    # Content is the authority. A "repair" that emitted a byte-identical
    # mesh IS the original, whatever the selector called it.
    base_kind = base_kind_check["verified"] or declared_kind.value

    surf = tifxyz.read(mesh)
    P64 = np.asarray(surf.points, np.float64)
    V_in = np.asarray(surf.valid, bool)
    nv, nu = surf.shape
    Q_in = retained_quads(P64, V_in, maxedge)
    areas = quad_area_grids(P64)
    a_in = fsum(areas["canonical"][Q_in].tolist())

    work = out / f"work_{hashlib.sha256(seg.encode()).hexdigest()[:12]}"
    mesh_out = out / f"{seg[:60]}_transformed.tifxyz"
    cert_path = out / f"{seg[:60]}_transform_certificate.json"

    R: dict = {
        "tool": "windcheck transform",
        "record_kind": "transformation certificate",
        "schema": "windcheck_transform/v1",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "segment": seg,
        "policy_version": FROZEN_POLICY_VERSION,
        "policy_hash": frozen_policy_hash(),
        "displacement": {
            "applied": base is not None,
            "base_mesh": str(mesh), "base_kind": base_kind,
            "base_kind_declared": declared_kind.value,
            "base_kind_verification": base_kind_check,
            "base_kind_rule": (
                "base_kind is decided by semantic manifest equality between "
                "the input mesh and the original published mesh, never by "
                "path"),
            "rule": ("bounded displacement first where a displacement-"
                     "repaired mesh exists for this segment, then ONE global "
                     "residual excision")},
        "provenance": release_provenance(
            policy_version=FROZEN_POLICY_VERSION,
            policy_hash=frozen_policy_hash()),
        "mesh_identity_rule": MANIFEST_RULE,
        "original_mesh": str(original),
        "original_mesh_hashes": original_manifest,
        "input_mesh": str(mesh),
        "input_mesh_hashes": input_manifest,
        "input_grid_shape": [int(nv), int(nu)],
        "input_n_valid_vertices": int(V_in.sum()),
        "input_n_retained_quads": int(Q_in.sum()),
        "input_area_canonical": a_in,
        "census_parameters": dict(params, diagonals=[0, 1]),
        "census_engine": str(ENGINE),
        "census_engine_sha256": sha(ENGINE),
        "constraint_reduction_rule": REDUCTION_RULE,
        "selection_status_rule": SELECTION_STATUS_RULE,
        "invalidation": HYBRID_INVALIDATION,
        "clean_definition": CLEAN_DEFINITION,
        "caveats": TRANSFORM_CAVEATS,
        "status": "started",
        "output_mesh": None,
        "claimed_clean": False,
    }

    def finish(**kw):
        R.update(kw)
        R["wall_seconds"] = round(time.time() - t0, 2)
        cert_path.write_text(json.dumps(R, indent=1, default=str) + "\n")
        return R

    if int(V_in.sum()) == 0 or int(Q_in.sum()) == 0:
        finish(status="triangle_empty_invalid",
               note=("the input carries no triangles, so no census, no cut "
                     "and no cleanliness claim is defined on it"))
        return {"segment": seg, "status": R["status"], "certificate": cert_path,
                "mesh_out": None, "wall_seconds": R["wall_seconds"], "out": out}

    before = census(mesh, f"{seg[:40]}_before", work, nv, nu, params)
    R["census_before"] = {f"d{d}": before["engine"][f"d{d}"] for d in (0, 1)}

    if before["clean"]:
        # Nothing to cut. The output is still emitted, so `transform` always
        # produces one artifact and a downstream consumer never has to branch
        # on whether a cut happened.
        emitted = emit_tifxyz(mesh, mesh_out, V_in,
                              np.zeros((nv, nu), bool))
        finish(status="already_clean", output_mesh=str(mesh_out),
               output_mesh_hashes=mesh_hashes(mesh_out),
               emission=emitted, n_removed_quads=0,
               n_invalidated_vertices=0,
               retained_area_fraction=1.0, excised_area_fraction=0.0,
               claimed_clean=True, geometry_status=GEOMETRY_STATUS_CLEAN,
               note=("the input was already transverse-clean; the output is "
                     "the input's coordinates, unchanged"))
        return {"segment": seg, "status": R["status"], "certificate": cert_path,
                "mesh_out": mesh_out, "removed_quads": 0,
                "retained_fraction": 1.0, "before": R["census_before"],
                "after": R["census_before"], "clean": True,
                "wall_seconds": R["wall_seconds"], "out": out}

    cons = coverage_from_rows(before["rows"], Q_in)
    sel = select_global_frozen(cons, P64, Q_in, (),
                               area_grid=areas["canonical"])
    R["selection"] = {k: sel[k] for k in
                      ("status", "selection_status", "policy_version",
                       "policy_hash", "method_mix", "achieved_area",
                       "combined_lower_bound", "combined_lower_bound_complete",
                       "ratio_achieved_over_bound",
                       "minimum_area_claim_admissible", "timings",
                       "scipy_version") if k in sel}
    R["selection"]["n_constraints"] = len(cons)
    R["selection"]["reduction"] = sel.get("reduction")
    if sel["status"] != "ok":
        finish(status=f"selection_{sel['status']}",
               note=("no feasible mask was found: NO output mesh and NO clean "
                     "claim. This is a recorded result, not a crash."))
        return {"segment": seg, "status": R["status"], "certificate": cert_path,
                "mesh_out": None, "wall_seconds": R["wall_seconds"], "out": out}

    X = np.zeros((nv, nu), bool)
    for v, u in sel["chosen"]:
        X[v, u] = True
    removed = (X[:-1, :-1] | X[1:, :-1] | X[:-1, 1:] | X[1:, 1:]) & Q_in
    kept = Q_in & ~removed
    V_out = V_in & ~X
    excised_cells = V_in & X
    a_kept = fsum(areas["canonical"][kept].tolist())
    a_removed = fsum(areas["canonical"][removed].tolist())

    emitted = emit_tifxyz(mesh, mesh_out, V_out, excised_cells)
    R["emission"] = emitted
    R["output_mesh"] = str(mesh_out)
    R["output_mesh_hashes"] = mesh_hashes(mesh_out)
    R["n_invalidated_vertices"] = int(X.sum())
    R["n_removed_quads"] = int(removed.sum())
    R["n_retained_quads_out"] = int(kept.sum())
    R["area"] = {"A_input": a_in, "A_retained": a_kept, "A_excised": a_removed,
                 "identity_residual": a_kept + a_removed - a_in,
                 "denominator": ("every fraction is relative to the "
                                 "PRE-EXCISION BASE, the mesh actually cut")}
    R["retained_area_fraction"] = (a_kept / a_in) if a_in else None
    R["excised_area_fraction"] = (a_removed / a_in) if a_in else None

    after = census(mesh_out, f"{seg[:40]}_after", work, nv, nu, params)
    R["census_after"] = {f"d{d}": after["engine"][f"d{d}"] for d in (0, 1)}
    residual = after["transverse_total"]
    submultiset = not (after["multiset"] - before["multiset"])
    R["contact_multiset_submultiset_of_input"] = bool(submultiset)
    clean = bool(after["clean"] and submultiset)
    finish(status=("clean" if clean else "residual_contacts"),
           residual_transverse=int(residual),
           claimed_clean=clean,
           geometry_status=(GEOMETRY_STATUS_CLEAN if clean else
                            "residual_transverse_contacts"),
           clean_claim_authority=(
               "the emitted arrays were RELOADED FROM DISK and recensused by "
               "engines/selfcross under both diagonals"),
           note=(None if clean else
                 "the output still carries transverse contacts; NO clean "
                 "claim is made for it"))

    return {"segment": seg, "status": R["status"], "certificate": cert_path,
            "mesh_out": mesh_out, "removed_quads": int(removed.sum()),
            "invalidated_vertices": int(X.sum()),
            "retained_fraction": R["retained_area_fraction"],
            "before": R["census_before"], "after": R["census_after"],
            "clean": clean, "displacement": base is not None,
            "selection_status": sel.get("selection_status"),
            "wall_seconds": R["wall_seconds"], "out": out}
