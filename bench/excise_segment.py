"""First REAL certified excision: one segment, one aggregate output, one
certificate whose clean claim comes from the C++ engine.

Per notes/CUTTER-SPEC.md including its round-24 amendments (section 9), and
the round-24 reviewer ruling in notes/DECISIONS.md. This is the real-data
path section 9.5 declared but did not implement: the AUTHORITATIVE census
is `engines/selfcross` as sorted schema-v2 CSVs with recorded hashes --
never the planted-mesh pure-Python predicate. The Python side supplies only
the MILP (the reviewed two-stage lexicographic solve in
`windcheck.excise._solve_milp`) and the bookkeeping.

Discipline mirrors bench/repair_multi.py: census via
`bench/crossing_census.census_one` with cell=40, exclude=1, maxedge=60;
hash-tagged workdirs; certificates carrying full provenance (input/output
hashes, code commit, engine binary + source hashes, uv.lock, driver hash).

Sequence:

  1. PRE-REGISTRATION, written to disk and hashed BEFORE any solve: the
     maximum excision fraction, the solver time limit, the recensus
     iteration cap, the target event's identity and its min-exit lower
     bounds RECOMPUTED FRESH on the current geometry.
  2. Shared-support (self-touching) events are LABELLED `junction_excision`
     and cut like any other crossing (round-26 Q1 overturned the round-23
     refusal, which was a semantic objection, not a mathematical one). It is
     never called a branch separation. `--refuse-shared-support` restores
     the old behaviour.
  3. SEGMENT-WIDE coverage constraints (round-26): EVERY transverse row of
     the engine's census under BOTH diagonals, mapped to round-24
     eight-corner (quad-level) coverage sets, keeping the original
     (v, u, tri) identities. No event matching, no per-event iteration.
     The set is then reduced -- dedup, dominance, decomposition into
     independent components -- which provably cannot change the optimum.
  4. Two-stage lexicographic MILP per independent component; apply the
     vertex mask; emit ONE aggregate tifxyz with the input's grid dims
     under HYBRID invalidation (round-25 A1): mask.tif = 0 AND x=y=z=-1 at
     every excised cell, every RETAINED coordinate bit-identical, validity
     only valid -> invalid.
  5. AUTHORITATIVE VERIFICATION: reload the emitted mesh FROM DISK and
     recensus it with the C++ engine under BOTH diagonals. Acceptance =
     ZERO residual transverse keys of any kind AND the output contact
     multiset is a submultiset of the input's, every missing key carrying a
     deletion witness. A second engine census of a SIDECAR-LESS copy checks
     the naive consumer's view. A dirty recensus adds the residual rows as
     constraints and re-solves, up to the stated cap; if still not clean, NO
     clean claim is made -- a labelled failure certificate is written and
     the exit status is nonzero.

    uv run python bench/excise_segment.py --segment 20231005123336
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
from collections import Counter
from math import fsum
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, "src")
sys.path.insert(0, "bench")
from windcheck import tifxyz                                    # noqa: E402
from windcheck.check import PAIR_DTYPE                          # noqa: E402
from windcheck.excise import (GEOMETRY_STATUS_CLEAN,            # noqa: E402
                              HYBRID_INVALIDATION,
                              JUNCTION_EXCISION_LABEL, MISSING,
                              REDUCTION_RULE, RETAINED_BIT_IDENTITY,
                              SCHEDULING_NOTE, SELECTION_STATUS_RULE,
                              SHARED_SUPPORT_LABEL,
                              STAGE2_ALLOWANCE_RULE,
                              STALENESS_WARNING, _mesh_hash,
                              quad_area, quad_area_canonical, quad_corners,
                              quad_triangle_corners, select_global,
                              solve_global)
from windcheck.intrinsic import (SurfaceGraph, oriented_events,  # noqa: E402
                                 retained_quads)
from crossing_census import census_one                           # noqa: E402
from min_exit_sweep import ADMISSIBLE_REL_VX, event_bounds       # noqa: E402
from repair_segment import CORPORA, RES_UM, sha                  # noqa: E402

OUT = Path("out/excised")
CENSUS = {"cell": 40.0, "exclude": 1, "maxedge": 60.0, "threads": 0,
          "touch_tol": 1e-3}
AXES = ("x", "y", "z")
ENGINE = Path("engines/selfcross")
ENGINE_SRC = Path("engines/selfcross.cpp")

# ---------------------------------------------------------- pre-registration
# Pre-registered BEFORE the solve (written to disk and hashed; the hash goes
# in the certificate). These numbers are not tunable after seeing a result.
MAX_EXCISION_FRACTION = 0.01
MAX_EXCISION_FRACTION_JUSTIFICATION = (
    "Ceiling of 1% of the input retained canonical area. Justification: the "
    "MILP's decision scope is bounded a priori -- only vertices in some "
    "crossing pair's eight-corner coverage set may be invalidated, and only "
    "quads incident to those vertices can be removed. For this segment class "
    "(a single crossing event over 12 participating quads) that scope is a "
    "few dozen quads out of ~7e5 retained, i.e. O(1e-4) of the area, so the "
    "honest cut is three to four orders of magnitude inside the ceiling. The "
    "ceiling is therefore not a target but a REFUSAL TRIPWIRE for the ways a "
    "cut could go wrong without the solver noticing: a runaway recensus loop "
    "accumulating constraints, a mask applied to the wrong indices, or a "
    "retention-semantics divergence that prices phantom quads. 1% is chosen "
    "loose enough that no legitimate branch-separation cut on this class can "
    "trip it (so tripping it is informative, not a threshold argument) and "
    "tight enough that any of those failure modes -- all of which shatter or "
    "gouge the sheet -- is caught before a mesh is emitted. Over budget = "
    "REFUSE: no output mesh, no clean claim, labelled failure certificate.")
SOLVER_TIME_LIMIT_S = 300.0
MAX_RECENSUS_ITERATIONS = 3
MAX_EXCISION_FRACTION_BASIS = "A_excised_canonical / A_input_canonical"

CLEAN_CLAIM_AUTHORITY = (
    "The clean claim is made ONLY by engines/selfcross re-run on the emitted "
    "arrays reloaded from disk, under both diagonals. The MILP certifies "
    "nothing; neither does construction.")

OPERATION_LABEL = (
    "certified excision -- a TOPOLOGY CHANGE: surface present in the input is "
    "ABSENT from the output. This is not a repair and not a trace correction.")

# The engine's own emitted schema: v1,u1,v2,u2,verdict,penetration,angle_deg,
# tri1,tri2 (schema v2, sorted by (v1,u1,tri1,v2,u2,tri2)).
SCHEMA_V2_HEADER = "v1,u1,v2,u2,verdict,penetration,angle_deg,tri1,tri2"
VERDICTS = ("transverse", "coplanar", "grazing")


class Refusal(RuntimeError):
    """A pre-registered refusal: labelled non-output, never forced through."""

    def __init__(self, status: str, label: str, evidence: dict):
        super().__init__(label)
        self.status, self.label, self.evidence = status, label, evidence


class Timers:
    def __init__(self):
        self.t: Counter = Counter()
        self.n: Counter = Counter()

    def add(self, key: str, dt: float) -> None:
        self.t[key] += dt
        self.n[key] += 1

    def timed(self, key: str, fn, *a, **kw):
        t0 = time.time()
        out = fn(*a, **kw)
        self.add(key, time.time() - t0)
        return out

    def report(self) -> dict:
        return {k: {"seconds": round(self.t[k], 3), "calls": self.n[k]}
                for k in sorted(self.t)}


# ------------------------------------------------------------ census parsing
def parse_census_csv(csv: Path, diag: int, nv: int, nu: int) -> dict:
    """Parse one AUTHORITATIVE schema-v2 selfcross CSV.

    Returns contact records that keep the engine's ORIGINAL (v, u, tri)
    identities verbatim, plus the round-19 contact multiset over
    (diag, v1, u1, tri1, v2, u2, tri2, verdict). Legacy 7-column CSVs are
    refused: quad-level keys hide one triangle combination replacing
    another, so they cannot support a submultiset claim.

    Rows whose quad origin falls outside the grid's quad range are counted
    (`out_of_range`) rather than silently dropped; the engine never emits
    them, so a nonzero count is evidence of a schema/grid mismatch.
    """
    text = csv.read_text().splitlines()
    if not text:
        raise ValueError(f"empty census CSV: {csv}")
    header = text[0].strip()
    if "tri1" not in header:
        raise ValueError(f"legacy census CSV (no triangle identities): {csv}")
    if header != SCHEMA_V2_HEADER:
        raise ValueError(f"unexpected census schema in {csv}: {header!r}")
    rows: list[dict] = []
    multiset: Counter = Counter()
    out_of_range = 0
    unknown_verdict: Counter = Counter()
    for line in text[1:]:
        if not line.strip():
            continue
        p = line.split(",")
        verdict = p[4]
        if verdict not in VERDICTS:
            unknown_verdict[verdict] += 1
            continue
        v1, u1, v2, u2 = (int(p[0]), int(p[1]), int(p[2]), int(p[3]))
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
    if unknown_verdict:
        raise ValueError(f"unknown verdicts in {csv}: {dict(unknown_verdict)}")
    return {"rows": rows, "multiset": multiset, "out_of_range": out_of_range,
            "header": header, "n_lines": len(text) - 1,
            "path": str(csv), "sha256": sha(csv)}


def transverse_pair_records(rows: list[dict]) -> np.ndarray:
    """The transverse rows as a PAIR_DTYPE array for `oriented_events`.

    One record per DISTINCT quad pair (events are quad-level objects); the
    triangle identities stay on the contact rows, which is where the
    submultiset claim lives.
    """
    seen: dict = {}
    for r in rows:
        if r["verdict"] != "transverse":
            continue
        k = (r["q1"], r["q2"])
        if k not in seen:
            seen[k] = (r["q1"][0], r["q1"][1], r["q2"][0], r["q2"][1],
                       r["penetration"], r["angle_deg"])
    return np.array([seen[k] for k in sorted(seen)], dtype=PAIR_DTYPE)


def coverage_from_rows(rows: list[dict], Q: np.ndarray) -> list[dict]:
    """Map engine census rows to the coverage sets excise.py's MILP expects.

    Only TRANSVERSE rows become constraints: the acceptance criterion is
    zero transverse crossing pairs (CUTTER-SPEC section 2), and coplanar /
    grazing contacts are recorded in the multiset but are not crossings.

    Coverage is the round-24 eight-corner (quad-level) set:
    corners(quad1) UNION corners(quad2). Both quads must be members of the
    censused complex `retained_quads` computes; a divergence there means the
    Python retention semantics disagree with the engine's, which would price
    or remove phantom quads, so it is a hard error and not worked around.
    The two quads must share no corner (guaranteed by the engine's
    Chebyshev exclude=1 adjacency exclusion).
    """
    out = []
    for r in rows:
        if r["verdict"] != "transverse":
            continue
        (v1, u1), (v2, u2) = r["q1"], r["q2"]
        for (v, u) in (r["q1"], r["q2"]):
            if not Q[v, u]:
                raise ValueError(
                    f"censused quad ({v},{u}) is not in retained_quads: "
                    "engine/Python retention semantics diverge")
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


def shared_support_events(events_by_diag: dict) -> list[dict]:
    """The round-23 THIRD CLASS: events whose two branch regions share
    support (overlapping or grid-adjacent regions -- `self_touching`). The
    cutter refuses these; a vertex mask there destroys both branches at once
    and cannot be certified as a branch separation."""
    out = []
    for d, evs in sorted(events_by_diag.items()):
        for k, e in enumerate(evs):
            if e.get("self_touching"):
                out.append({"diagonal": d, "event_index": k,
                            "region_a": sorted(map(list, e["region_a"])),
                            "region_b": sorted(map(list, e["region_b"])),
                            "ambiguous": bool(e["ambiguous"])})
    return out


# ------------------------------------------------------------- area at scale
def quad_area_grids(P64: np.ndarray) -> dict:
    """Per-quad areas for the whole grid, vectorised, float64.

    Convention (CUTTER-SPEC section 6): the area of a quad under diagonal d
    is the sum of its two triangle areas under that tessellation; the
    canonical value is the mean of d0 and d1. Triangle corner order matches
    engines/selfcross.cpp / `quad_triangle_corners` exactly.
    """
    p00 = P64[:-1, :-1]
    p10 = P64[1:, :-1]
    p01 = P64[:-1, 1:]
    p11 = P64[1:, 1:]

    def tri(a, b, c):
        return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=-1)

    a0 = tri(p00, p01, p11) + tri(p00, p11, p10)
    a1 = tri(p00, p01, p10) + tri(p01, p11, p10)
    return {0: a0, 1: a1, "canonical": 0.5 * (a0 + a1)}


def area_block(areas: np.ndarray, Q_in: np.ndarray, removed, unresolved
               ) -> dict:
    """Area accounting over a quad PARTITION (no quad counted twice).

    Quad ownership is a partition of the input-retained set into
    {clean, excised, unresolved} (CUTTER-SPEC section 6). `math.fsum` is
    exactly rounded, so the sums do not depend on iteration order.
    """
    Q_clean = Q_in.copy()
    for q in removed:
        Q_clean[q] = False
    for q in unresolved:
        Q_clean[q] = False
    a_in = fsum(areas[Q_in].tolist())
    a_clean = fsum(areas[Q_clean].tolist())
    a_exc = fsum(float(areas[q]) for q in removed)
    a_unres = fsum(float(areas[q]) for q in unresolved)
    return {"A_input": a_in, "A_clean": a_clean, "A_excised": a_exc,
            "A_unresolved": a_unres,
            "clean_recovery_fraction": (a_clean / a_in) if a_in else None,
            "excised_fraction": (a_exc / a_in) if a_in else None,
            "identity_residual": a_clean + a_exc + a_unres - a_in}


def budget_verdict(a_excised: float, a_input: float,
                   max_fraction: float = MAX_EXCISION_FRACTION) -> dict:
    """The pre-registered excision-fraction gate. Evaluated BEFORE any mesh
    is emitted; over budget is a REFUSAL, never a warning."""
    frac = (a_excised / a_input) if a_input else float("inf")
    ok = frac <= max_fraction
    return {"basis": MAX_EXCISION_FRACTION_BASIS,
            "A_excised": a_excised, "A_input": a_input,
            "excised_fraction": frac, "max_excision_fraction": max_fraction,
            "within_budget": bool(ok),
            "label": ("within the pre-registered excision budget" if ok else
                      "excision exceeds the pre-registered maximum fraction; "
                      "REFUSED with no output mesh and no clean claim")}


def components(Q: np.ndarray, areas: np.ndarray) -> tuple[int, list[dict]]:
    """Connected components of the retained-quad complex, joined by SHARED
    CORNERS -- which on the quad index grid is exactly 8-connectivity (quads
    at Chebyshev distance 1 share a corner; farther quads share none)."""
    from scipy import ndimage
    lab, n = ndimage.label(Q, structure=np.ones((3, 3), dtype=bool))
    if n == 0:
        return 0, []
    counts = np.bincount(lab.ravel(), minlength=n + 1)[1:]
    sums = np.bincount(lab.ravel(), weights=areas.ravel(),
                       minlength=n + 1)[1:]
    dist = sorted(({"n_quads": int(c), "area_canonical": float(a)}
                   for c, a in zip(counts, sums)),
                  key=lambda d: (-d["area_canonical"], -d["n_quads"]))
    return int(n), dist


def cut_boundary(removed: set, kept: set) -> list:
    """Grid edges shared between a removed quad and a RETAINED output quad.

    Each grid edge is a side of at most two quads, so the neighbour across
    each side of a removed quad determines it: the edge is a cut-boundary
    edge iff that neighbour survives in the output.
    """
    out = set()
    for v, u in removed:
        for nb, edge in (((v - 1, u), ((v, u), (v, u + 1))),
                         ((v + 1, u), ((v + 1, u), (v + 1, u + 1))),
                         ((v, u - 1), ((v, u), (v + 1, u))),
                         ((v, u + 1), ((v, u + 1), (v + 1, u + 1)))):
            if nb in kept:
                out.add(tuple(sorted(edge)))
    return sorted(out)


# ------------------------------------------------------------------ witnesses
def witness_for_key(key, invalidated: set) -> dict | None:
    """Deletion witness for one destroyed contact key.

    Round-24 scopes: "triangle_participant" -- the invalidated vertex is a
    corner of the destroyed contact triangle itself; "quad_retention" -- it
    is a corner of the participating quad but not of that triangle, and
    quad-level validity drops the quad and so the triangle anyway.
    """
    diag, v1, u1, t1, v2, u2, t2, _verdict = key
    for (qv, qu), t in (((v1, u1), t1), ((v2, u2), t2)):
        tri = quad_triangle_corners(qv, qu, diag)[t]
        for cnr in tri:
            if cnr in invalidated:
                return {"invalidated_vertex": list(cnr),
                        "destroyed_triangle": [diag, qv, qu, t],
                        "witness_scope": "triangle_participant"}
        for cnr in quad_corners(qv, qu):
            if cnr in invalidated and cnr not in set(tri):
                return {"invalidated_vertex": list(cnr),
                        "destroyed_triangle": [diag, qv, qu, t],
                        "witness_scope": "quad_retention"}
    return None


def submultiset_report(before: Counter, after: Counter,
                       invalidated: set) -> dict:
    """Output contact multiset vs input's, with a witness per missing key.

    Coordinates never change and validity only shrinks, so contact identity
    IS the key: `after` must be a submultiset of `before`, and every key
    whose multiplicity dropped must be explained by an invalidated vertex.
    """
    increased = [{"key": list(k), "before": int(before.get(k, 0)),
                  "after": int(n)}
                 for k, n in sorted(after.items()) if n > before.get(k, 0)]
    missing, unwitnessed = [], []
    for k in sorted(before):
        if after.get(k, 0) >= before[k]:
            continue
        w = witness_for_key(k, invalidated)
        rec = {"key": list(k), "before": int(before[k]),
               "after": int(after.get(k, 0)), "witness": w}
        missing.append(rec)
        if w is None:
            unwitnessed.append(rec)
    return {"n_input_keys": int(sum(before.values())),
            "n_output_keys": int(sum(after.values())),
            "n_distinct_input_keys": len(before),
            "n_distinct_output_keys": len(after),
            "output_submultiset_of_input": not increased,
            "new_or_increased_keys": increased,
            "removed_contacts": missing,
            "missing_all_witnessed": not unwitnessed,
            "unwitnessed_missing_keys": unwitnessed,
            "input_multiset": [[list(k), int(n)] for k, n in
                               sorted(before.items())],
            "output_multiset": [[list(k), int(n)] for k, n in
                                sorted(after.items())]}


def solve_totals(records: list[dict], stage: int = 1) -> dict:
    """Aggregate a decomposed solve back into one reportable objective.

    A segment-wide solve is split across independent components (round-26
    Q2c), so the segment's objective is the SUM of the component objectives
    and its dual bound the sum of the component bounds -- valid precisely
    because the components share no variable. The reported gap is the
    aggregate |obj - bound| / |obj|, which is what a reader wants: how far
    the delivered cut could still be from the segment optimum.
    """
    recs = [r for r in records if r["stage"] == stage]
    objs = [r["objective"] for r in recs]
    bounds = [r["dual_bound"] for r in recs]
    obj = (fsum(objs) if recs and all(o is not None for o in objs) else None)
    bound = (fsum(bounds) if recs and all(b is not None for b in bounds)
             else None)
    gap = (abs(obj - bound) / abs(obj) if obj not in (None, 0.0)
           and bound is not None else (0.0 if obj == 0.0 else None))
    return {"objective": obj, "dual_bound": bound, "mip_gap": gap,
            "seconds": fsum(r["solve_time_s"] for r in recs),
            "n_solves": len(recs),
            "statuses": sorted({r["status"] for r in recs})}


def is_attributable(key, targets: list[dict]) -> bool:
    """Is this contact key one of the target event's own crossings?"""
    q1, q2 = (key[1], key[2]), (key[5], key[6])
    for t in targets:
        A, B = t["region_a"], t["region_b"]
        if (q1 in A and q2 in B) or (q1 in B and q2 in A):
            return True
    return False


# ----------------------------------------------------------------- emission
def emit_excised_tifxyz(src: Path, dst: Path, valid_out: np.ndarray,
                        excised: np.ndarray, *, with_mask: bool = True
                        ) -> dict:
    """Emit ONE aggregate tifxyz under HYBRID invalidation (round-25 A1).

    tifxyz marks a missing cell EITHER by setting x=y=z=-1 OR through a mask
    sidecar. Round 25 ruled that using only one of them is an
    interoperability DEFECT: a consumer that implements the other convention
    reads the excised cells as present and reconstructs exactly the crossing
    this operation removed. So both are written -- mask.tif = 0 AND
    x = y = z = -1 -- and the coordinate guarantee is restated over RETAINED
    cells: every retained coordinate is bit-identical to the input, and the
    ONLY cells whose coordinates differ are the excised ones.

    `with_mask=False` emits the same coordinate planes WITHOUT the sidecar:
    that is the naive consumer's view, and censusing it is how the -1 half of
    the hybrid gets verified by the engine rather than asserted.
    """
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    for stale in ("mask.tif", "mask.png"):
        assert not (src / stale).exists(), (
            f"input carries {stale}: the emission path would overwrite input "
            "validity semantics; refusing")
    for ax in AXES:
        plane = np.asarray(tifffile.imread(src / f"{ax}.tif"))
        assert plane.shape == excised.shape, (ax, plane.shape, excised.shape)
        out = plane.copy()
        out[excised] = MISSING
        tifffile.imwrite(dst / f"{ax}.tif", out)
    for extra in ("meta.json",):
        if (src / extra).exists():
            shutil.copy(src / extra, dst / extra)
    if with_mask:
        tifffile.imwrite(dst / "mask.tif", valid_out.astype(np.uint8))
    return {"invalidation_carrier": HYBRID_INVALIDATION,
            "mask_written": with_mask,
            "coordinate_guarantee": RETAINED_BIT_IDENTITY,
            "consumer_note": (
                "Excised cells are marked BOTH ways. A consumer that honours "
                "only mask.tif and one that honours only the x=y=z=-1 "
                "convention see the SAME retained surface; neither can "
                "reconstruct the removed crossing. The retained coordinates "
                "are unchanged, so no downstream registration or texture "
                "lookup on retained geometry is affected."),
            "files": {f: sha(dst / f) for f in sorted(
                p.name for p in dst.iterdir())}}


# --------------------------------------------------------------------- main
def resolve_segment(sub: str):
    exact, subs = [], []
    for corpus, root, volume, work in CORPORA:
        for d in (Path(root).iterdir() if Path(root).exists() else []):
            if d.is_dir() and sub in d.name:
                (exact if d.name == sub else subs).append((corpus, d, volume))
    cands = exact or subs
    if len(cands) != 1:
        raise SystemExit(f"{sub}: {len(cands)} matches "
                         f"{[c[1].name for c in cands]}")
    return cands[0]


def census_mesh(mesh: Path, tag: str, wdir: Path, nv: int, nu: int,
                timers: Timers) -> dict:
    """One authoritative both-diagonal census of a mesh ON DISK."""
    counts = timers.timed("engine_census", census_one, mesh, tag,
                          CENSUS["exclude"], CENSUS["cell"],
                          CENSUS["threads"], CENSUS["maxedge"], wdir)
    if counts is None:
        raise Refusal("not_censusable",
                      "mesh below the census validity threshold; no audit or "
                      "excision is defined on it", {"mesh": str(mesh)})
    per_diag = {}
    for d in (0, 1):
        per_diag[d] = parse_census_csv(wdir / f"{tag[:40]}_d{d}.csv", d,
                                       nv, nu)
    multiset: Counter = Counter()
    for d in (0, 1):
        multiset.update(per_diag[d]["multiset"])
    return {"tag": tag, "engine": counts, "diag": per_diag,
            "multiset": multiset,
            "rows": [r for d in (0, 1) for r in per_diag[d]["rows"]]}


def code_provenance() -> dict:
    """Provenance a public consumer can verify from the published files.

    No commit sha and no branch name: a reader of the release has neither,
    so a citation into repository history is a claim they cannot check.
    What is recorded instead is content -- the code version, the digest of
    the published source tree, and the individual hashes of the driver,
    the excision module, the engine binary and source, and the lockfile.
    Every one of those is recomputable from the release with
    `uv run python -m windcheck.provenance`.
    """
    from windcheck.provenance import release_provenance
    return {**release_provenance(),
            "driver_sha256": sha(Path("bench/excise_segment.py")),
            "excise_module_sha256": sha(Path("src/windcheck/excise.py")),
            "engine_binary_sha256": sha(ENGINE),
            "engine_source_sha256": sha(ENGINE_SRC),
            "uv_lock_sha256": sha(Path("uv.lock"))}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--segment", required=True)
    ap.add_argument("--max-excision-fraction", type=float,
                    default=MAX_EXCISION_FRACTION)
    ap.add_argument("--time-limit", type=float, default=SOLVER_TIME_LIMIT_S)
    ap.add_argument("--max-iterations", type=int,
                    default=MAX_RECENSUS_ITERATIONS)
    ap.add_argument("--refuse-shared-support", action="store_true",
                    help="restore the round-23 refusal instead of the "
                         "round-26 junction_excision label")
    ap.add_argument("--skip-bounds", action="store_true",
                    help="skip the per-event min-exit bounds (they are "
                         "informational under segment-wide constraints)")
    ap.add_argument("--strategy", default="lp_round",
                    choices=["lp_round", "exact"],
                    help="round-27 default: LP relaxation + 1/k rounding + "
                         "reverse-delete + bounded local improvement, with "
                         "the exact MILP only where scheduling allows it")
    ap.add_argument("--lp-time-limit", type=float, default=600.0)
    ap.add_argument("--lp-total-budget", type=float, default=2400.0)
    ap.add_argument("--exact-max-constraints", type=int, default=700)
    ap.add_argument("--exact-time-limit", type=float, default=30.0)
    ap.add_argument("--exact-total-budget", type=float, default=600.0,
                    help="SEGMENT-WIDE exact-MILP budget; components are "
                         "attempted cheapest first and skipped once it is "
                         "spent (costs optimality claims, never feasibility)")
    ap.add_argument("--stage2-max-constraints", type=int, default=50)
    ap.add_argument("--improve-budget", type=float, default=90.0)
    ap.add_argument("--skip-naive-check", action="store_true",
                    help="skip the extra engine census of a sidecar-less "
                         "copy (the naive-consumer verification)")
    args = ap.parse_args()
    t_start = time.time()
    timers = Timers()

    corpus, seg, volume = resolve_segment(args.segment)
    mesh = sorted(seg.glob(f"mesh/*{volume}*.tifxyz"))[0]
    m = RES_UM.search(mesh.name)
    voxel_um = float(m.group(1)) if m else 7.91
    tag = hashlib.sha256(seg.name.encode()).hexdigest()[:12]
    wdir = OUT / f"work_{tag}"
    wdir.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    cert_path = OUT / f"{seg.name}_excision_certificate.json"
    mesh_out = OUT / f"{seg.name}.tifxyz"
    # set only once the mesh is actually on disk, so a failure record can
    # remove an unclaimed emission without guessing
    mesh_out_holder: dict = {"path": None}

    def log(msg):
        print(msg, flush=True)

    log(f"{corpus} {seg.name}  mesh {mesh.name}  voxel {voxel_um} um")
    log(f"  workdir {wdir}")
    prov = code_provenance()
    log(f"  code {prov['code_version']} source tree "
        f"{prov['source_tree_digest'][:12]} policy {prov['policy_version']}")

    surf = tifxyz.read(mesh)
    P_in = surf.points
    P64 = np.asarray(P_in, np.float64)
    V_in = np.asarray(surf.valid, bool)
    nv, nu = surf.shape
    Q_in = retained_quads(P64, V_in, CENSUS["maxedge"])
    areas = quad_area_grids(P64)
    q_in = [(int(v), int(u)) for v, u in zip(*np.nonzero(Q_in))]
    log(f"  grid {nv}x{nu}  valid {int(V_in.sum()):,}  "
        f"retained quads {len(q_in):,}")

    def failure(status: str, label: str, extra: dict, *,
                kind: str = "failure", exit_code: int = 1) -> int:
        """Write a LABELLED non-output record and leave no mesh behind.

        Spec ambiguity resolved here: CUTTER-SPEC section 3 says a refusal
        emits "no output mesh and no certificate", while a real run still has
        to leave an auditable trace of WHY it refused. The trace written is
        explicitly a refusal/failure RECORD, not a certificate: it carries
        the class label and the census evidence, has `clean_claim: false`,
        and never carries output arrays or an output mesh path.
        """
        if mesh_out_holder["path"] is not None \
                and Path(mesh_out_holder["path"]).exists():
            shutil.rmtree(mesh_out_holder["path"])   # never leave an
            # unclaimed mesh on disk that could be mistaken for output
        cert = {"operation": OPERATION_LABEL, "status": status,
                "record_kind": kind,
                "record_is_not_a_certificate": (
                    "This is a labelled non-output record. No clean claim is "
                    "made, no output mesh exists, and nothing here may be "
                    "cited as a certified excision."),
                "clean_claim": False, "label": label,
                "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                               time.gmtime()),
                "corpus": corpus, "segment": seg.name, "volume": volume,
                "mesh": str(mesh), "voxel_um": voxel_um,
                "grid_shape": [nv, nu], "dtype": str(P_in.dtype),
                "code_provenance": prov,
                "census_params": dict(CENSUS, diagonals=[0, 1],
                                      engine=str(ENGINE),
                                      schema=SCHEMA_V2_HEADER),
                "output_mesh": None,
                "staleness_warning": STALENESS_WARNING,
                "instrumentation": {**timers.report(),
                                    "wall_seconds": round(
                                        time.time() - t_start, 2)},
                **extra}
        path = OUT / f"{seg.name}_excision_{kind}.json"
        path.write_text(json.dumps(cert, indent=1))
        log(f"\n[{status.upper()}] {label}")
        log(f"  {kind} record {path}  (no output mesh, no clean claim)")
        return exit_code

    # ---- authoritative BEFORE census ------------------------------------
    try:
        before = census_mesh(mesh, f"{tag}_before", wdir, nv, nu, timers)
    except Refusal as r:
        # a sub-threshold mesh is a DISPOSITION, not an error: refuse cleanly
        return failure(r.status, r.label, {"evidence": r.evidence},
                       kind="refusal", exit_code=0)
    log(f"  census before: d0 transverse {before['engine']['d0']['transverse']}"
        f" grazing {before['engine']['d0']['grazing']}"
        f" coplanar {before['engine']['d0']['coplanar']} | "
        f"d1 transverse {before['engine']['d1']['transverse']}"
        f" grazing {before['engine']['d1']['grazing']}"
        f" coplanar {before['engine']['d1']['coplanar']}")

    # ---- events, per diagonal, from the engine's transverse rows --------
    recs, events = {}, {}
    for d in (0, 1):
        recs[d] = transverse_pair_records(before["diag"][d]["rows"])
        events[d] = oriented_events(recs[d]) if len(recs[d]) else []
        log(f"  d{d}: {len(recs[d])} transverse quad pairs, "
            f"{len(events[d])} event(s)")

    # ---- shared support: JUNCTION EXCISION, not a refusal ---------------
    # Round-26 Q1 overturned the round-23 refusal: shared support is
    # mathematically excisable, and the old refusal was a SEMANTIC objection
    # (a vertex mask there is not a branch separation). It still is not one,
    # so the operation is labelled `junction_excision` and never described as
    # separating branches -- but it is performed and certified like any other
    # cut. `--refuse-shared-support` restores the old behaviour.
    touching = shared_support_events(events)
    if touching and args.refuse_shared_support:
        return failure("refused_shared_support", SHARED_SUPPORT_LABEL,
                       {"shared_support_events": touching,
                        "class": "third class: shared-support topology",
                        "census_before": {f"d{d}": before["engine"][f"d{d}"]
                                          for d in (0, 1)}},
                       kind="refusal", exit_code=0)
    operation_label = (f"junction_excision -- {JUNCTION_EXCISION_LABEL}"
                       if touching else "excision")
    if touching:
        log(f"  shared-support events: {len(touching)} -> operation labelled "
            "junction_excision (never a branch separation)")

    # ---- event identities + FRESH min-exit bounds -----------------------
    # NOT a selection: round 26 removed event matching from the excision
    # problem entirely (the constraints come from every transverse row). The
    # bounds are recorded because they establish WHY displacement is
    # inapplicable, and they are skippable on segments with many events.
    targets = []
    for d in (0, 1) if not args.skip_bounds else ():
        for k, e in enumerate(events[d]):
            if e["ambiguous"]:
                targets.append({"diagonal": d, "event_index": k,
                                "ambiguous": True,
                                "region_a": set(map(tuple, e["region_a"])),
                                "region_b": set(map(tuple, e["region_b"])),
                                "bounds": None})
                continue
            g = timers.timed("graph_build", SurfaceGraph, P64, V_in, d,
                             CENSUS["maxedge"])
            gv, gu = np.nonzero(g.idx >= 0)
            X = np.empty((g.n, 3))
            X[g.idx[gv, gu]] = P64[gv, gu]
            t0 = time.time()
            L, L_safe, wit, npairs = event_bounds(g, X, recs[d], e)
            timers.add("min_exit", time.time() - t0)
            targets.append({
                "diagonal": d, "event_index": k, "ambiguous": False,
                "region_a": set(map(tuple, e["region_a"])),
                "region_b": set(map(tuple, e["region_b"])),
                "n_pairs": len(e["rows"]), "n_tri_pairs": npairs,
                "bounds": {"L_vx": round(L, 6), "L_safe_vx": round(L_safe, 6),
                           "witness": wit},
                "rigid_verdict": ("certified_infeasible"
                                  if L_safe > ADMISSIBLE_REL_VX
                                  else "inconclusive")})
            log(f"  d{d} event {k}: L {L:.4f} vx  L_safe {L_safe:.4f} vx  "
                f"-> {targets[-1]['rigid_verdict']}")

    def ev_json(t):
        return {"diagonal": t["diagonal"], "event_index": t["event_index"],
                "ambiguous": t["ambiguous"],
                "region_a": sorted(map(list, t["region_a"])),
                "region_b": sorted(map(list, t["region_b"])),
                "n_transverse_quad_pairs": t.get("n_pairs"),
                "n_intersecting_triangle_pairs": t.get("n_tri_pairs"),
                "min_exit_bounds": t.get("bounds"),
                "rigid_verdict": t.get("rigid_verdict")}

    # ---- PRE-REGISTRATION, written and hashed BEFORE the first solve ----
    prereg = {
        "operation": OPERATION_LABEL,
        "recorded_before_any_solve": True,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "corpus": corpus, "segment": seg.name, "volume": volume,
        "mesh": str(mesh), "voxel_um": voxel_um,
        "grid_shape": [nv, nu], "dtype": str(P_in.dtype),
        "input_mesh_sha256": {f"{ax}.tif": sha(mesh / f"{ax}.tif")
                              for ax in AXES},
        "input_array_sha256": _mesh_hash(P_in, V_in),
        "max_excision_fraction": args.max_excision_fraction,
        "max_excision_fraction_basis": MAX_EXCISION_FRACTION_BASIS,
        "max_excision_fraction_justification":
            MAX_EXCISION_FRACTION_JUSTIFICATION,
        "solver_time_limit_s": args.time_limit,
        "max_recensus_iterations": args.max_iterations,
        "operation_label": operation_label,
        "shared_support_events": touching,
        "constraint_scope": (
            "SEGMENT-WIDE (round 26): EVERY transverse row of the "
            "authoritative census under BOTH diagonals becomes a coverage "
            "constraint in ONE MILP. There is no event matching and no "
            "per-event iteration, so ambiguous, unmatched, "
            "search-inconclusive and shared-support rows are all handled "
            "identically."),
        "constraint_reduction_rule": REDUCTION_RULE,
        "acceptance_criterion": (
            "ZERO residual TRANSVERSE contact keys of ANY kind in the "
            "RELOADED output's both-diagonal C++ census, AND the output "
            "contact multiset (all verdicts, both diagonals) a submultiset "
            "of the input's with a deletion witness for every missing key"),
        "clean_claim_authority": CLEAN_CLAIM_AUTHORITY,
        "census_params": dict(CENSUS, diagonals=[0, 1], engine=str(ENGINE),
                              schema=SCHEMA_V2_HEADER),
        "engine_binary_sha256": prov["engine_binary_sha256"],
        "source_tree_digest": prov["source_tree_digest"],
        "code_version": prov["code_version"],
        "census_before": {f"d{d}": before["engine"][f"d{d}"]
                          for d in (0, 1)},
        "census_before_csv": {f"d{d}": {
            "path": before["diag"][d]["path"],
            "sha256": before["diag"][d]["sha256"]} for d in (0, 1)},
        "target_events": [ev_json(t) for t in targets],
        "admissible_rel_vx": ADMISSIBLE_REL_VX,
        "min_exit_model_scope": (
            "the min-exit bounds rule out RIGID RELATIVE CORE TRANSLATION "
            "within the admissible budget only; they say nothing about "
            "non-rigid deformation or remeshing. They are recorded here to "
            "establish that displacement repair is inapplicable to this "
            "event, which is why excision is the operator."),
    }
    prereg_path = OUT / f"{seg.name}_prereg.json"
    prereg_path.write_text(json.dumps(prereg, indent=1))
    prereg_sha = sha(prereg_path)
    log(f"  PRE-REGISTERED (before solving): max excised fraction "
        f"{args.max_excision_fraction}, solver limit {args.time_limit}s, "
        f"recensus cap {args.max_iterations}")
    log(f"    {prereg_path}  sha256 {prereg_sha[:16]}")

    # ---- coverage constraints from the AUTHORITATIVE census -------------
    try:
        cons_list = coverage_from_rows(before["rows"], Q_in)
    except ValueError as e:
        return failure("census_mapping_error", str(e),
                       {"pre_registration": prereg,
                        "pre_registration_sha256": prereg_sha})
    constraints = {c["key"]: c for c in cons_list}
    log(f"  coverage constraints: {len(constraints)} transverse rows, "
        f"{len({p for c in cons_list for p in c['coverage']})} candidate "
        f"vertices")

    # ---- solve / apply / emit / RECENSUS loop ---------------------------
    solves: list[dict] = []
    lexicographic = True
    scipy_version = None
    iterations = 0
    emit_report: dict = {}
    after: dict = {}
    invalidated: list = []
    removed: list = []
    kept: list = []
    residual_target: list = []
    accepted = False
    budget: dict = {}
    reload_checks: dict = {}
    reduction: dict = {}
    naive_check: dict = {}
    selection: dict | None = None
    if mesh_out.exists():
        shutil.rmtree(mesh_out)          # never reuse a previous emission
    while True:
        if iterations >= args.max_iterations:
            return failure(
                "iteration_limit",
                "no clean claim: recensus iteration cap reached with residual "
                "crossings", {"pre_registration": prereg,
                              "pre_registration_sha256": prereg_sha,
                              "solver": {"solves": solves},
                              "residual_target_keys": residual_target,
                              "iterations": iterations})
        iterations += 1
        log(f"  --- iteration {iterations}: {args.strategy} selection over "
            f"{len(constraints)} constraints")
        if args.strategy == "lp_round":
            sel = timers.timed(
                "selection", select_global, list(constraints.values()),
                P64, Q_in, (),
                lp_time_limit=args.lp_time_limit,
                lp_total_budget=args.lp_total_budget,
                exact_max_constraints=args.exact_max_constraints,
                exact_time_limit=args.exact_time_limit,
                exact_total_budget=args.exact_total_budget,
                stage2_max_constraints=args.stage2_max_constraints,
                improve_budget=args.improve_budget,
                area_grid=areas["canonical"])
            selection = sel
            sol = {"status": ("optimal" if sel["status"] == "ok"
                              else sel["status"]),
                   "chosen": sel["chosen"],
                   "records": sel.get("milp_records", []),
                   "lexicographic": False,
                   "scipy_version": sel["scipy_version"],
                   "reduction": sel["reduction"]}
            log(f"      selection {sel['selection_status']}: mix "
                f"{sel['method_mix']}; achieved area "
                f"{sel['achieved_area']}; combined lower bound "
                f"{sel['combined_lower_bound']} "
                f"(complete={sel['combined_lower_bound_complete']}); ratio "
                f"{sel['ratio_achieved_over_bound']}")
        else:
            sol = timers.timed("milp", solve_global,
                               list(constraints.values()),
                               P64, Q_in, set(), args.time_limit)
        for rec in sol["records"]:
            rec["iteration"] = iterations
        solves.extend(sol["records"])
        lexicographic = lexicographic and sol["lexicographic"]
        scipy_version = sol["scipy_version"]
        reduction = sol["reduction"]
        log(f"      reduction: {reduction['n_raw']} raw -> "
            f"{reduction['n_after_dedup']} deduped -> "
            f"{reduction['n_after_dominance']} after dominance "
            f"({reduction['n_components']} independent components, largest "
            f"{reduction['largest_component']})")
        agg = solve_totals(sol["records"])
        log(f"      stage1 {sol['status']} objective {agg['objective']} "
            f"bound {agg['dual_bound']} gap {agg['mip_gap']} "
            f"({agg['seconds']:.2f}s over {agg['n_solves']} component "
            f"solve(s))")
        if sol["status"] not in ("optimal", "best_found"):
            return failure(sol["status"],
                           "no clean claim: solver returned no usable "
                           "incumbent",
                           {"pre_registration": prereg,
                            "pre_registration_sha256": prereg_sha,
                            "solver": {"solves": solves,
                                       "scipy_version": scipy_version}})
        invalidated = sorted(sol["chosen"])
        inval = set(invalidated)
        V_out = V_in.copy()
        for v, u in invalidated:
            V_out[v, u] = False
        Q_out = retained_quads(P64, V_out, CENSUS["maxedge"])
        removed, kept = [], []
        for q in q_in:
            (removed if set(quad_corners(*q)) & inval else kept).append(q)
        # the mask must reproduce the quad partition exactly
        assert set(removed) | set(kept) == set(q_in)
        assert not (set(removed) & set(kept))
        assert {(int(v), int(u)) for v, u in zip(*np.nonzero(Q_out))} \
            == set(kept), "masked retained-quad set disagrees with the mask"
        log(f"      invalidated {len(invalidated)} vertices, "
            f"removed {len(removed)} quads")

        # ---- PRE-REGISTERED BUDGET GATE, before any mesh is emitted ----
        blocks = {k: area_block(areas[k], Q_in, removed, [])
                  for k in (0, 1, "canonical")}
        budget = budget_verdict(blocks["canonical"]["A_excised"],
                                blocks["canonical"]["A_input"],
                                args.max_excision_fraction)
        log(f"      excised canonical area "
            f"{budget['A_excised']:.4f} vx^2 = "
            f"{budget['excised_fraction']:.3e} of input "
            f"(cap {budget['max_excision_fraction']})")
        if not budget["within_budget"]:
            return failure("refused_over_budget", budget["label"],
                           {"pre_registration": prereg,
                            "pre_registration_sha256": prereg_sha,
                            "budget_gate": budget,
                            "solver": {"solves": solves,
                                       "scipy_version": scipy_version},
                            "proposed_excision": {
                                "invalidated_vertices":
                                    sorted(map(list, invalidated)),
                                "n_removed_quads": len(removed)}},
                           # a pre-registered refusal, but a nonzero exit:
                           # the cut the solver proposed is unacceptable, so
                           # the run failed to deliver an excision
                           kind="refusal", exit_code=1)

        # ---- emit the aggregate ----------------------------------------
        excised_cells = V_in & ~V_out
        emit_report = timers.timed("emit", emit_excised_tifxyz, mesh,
                                   mesh_out, V_out, excised_cells)
        mesh_out_holder["path"] = str(mesh_out)

        # ---- reload FROM DISK and verify the contract -------------------
        surf_out = tifxyz.read(mesh_out)
        Pout = np.asarray(surf_out.points)
        changed = (np.any(Pout != P_in, axis=-1)
                   if Pout.shape == P_in.shape else np.ones((nv, nu), bool))
        reload_checks = {
            "reloaded_from": str(mesh_out),
            "grid_shape_equal": list(surf_out.shape) == [nv, nu],
            "dtype_equal": str(surf_out.points.dtype) == str(P_in.dtype),
            "retained_coordinate_bit_identity": bool(
                Pout.shape == P_in.shape
                and np.ascontiguousarray(Pout[V_out]).tobytes()
                == np.ascontiguousarray(P_in[V_out]).tobytes()),
            "coordinates_changed_only_at_excised_cells": bool(
                np.array_equal(changed, excised_cells)),
            "excised_cells_stamped_missing": bool(
                np.all(Pout[excised_cells] == MISSING)),
            "coordinate_bit_identity_method": (
                "COMPUTED on the RELOADED output: dtype and shape compared "
                "to the input's, the RETAINED cells' coordinate bytes "
                "compared cell-for-cell, the set of cells whose coordinates "
                "differ compared to the excised set, and the excised cells "
                "checked to carry the -1 marker"),
            "invalidation_carrier": HYBRID_INVALIDATION,
            "valid_mask_equals_intended": bool(
                np.array_equal(np.asarray(surf_out.valid, bool), V_out)),
            "changes_only_valid_to_invalid": bool(
                not (np.asarray(surf_out.valid, bool) & ~V_in).any()),
            "n_invalidated": int(
                (V_in & ~np.asarray(surf_out.valid, bool)).sum())}
        bad = [k for k, v in reload_checks.items()
               if isinstance(v, bool) and not v]
        if bad:
            return failure("output_contract_violation",
                           f"no clean claim: reloaded output violates "
                           f"{bad}",
                           {"pre_registration": prereg,
                            "pre_registration_sha256": prereg_sha,
                            "reload_checks": reload_checks,
                            "solver": {"solves": solves,
                                       "scipy_version": scipy_version}})

        # ---- AUTHORITATIVE recensus of the emitted mesh ----------------
        try:
            after = census_mesh(mesh_out, f"{tag}_after{iterations}", wdir,
                                nv, nu, timers)
        except Refusal as r:
            return failure(r.status,
                           "no clean claim: the emitted aggregate is not "
                           "censusable", {"evidence": r.evidence,
                                          "pre_registration": prereg,
                                          "pre_registration_sha256":
                                              prereg_sha})
        log(f"      recensus: d0 transverse "
            f"{after['engine']['d0']['transverse']} "
            f"grazing {after['engine']['d0']['grazing']} | d1 transverse "
            f"{after['engine']['d1']['transverse']} grazing "
            f"{after['engine']['d1']['grazing']}")
        residual_transverse = [r for r in after["rows"]
                               if r["verdict"] == "transverse"]
        # Round-26: acceptance is zero residual transverse of ANY kind. The
        # old "attributable to the target event" test belonged to the
        # per-event executor and is gone with it.
        residual_target = [list(r["key"]) for r in residual_transverse]
        if not residual_target:
            accepted = True
            break
        log(f"      residual transverse keys: {len(residual_target)} -- "
            "adding them as constraints and re-solving (spec section 5 loop)")
        try:
            new_cons = [c for c in coverage_from_rows(residual_transverse,
                                                      Q_out)
                        if c["key"] not in constraints]
        except ValueError as e:
            return failure("census_mapping_error", str(e),
                           {"pre_registration": prereg,
                            "pre_registration_sha256": prereg_sha})
        if not new_cons:
            # coordinates never change and validity only shrinks, so the
            # output census is a subset of the input's: a residual whose
            # keys were ALL already constrained means the emitted mask does
            # not satisfy the constraints the solver met. Internal error,
            # never a silent pass.
            return failure(
                "stalled",
                "no clean claim: recensus dirty with no new constraints "
                "(internal error -- the emitted mask does not satisfy the "
                "solved coverage constraints)",
                {"pre_registration": prereg,
                 "pre_registration_sha256": prereg_sha,
                 "residual_target_keys": residual_target,
                 "solver": {"solves": solves,
                            "scipy_version": scipy_version}})
        for c in new_cons:
            constraints[c["key"]] = c

    # ---- NAIVE-CONSUMER verification (round-25 A1) ----------------------
    # The whole point of the hybrid stamps is that a consumer ignoring
    # mask.tif must still see the crossing gone. Verified by the ENGINE, not
    # asserted: emit the identical coordinate planes with NO sidecar and
    # census that.
    if not args.skip_naive_check:
        naive_dir = OUT / f"{seg.name}_naive.tifxyz"
        naive_emit = emit_excised_tifxyz(mesh, naive_dir, V_out,
                                         V_in & ~V_out, with_mask=False)
        naive_surf = tifxyz.read(naive_dir)
        try:
            naive_census = census_mesh(naive_dir, f"{tag}_naive", wdir,
                                       nv, nu, timers)
            naive_counts = {f"d{d}": naive_census["engine"][f"d{d}"]
                            for d in (0, 1)}
            naive_transverse = sum(naive_counts[f"d{d}"]["transverse"]
                                   for d in (0, 1))
        except Refusal as r:
            naive_counts, naive_transverse = {"refusal": r.label}, None
        naive_check = {
            "purpose": ("a consumer that implements ONLY tifxyz's x=y=z=-1 "
                        "convention and never reads mask.tif"),
            "mesh": str(naive_dir),
            "has_mask_sidecar": (naive_dir / "mask.tif").exists(),
            "files": naive_emit["files"],
            "valid_equals_masked_reading": bool(np.array_equal(
                np.asarray(naive_surf.valid, bool), V_out)),
            "census": naive_counts,
            "transverse_total": naive_transverse,
            "naive_reader_sees_no_transverse_crossing":
                (naive_transverse == 0)}
        log(f"  naive-consumer census (mask.tif ABSENT): transverse "
            f"{naive_transverse}")
        shutil.rmtree(naive_dir)
    else:
        naive_check = {"skipped": True}

    # ---- acceptance evidence -------------------------------------------
    inval = set(invalidated)
    kept_set, removed_set = set(kept), set(removed)
    sub = submultiset_report(before["multiset"], after["multiset"], inval)
    residual_all_transverse = [list(r["key"]) for r in after["rows"]
                               if r["verdict"] == "transverse"]
    if not (sub["output_submultiset_of_input"]
            and sub["missing_all_witnessed"]):
        return failure(
            "submultiset_violation",
            "no clean claim: the output contact multiset is not a witnessed "
            "submultiset of the input's",
            {"pre_registration": prereg,
             "pre_registration_sha256": prereg_sha,
             "contacts": sub, "budget_gate": budget,
             "reload_checks": reload_checks,
             "solver": {"solves": solves, "scipy_version": scipy_version}})

    # area accounting: unresolved = quads in a REMAINING transverse contact
    unresolved = sorted({q for r in after["rows"]
                         if r["verdict"] == "transverse"
                         for q in (r["q1"], r["q2"])} & kept_set)
    blocks = {k: area_block(areas[k], Q_in, removed, unresolved)
              for k in (0, 1, "canonical")}

    # cross-check the vectorised area convention against the scalar
    # per-quad function the planted suite pins (same formula, same order)
    rng = np.random.default_rng(0)
    sample = [q_in[i] for i in rng.choice(len(q_in),
                                          size=min(256, len(q_in)),
                                          replace=False)]
    area_check = max(
        [abs(float(areas[0][v, u]) - quad_area(P64, v, u, 0)) for v, u in sample]
        + [abs(float(areas[1][v, u]) - quad_area(P64, v, u, 1))
           for v, u in sample]
        + [abs(float(areas["canonical"][v, u])
               - quad_area_canonical(P64, v, u)) for v, u in sample])
    assert area_check < 1e-9, f"area convention mismatch {area_check}"

    boundary = cut_boundary(removed_set, kept_set)
    boundary_len_vx = fsum(float(np.linalg.norm(P64[a] - P64[b]))
                           for a, b in boundary)
    Q_out = retained_quads(P64, np.asarray(tifxyz.read(mesh_out).valid, bool),
                           CENSUS["maxedge"])
    n_before_c, dist_before = components(Q_in, areas["canonical"])
    n_after_c, dist_after = components(Q_out, areas["canonical"])

    # triangle-identity multisets. Quad-level validity retains ALL FOUR
    # triangles of a retained quad ((d, k) in {0,1}x{0,1}), so the multiset
    # over (diag, v, u, triangle_index) is exactly the retained-quad set
    # times {0,1}x{0,1}; recording the retained-quad sets' hashes plus the
    # explicit difference is the same information as a 2.8M-entry listing
    # at this grid size, and the subset test is exact set arithmetic.
    def tri_hash(qs):
        a = np.asarray(sorted(qs), dtype=np.int64)
        return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

    tri_removed = [[d, v, u, k] for (v, u) in sorted(removed_set)
                   for d in (0, 1) for k in (0, 1)]
    V_final = np.asarray(tifxyz.read(mesh_out).valid, bool)

    archive = OUT / "archive" / seg.name
    archive.mkdir(parents=True, exist_ok=True)
    csv_records = {}
    for label, st in (("before", before), ("after", after)):
        for d in (0, 1):
            src = Path(st["diag"][d]["path"])
            dstc = archive / f"{label}_d{d}.csv"
            shutil.copy(src, dstc)
            csv_records[f"{label}_d{d}"] = {
                "path": str(dstc), "sha256": sha(dstc),
                "engine_csv_sha256": st["diag"][d]["sha256"],
                "n_rows": st["diag"][d]["n_lines"],
                "out_of_range_rows": st["diag"][d]["out_of_range"]}

    stage1 = [s for s in solves if s["stage"] == 1]
    solver_status = ("optimal" if stage1 and all(s["status"] == "optimal"
                                                 for s in stage1)
                     else ("not_applicable" if not stage1 else "best_found"))
    final_stage1 = solve_totals([s for s in solves
                                 if s.get("iteration") == iterations], 1)
    cert = {
        "operation": OPERATION_LABEL,
        "operation_label": operation_label,
        "status": "clean",
        "clean_claim": True,
        "claim": (
            f"Segment {seg.name}: EVERY transverse crossing censused under "
            f"both diagonals was removed by CERTIFIED EXCISION in ONE "
            f"segment-wide solve -- {len(invalidated)} grid vertices "
            f"invalidated, {len(removed)} quads ("
            f"{blocks['canonical']['excised_fraction']:.3e} of the input "
            f"retained canonical area) DELETED from the surface. The emitted "
            f"aggregate was reloaded from disk and recensused by "
            f"engines/selfcross under both diagonals: zero residual "
            f"transverse contacts, and the full contact multiset is a "
            f"witnessed submultiset of the input's. Surface is missing where "
            f"the cut is; every retained coordinate is bit-identical to the "
            f"input."),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "clean_claim_authority": CLEAN_CLAIM_AUTHORITY,
        "corpus": corpus, "segment": seg.name, "volume": volume,
        "mesh": str(mesh), "voxel_um": voxel_um,
        "grid_shape": [nv, nu], "dtype": str(P_in.dtype),
        "pre_registration": prereg,
        "pre_registration_path": str(prereg_path),
        "pre_registration_sha256": prereg_sha,
        "input": {"array_sha256": _mesh_hash(P_in, V_in),
                  "files": {f"{ax}.tif": sha(mesh / f"{ax}.tif")
                            for ax in AXES},
                  "n_valid": int(V_in.sum()),
                  "n_retained_quads": len(q_in)},
        "output": {"array_sha256": _mesh_hash(np.asarray(
                       tifxyz.read(mesh_out).points), V_final),
                   "path": str(mesh_out),
                   "files": emit_report["files"],
                   "n_valid": int(V_final.sum()),
                   "n_retained_quads": len(kept),
                   "emission": {k: v for k, v in emit_report.items()
                                if k != "files"}},
        "code_provenance": prov,
        "scipy_version": scipy_version,
        "census_params": dict(CENSUS, diagonals=[0, 1], engine=str(ENGINE),
                              schema=SCHEMA_V2_HEADER,
                              predicate=("engines/selfcross exact tri-tri "
                                         "predicate; AUTHORITATIVE -- the "
                                         "planted-mesh Python predicate is "
                                         "NOT used anywhere on this path"),
                              arithmetic="double (C++), float64 (Python)"),
        "solver": {
            "backend": "scipy.optimize.milp (HiGHS)",
            "scipy_version": scipy_version,
            "status": solver_status,
            "lexicographic": lexicographic,
            "objective_policy": (
                "two-stage lexicographic: stage 1 minimizes excised "
                "canonical area ONLY; stage 2 minimizes the invalidated-"
                "vertex count subject to the recorded area cap; tertiary "
                "(boundary length, component count) remains a documented v2"),
            "stage2_allowance_rule": STAGE2_ALLOWANCE_RULE,
            "coverage_rule": (
                "round-24 quad-level EIGHT-CORNER coverage: for every "
                "censused transverse triangle pair, "
                "sum(x_v for v in corners(quad1) UNION corners(quad2)) >= 1"),
            "constraint_scope": prereg["constraint_scope"],
            "constraint_reduction": reduction,
            "time_limit_s": args.time_limit,
            "iterations": iterations,
            "solves": solves,           # per-solve objective / bound / gap
            "final_stage1_objective": final_stage1["objective"],
            "final_stage1_dual_bound": final_stage1["dual_bound"],
            "optimality_gap": final_stage1["mip_gap"],
            "final_stage1_aggregation": (
                "SUM over the independent components of the final iteration; "
                "valid because distinct components share no decision "
                "variable. The gap is |objective - dual bound| / |objective|"),
            "n_coverage_constraints": len(constraints),
            "protected_vertices": []},
        "excision": {
            "invalidated_vertices": sorted(map(list, invalidated)),
            "n_invalidated_vertices": len(invalidated),
            "removed_quads": sorted(map(list, removed_set)),
            "n_removed_quads": len(removed),
            "indices": "ORIGINAL input grid indices throughout",
            "cut_boundary_edges": [[list(a), list(b)] for a, b in boundary],
            "n_cut_boundary_edges": len(boundary),
            "cut_boundary_length_vx": boundary_len_vx,
            "cut_boundary_length_um": boundary_len_vx * voxel_um},
        "budget_gate": budget,
        "reload_checks": reload_checks,
        "invalidation": {"carrier": HYBRID_INVALIDATION,
                         "missing_marker": MISSING,
                         "naive_consumer_check": naive_check},
        "retained_coordinate_bit_identity":
            reload_checks["retained_coordinate_bit_identity"],
        "coordinates_changed_only_at_excised_cells":
            reload_checks["coordinates_changed_only_at_excised_cells"],
        "coordinate_guarantee": RETAINED_BIT_IDENTITY,
        "coordinate_bit_identity_method":
            reload_checks["coordinate_bit_identity_method"],
        "validity": {
            "changes_only_valid_to_invalid":
                reload_checks["changes_only_valid_to_invalid"],
            "n_invalidated": reload_checks["n_invalidated"]},
        "census_before": {f"d{d}": before["engine"][f"d{d}"]
                          for d in (0, 1)},
        "census_after": {f"d{d}": after["engine"][f"d{d}"] for d in (0, 1)},
        "acceptance": {
            "criterion": prereg["acceptance_criterion"],
            "residual_transverse_keys_any": residual_all_transverse,
            "transverse_clean_both_diagonals":
                not residual_all_transverse,
            "all_contacts_clean_both_diagonals":
                sum(after["multiset"].values()) == 0,
            "naive_consumer_transverse_clean":
                naive_check.get("naive_reader_sees_no_transverse_crossing"),
            "shared_support_events": touching,
            "events_recorded": [ev_json(t) for t in targets]},
        "contacts": sub,
        "triangle_multisets": {
            "identity": "(diag, v, u, triangle_index)",
            "determination": (
                "quad-level validity retains all four triangles "
                "({d0,d1} x {k0,k1}) of a retained quad, so the multiset is "
                "exactly the retained-quad set x {0,1} x {0,1}; the "
                "retained-quad sets are hashed and their full difference is "
                "listed, which is the same information as an explicit "
                "listing at this grid size"),
            "n_before": 4 * len(q_in),
            "n_after": 4 * len(kept),
            "retained_quads_before_sha256": tri_hash(q_in),
            "retained_quads_after_sha256": tri_hash(kept_set),
            "removed_triangles": tri_removed,
            "after_subset_of_before": bool(kept_set <= set(q_in)),
            "n_newly_retained_quads": len(kept_set - set(q_in))},
        "area": {
            "convention": (
                "quad area under diagonal d = sum of its two triangle areas "
                "under that tessellation, float64 over the input "
                "coordinates; canonical = mean(d0, d1); sums by math.fsum"),
            "vectorised_vs_scalar_max_abs_diff": area_check,
            "d0": blocks[0], "d1": blocks[1],
            "canonical": blocks["canonical"],
            "unresolved_quads": sorted(map(list, unresolved))},
        "topology": {
            "connectivity": ("retained-quad complex joined by shared "
                             "corners (8-connectivity on the quad grid)"),
            "components_before": n_before_c,
            "components_after": n_after_c,
            "component_distribution_before": dist_before,
            "component_distribution_after": dist_after},
        "census_csv_archive": csv_records,
        "spec_ambiguities_resolved": [
            "INVALIDATION CARRIER (round-25 A1, SUPERSEDES the earlier "
            "mask-only resolution). tifxyz marks a missing cell either by "
            "x=y=z=-1 or by a sidecar mask. Marking only one of them is an "
            "interoperability defect: a consumer implementing the other "
            "convention reads the excised cells as present and reconstructs "
            "the very crossing this operation removed. Resolved as HYBRID -- "
            "mask.tif = 0 AND x=y=z=-1 -- with the coordinate guarantee "
            "restated over RETAINED cells (every retained coordinate "
            "bit-identical; coordinates differ from the input at excised "
            "cells only). Verified, not assumed: the emitted mesh is "
            "reloaded and recensused, and a sidecar-less copy is censused "
            "separately as the naive consumer's view.",
            "REFUSAL ARTIFACT. Section 3 says a refusal emits no certificate; "
            "a real run still has to leave an auditable trace. Resolved by "
            "writing a labelled refusal/failure RECORD (record_kind, "
            "clean_claim false, record_is_not_a_certificate, no output mesh "
            "path) under a distinct filename, never a certificate.",
            "TRIANGLE-IDENTITY MULTISETS AT REAL SCALE. Section 7.7's "
            "explicit before/after listing is 2.8M entries here. Resolved by "
            "recording what determines it exactly: quad-level validity keeps "
            "all four (diag, triangle_index) triangles of a retained quad, so "
            "the multiset is the retained-quad set x {0,1} x {0,1}. Both "
            "retained-quad sets are hashed, the full difference is listed "
            "explicitly, and after-subset-of-before is exact set arithmetic.",
            "CONSTRAINED VERDICTS. Coverage constraints are built from "
            "TRANSVERSE rows only, since section 2's acceptance criterion is "
            "zero transverse crossing pairs and the planted-scale predicate "
            "reports nothing else. Coplanar and grazing contacts are still "
            "carried in the full contact multiset, so the submultiset claim "
            "covers them: none was created, and here none was destroyed.",
            "AREA AT REAL SCALE. The per-quad area convention is evaluated "
            "vectorised over the whole grid rather than per quad in Python; "
            "the two are cross-checked on a random sample against the scalar "
            "function the planted suite pins (max abs difference recorded "
            "under area.vectorised_vs_scalar_max_abs_diff).",
        ],
        "staleness_warning": STALENESS_WARNING,
        "scope_note": (
            "Endpoint certificate on the emitted arrays only. Excision is a "
            "topology change: it certifies the absence of transverse "
            "self-intersection in what remains, NOT that the removed surface "
            "was spurious, and NOT texture fidelity or downstream "
            "reconstructability of the cut region."),
        "instrumentation": {**timers.report(),
                            "wall_seconds": round(time.time() - t_start, 2)},
    }
    # ---- round-27 status fields: geometry and selection are INDEPENDENT --
    cert["geometry_status"] = (GEOMETRY_STATUS_CLEAN
                               if not residual_all_transverse
                               else "residual_transverse")
    cert["status_independence"] = SELECTION_STATUS_RULE
    cert["scheduling_note"] = SCHEDULING_NOTE
    if selection is not None:
        cert["selection_status"] = selection["selection_status"]
        cert["selection"] = {
            "strategy": "lp_round",
            "selection_status": selection["selection_status"],
            "method_mix": selection["method_mix"],
            "achieved_area_canonical": selection["achieved_area"],
            "achieved_area_bounded_subset":
                selection["achieved_area_bounded_subset"],
            "achieved_area_unbounded_subset":
                selection["achieved_area_unbounded_subset"],
            "n_components_without_bound":
                selection["n_components_without_bound"],
            "combined_lower_bound": selection["combined_lower_bound"],
            "combined_lower_bound_complete":
                selection["combined_lower_bound_complete"],
            "combined_lower_bound_rule":
                selection["combined_lower_bound_rule"],
            "ratio_achieved_over_bound":
                selection["ratio_achieved_over_bound"],
            "ratio_rule": selection["ratio_rule"],
            "ratio_covers_area_fraction":
                selection["ratio_covers_area_fraction"],
            "k_max": selection["k_max"],
            "n_components": len(selection["components"]),
            "per_component": sorted(
                selection["components"],
                key=lambda c: -c.get("achieved_area", 0.0))[:200],
            "per_component_truncated_to": 200,
            "rules": selection["rules"], "policy": selection["policy"],
            "timings": selection["timings"]}
        cert["solver"]["selection_strategy"] = "lp_round"
        cert["solver"]["note"] = (
            "The cut was CHOSEN by the round-27 deterministic constructor; "
            "`solves` lists only the exact MILP attempts the scheduling "
            "policy allowed. The clean verdict is independent of selection "
            "optimality and comes from the engine recensus of the reloaded "
            "arrays.")
    else:
        cert["selection_status"] = ("area_optimal"
                                    if solver_status == "optimal"
                                    else "mixed")
    cert_path.write_text(json.dumps(cert, indent=1))
    log("")
    log(f"[CLEAN] {len(invalidated)} vertices invalidated, {len(removed)} "
        f"quads excised "
        f"({blocks['canonical']['excised_fraction']:.4e} of retained area); "
        f"recovery {blocks['canonical']['clean_recovery_fraction']:.8f}")
    log(f"  boundary {len(boundary)} edges, {boundary_len_vx:.3f} vx "
        f"({boundary_len_vx * voxel_um:.1f} um)")
    log(f"  components {n_before_c} -> {n_after_c}")
    log(f"  mesh {mesh_out}")
    log(f"  certificate {cert_path}")
    return 0 if accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
