"""SHADOW GLOBAL EXCISION: measure the cost of the August operator.

This is the W1 prove-or-kill instrument for round-26 Q2. It runs the real
segment-wide excision on a real segment and records what it cost --
WITHOUT needing the result to be good. A timeout, a solver failure, an
unusable area price and a clean cheap cut are all equally valid outputs; the
only failure mode this script has is not recording a number.

For one segment it:

  1. Loads the DISPLACEMENT-REPAIRED mesh from out/repaired/multi/meshes/ if
     one exists, else the original published mesh. (The corpus play is
     "displacement first, then one global residual excision", so the
     repaired mesh is the honest input.)
  2. Censuses it with engines/selfcross under BOTH diagonals -- authoritative.
  3. Builds ONE segment-wide MILP from EVERY transverse row (round 26): no
     event matching, no per-event iteration.
  4. Reduces the constraint set (dedup, dominance, independent components)
     and records the ratios.
  5. Solves, applies the mask, emits ONE aggregate tifxyz under HYBRID
     invalidation, reloads it from disk and recensuses it.
  6. Records: constraints before/after reduction, achieved area, solver
     status/objective/dual bound/gap, invalidated vertices, removed quads,
     retained fraction of canonical area per d0/d1/canonical, components
     before/after, per original component the retained fraction and the
     LARGEST-DESCENDANT RECOVERY R_main, descendant counts and area
     distributions, fragments >1%, cut-boundary length in edges and um,
     MILP time, total wall time, and the post-reload census transverse
     counts.

The record is snapshotted to disk after every stage, so a `timeout` kill
still leaves the numbers gathered up to that point -- with `status` naming
the stage that did not finish.

    uv run python bench/excise_shadow.py --segment 20231005123336

PRODUCTION CERTIFICATE MODE (round-28 Q3/Q4). The same pipeline, the same
frozen policy, the same thresholds -- only the provenance and the artifact
handling change:

    uv run python bench/excise_shadow.py --segment <name> --certificate \
        --base-manifest out/corpus_bases.json --out-root out/excised/corpus

  * `--base-manifest` pins the input mesh and the voxel scale from the
    hash-verified base manifest instead of rediscovering them, and the
    record carries that provenance (base kind/mesh/hashes, original mesh
    and hashes, originating repair-certificate hash, geometry key,
    duplicate_of, is_canonical, and the manifest file's own sha256).
  * `--certificate` names the record an EXCISION CERTIFICATE, writes it as
    <segment>_excision_certificate.json with the mesh as
    <segment>_excised.tifxyz, and ALWAYS keeps the emitted mesh, because
    the mesh is the production artifact.
  * Every exit records a `terminal_disposition`: there is no silent skip.

With none of those flags the script behaves exactly as before.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
import traceback
from math import fsum
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
sys.path.insert(0, "bench")
from windcheck import tifxyz                                    # noqa: E402
from windcheck.excise import (FROZEN_GREEDY_FIRST_RULE,          # noqa: E402
                              FROZEN_POLICY, FROZEN_POLICY_VERSION,
                              GEOMETRY_STATUS_CLEAN,
                              HYBRID_INVALIDATION, MISSING,
                              REDUCTION_RULE, SCHEDULING_NOTE,
                              SELECTION_STATUS_RULE, frozen_policy_hash,
                              reduce_constraints, select_global,
                              select_global_frozen, solve_global)
from windcheck.intrinsic import retained_quads                  # noqa: E402
from excise_segment import (AXES, CENSUS, Refusal, Timers,      # noqa: E402
                            census_mesh, code_provenance,
                            coverage_from_rows, emit_excised_tifxyz,
                            quad_area_grids, resolve_segment, solve_totals)
from repair_segment import RES_UM, sha                          # noqa: E402

OUT = Path("out/shadow")
REPAIRED = Path("out/repaired/multi/meshes")
FRAGMENT_FRACTION = 0.01        # a "fragment" is >1% of its parent's area
REPORT_COMPONENTS = 50          # per-component detail rows kept in the JSON

# ------------------------------------------------- production certificates
CERTIFICATE_RECORD_KIND = "excision certificate"
SHADOW_RECORD_KIND = "shadow measurement (NOT a certificate, NOT a claim)"

# Round-28 Q4(2): EVERY manifest entry ends on exactly one of these. There
# is no silent skip: a mesh nobody can census, and a mesh with no triangles
# to census, each get a named disposition carrying its evidence.
TERMINAL_DISPOSITIONS = (
    "transformed",              # cut, emitted, reloaded, recensused 0/0
    "already_clean",            # transverse-clean at input; NO cut made
    "duplicate_alias",          # byte-identical alias (written by the driver)
    "triangle_empty_invalid",   # no valid vertices / no retained quads
    "not_censusable",           # engine refusal, input or output
    "residual_transverse",      # cut and emitted, but contacts survive
    "error",                    # crash, timeout, or no feasible mask
)


# ------------------------------------------------------------- topology
CORE_AREA_FRACTION = 0.999      # round-28 Q1: the 99.9%-AREA CORE
CORE_R_MAIN_GATE = 0.90


def area_core(areas, rmains, core_fraction: float = CORE_AREA_FRACTION,
              gate: float = CORE_R_MAIN_GATE) -> dict:
    """Round-28 Q1 fragmentation gate, defined on the INPUT alone.

    Sort the segment's INPUT connected components by canonical area; the
    99.9%-AREA CORE is the SMALLEST PREFIX whose cumulative area reaches
    `core_fraction` of the input retained area. GATE: every core component
    must have R_main >= 0.90.

    The remainder is the SUB-0.1% COMPONENT TAIL. It is never called
    "debris": geometry has not established that those components are
    debris, only that they are small. The old "component >= 1% of segment
    area" rule had a loophole -- fifty genuine 0.9% components could hold
    45% of a segment and every one of them would escape the gate.
    """
    order = sorted(range(len(areas)), key=lambda i: (-areas[i], i))
    total = fsum(areas)
    core, cum = [], 0.0
    for i in order:
        if total > 0 and cum >= core_fraction * total:
            break
        core.append(i)
        cum += areas[i]
    tail = [i for i in order if i not in set(core)]
    core_r = [rmains[i] for i in core if rmains[i] is not None]
    return {
        "definition": ("the 99.9%-AREA CORE is the smallest prefix of the "
                       "INPUT components, sorted by canonical area, whose "
                       "cumulative area reaches 99.9% of the input retained "
                       "area; the remainder is the SUB-0.1% COMPONENT TAIL "
                       "(never 'debris')"),
        "core_fraction_target": core_fraction,
        "gate_R_main": gate,
        "n_input_components": len(areas),
        "n_core_components": len(core),
        "core_area": cum,
        "core_area_fraction": (cum / total) if total else None,
        "n_tail_components": len(tail),
        "tail_area": fsum(areas[i] for i in tail),
        "tail_area_fraction": ((fsum(areas[i] for i in tail) / total)
                               if total else None),
        "min_R_main_core": (min(core_r) if core_r else None),
        "n_core_components_below_gate": sum(1 for v in core_r if v < gate),
        "core_gate_pass": bool(core_r and min(core_r) >= gate) or not core_r,
        "core_component_indices": core[:REPORT_COMPONENTS],
    }


def component_recovery(Q_in: np.ndarray, Q_out: np.ndarray,
                       areas: np.ndarray) -> dict:
    """Round-26 Q3: area alone is inadequate, because a hairline cut that
    splits a sheet 50/50 scores 99.99% on area. The primary statistic is

        R_main(C) = area(largest OUTPUT descendant of C) / area(C)

    for every INPUT component C. Output quads are a subset of input quads
    and connectivity only shrinks, so every output component lies inside
    exactly one input component and the parent map is well defined.
    """
    from scipy import ndimage
    st = np.ones((3, 3), bool)          # shared corners == 8-connectivity
    lab_in, n_in = ndimage.label(Q_in, structure=st)
    lab_out, n_out = ndimage.label(Q_out, structure=st)
    a_in = np.bincount(lab_in[Q_in], weights=areas[Q_in],
                       minlength=n_in + 1)
    q_in = np.bincount(lab_in[Q_in], minlength=n_in + 1)
    if n_out:
        li, lo = lab_in[Q_out], lab_out[Q_out]
        a_out = np.bincount(lo, weights=areas[Q_out], minlength=n_out + 1)
        parent = np.zeros(n_out + 1, np.int64)
        parent[lo] = li                 # constant within an output component
    else:
        a_out = np.zeros(1)
        parent = np.zeros(1, np.int64)

    retained = np.zeros(n_in + 1)
    largest = np.zeros(n_in + 1)
    ndesc = np.zeros(n_in + 1, np.int64)
    frags = np.zeros(n_in + 1, np.int64)
    for lo_lab in range(1, n_out + 1):
        p = int(parent[lo_lab])
        retained[p] += a_out[lo_lab]
        ndesc[p] += 1
        if a_out[lo_lab] > largest[p]:
            largest[p] = a_out[lo_lab]
        if a_out[lo_lab] > FRAGMENT_FRACTION * a_in[p]:
            frags[p] += 1

    order = np.argsort(-a_in[1:]) + 1 if n_in else np.array([], np.int64)
    rows = []
    for i in order[:REPORT_COMPONENTS]:
        i = int(i)
        desc = sorted((float(a_out[j]) for j in range(1, n_out + 1)
                       if int(parent[j]) == i), reverse=True)
        rows.append({
            "input_component": i,
            "n_quads": int(q_in[i]),
            "area_canonical": float(a_in[i]),
            "area_fraction_of_segment": (float(a_in[i] / a_in[1:].sum())
                                         if n_in else None),
            "retained_fraction": (float(retained[i] / a_in[i])
                                  if a_in[i] else None),
            "R_main": (float(largest[i] / a_in[i]) if a_in[i] else None),
            "n_descendants": int(ndesc[i]),
            "descendant_areas_top10": desc[:10],
            "fragments_over_1pct": int(frags[i])})

    live = [(int(i), float(a_in[i])) for i in range(1, n_in + 1)
            if a_in[i] > 0]
    total = fsum(a for _, a in live)
    rmain = {i: (largest[i] / a_in[i]) for i, _ in live}
    big = [i for i, a in live if a >= FRAGMENT_FRACTION * total]
    core = area_core([a for _, a in live],
                     [float(rmain[i]) for i, _ in live])
    destroyed = [(i, a) for i, a in live if retained[i] == 0.0]
    desc_counts = [int(ndesc[i]) for i, _ in live]
    return {
        # ---- round-28 Q1: the INPUT-ONLY 99.9%-area core gate ------------
        "core_gate": core,
        "unthresholded": {
            "note": ("reported for EVERY input component, with no size "
                     "threshold of any kind"),
            "area_weighted_R_main": (fsum(rmain[i] * a for i, a in live)
                                     / total if total else None),
            "min_R_main_all_components": (min(rmain.values()) if rmain
                                          else None),
            "n_input_components_fully_destroyed": len(destroyed),
            "area_of_fully_destroyed_input_components":
                fsum(a for _, a in destroyed),
            "area_fraction_of_fully_destroyed_input_components":
                (fsum(a for _, a in destroyed) / total) if total else None,
            "n_descendants_total": int(sum(desc_counts)),
            "n_descendants_max": (max(desc_counts) if desc_counts else 0),
            "descendant_area_distribution_top20": sorted(
                (float(a_out[j]) for j in range(1, n_out + 1)),
                reverse=True)[:20],
            "n_output_descendants": int(n_out)},
        "definition": ("R_main(C) = area(largest output descendant of input "
                       "component C) / area(C); reported per INPUT component "
                       "because a fragmenting cut can score ~100% on area"),
        "components_before": int(n_in),
        "components_after": int(n_out),
        "min_R_main_any_component": (min(rmain.values()) if rmain else None),
        "min_R_main_components_over_1pct_of_area":
            (min(rmain[i] for i in big) if big else None),
        "n_components_over_1pct_of_area": len(big),
        "n_components_R_main_below_0.9": sum(1 for v in rmain.values()
                                             if v < 0.9),
        "n_components_over_1pct_R_main_below_0.9":
            sum(1 for i in big if rmain[i] < 0.9),
        "area_weighted_R_main": (fsum(rmain[i] * a for i, a in live) / total
                                 if total else None),
        "n_input_components_fully_destroyed": sum(1 for i, _ in live
                                                  if retained[i] == 0.0),
        "per_component": rows,
        "per_component_truncated_to": REPORT_COMPONENTS}


def cut_boundary_arrays(P64: np.ndarray, removed: np.ndarray,
                        kept: np.ndarray) -> tuple[int, float]:
    """Cut-boundary edges and their physical length, vectorised.

    A cut-boundary edge is a grid edge shared between a REMOVED quad and a
    RETAINED output quad. Quads (a,u) and (a+1,u) share the row edge
    ((a+1,u),(a+1,u+1)); quads (v,b) and (v,b+1) share the column edge
    ((v,b+1),(v+1,b+1)).
    """
    row_len = np.linalg.norm(P64[:, 1:] - P64[:, :-1], axis=-1)   # (nv, nu-1)
    col_len = np.linalg.norm(P64[1:, :] - P64[:-1, :], axis=-1)   # (nv-1, nu)
    m_row = ((removed[:-1, :] & kept[1:, :])
             | (kept[:-1, :] & removed[1:, :]))                   # at row a+1
    m_col = ((removed[:, :-1] & kept[:, 1:])
             | (kept[:, :-1] & removed[:, 1:]))                   # at col b+1
    nvq, nuq = removed.shape
    n = int(m_row.sum() + m_col.sum())
    length = (fsum(row_len[1:nvq, :nuq][m_row].tolist())
              + fsum(col_len[:nvq, 1:nuq][m_col].tolist()))
    return n, length


def headline_denominator(original: Path, removed: np.ndarray,
                         nv: int, nu: int, maxedge: float) -> dict:
    """Round-28 Q3 HEADLINE retention, the second area denominator.

    Operational retention is measured against the PRE-EXCISION BASE (which
    may be a displacement-repaired mesh). The HEADLINE number must instead
    price the REMOVED ORIGINAL QUAD INDICES using the ORIGINAL PUBLISHED
    COORDINATES, so that displacement-induced area changes cannot
    contaminate the preservation percentage in either direction.

    A removed quad that is not a retained quad of the ORIGINAL mesh
    contributes nothing: it was not part of the originally represented
    surface. Those are counted separately, never silently dropped.
    """
    surf = tifxyz.read(original)
    P0 = np.asarray(surf.points, np.float64)
    V0 = np.asarray(surf.valid, bool)
    if [int(s) for s in surf.shape] != [int(nv), int(nu)]:
        return {"measurable": False,
                "reason": ("the original published mesh has a different grid "
                           "shape, so quad indices are not comparable"),
                "original_grid_shape": [int(s) for s in surf.shape]}
    Q0 = retained_quads(P0, V0, maxedge)
    A0 = quad_area_grids(P0)["canonical"]
    a_total = fsum(A0[Q0].tolist())
    hit = removed & Q0
    a_removed = fsum(A0[hit].tolist())
    return {
        "measurable": True,
        "definition": ("HEADLINE retention = 1 - (canonical area of the "
                       "removed ORIGINAL quad indices, priced with the "
                       "ORIGINAL PUBLISHED coordinates) / (canonical area "
                       "of the ORIGINAL mesh's retained quads)"),
        "original_mesh": str(original),
        "A_original_canonical": a_total,
        "A_removed_priced_on_original": a_removed,
        "headline_retained_fraction": (1.0 - a_removed / a_total) if a_total
                                      else None,
        "headline_excised_fraction": (a_removed / a_total) if a_total
                                     else None,
        "n_removed_quads": int(removed.sum()),
        "n_removed_quads_in_original": int(hit.sum()),
        "n_removed_quads_absent_from_original": int((removed & ~Q0).sum()),
        "n_original_retained_quads": int(Q0.sum())}


# ------------------------------------------------------------------ main
def pick_input_mesh(seg_name: str, original: Path) -> dict:
    """The displacement-repaired mesh if one exists, else the original."""
    cand = REPAIRED / f"{seg_name}_repaired.tifxyz"
    if cand.is_dir() and (cand / "x.tif").exists():
        return {"path": cand, "source": "displacement_repaired",
                "note": ("out/repaired/multi/meshes: the corpus play is "
                         "bounded displacement first, then ONE global "
                         "residual excision")}
    return {"path": original, "source": "original_published",
            "note": "no displacement-repaired mesh exists for this segment"}


def manifest_entry(path: Path, segment: str) -> dict:
    """The round-28 Q3 BASE MANIFEST entry for one segment.

    In production the input mesh is NOT rediscovered by globbing: it is
    whatever the hash-verified manifest pinned, so a certificate can be
    replayed against the exact bytes it was cut from. The manifest file's
    own sha256 travels into the record for the same reason.
    """
    path = Path(path)
    doc = json.loads(path.read_text())
    for e in (doc.get("entries") or []):
        if e.get("segment") == segment:
            return {**e,
                    "base_manifest": str(path),
                    "base_manifest_schema": doc.get("schema"),
                    "base_manifest_sha256": sha(path)}
    raise SystemExit(f"{segment}: no entry in base manifest {path}")


def base_provenance(entry: dict) -> dict:
    """The subset of the manifest entry a certificate must carry."""
    keys = ("base_kind", "base_mesh", "base_hashes", "original_mesh",
            "original_hashes", "repair_certificate",
            "repair_certificate_sha256", "geometry_key",
            "original_geometry_key", "duplicate_of", "is_canonical",
            "base_manifest", "base_manifest_schema", "base_manifest_sha256")
    return {k: entry.get(k) for k in keys}


def mesh_plane_hashes(mesh: Path) -> dict:
    return {ax: sha(mesh / f"{ax}.tif") for ax in AXES}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--segment", required=True)
    ap.add_argument("--time-limit", type=float, default=1800.0,
                    help="TOTAL MILP budget in seconds, across components "
                         "(--strategy exact only)")
    ap.add_argument("--strategy", default="greedy_first",
                    choices=["greedy_first", "lp_round", "exact"],
                    help="round-28 FROZEN default `greedy_first`: an "
                         "area-aware greedy feasible incumbent for EVERY "
                         "component first, reverse-deleted immediately, then "
                         "OPTIONAL LP/exact improvement under segment-wide "
                         "budgets. Its policy constants are FROZEN and are "
                         "deliberately NOT exposed as flags. `lp_round` is "
                         "the round-27 prototype, kept for comparison")
    ap.add_argument("--lp-time-limit", type=float, default=600.0,
                    help="per-component LP cap; on expiry that component "
                         "falls back to the area-aware greedy")
    ap.add_argument("--lp-total-budget", type=float, default=2400.0,
                    help="segment-wide LP budget across all components")
    ap.add_argument("--exact-max-constraints", type=int, default=700)
    ap.add_argument("--exact-time-limit", type=float, default=30.0)
    ap.add_argument("--exact-total-budget", type=float, default=600.0,
                    help="SEGMENT-WIDE exact-MILP budget; components are "
                         "attempted cheapest first and skipped once it is "
                         "spent (costs optimality claims, never feasibility)")
    ap.add_argument("--stage2-max-constraints", type=int, default=50)
    ap.add_argument("--improve-budget", type=float, default=90.0,
                    help="SEGMENT-WIDE local-improvement budget in seconds")
    ap.add_argument("--label", default="excision",
                    help="operation label, e.g. junction_excision")
    ap.add_argument("--keep-mesh", action="store_true")
    ap.add_argument("--out-root", type=Path, default=OUT,
                    help="records, meshes and census workdirs go here "
                         f"(default {OUT}, the shadow location)")
    ap.add_argument("--certificate", action="store_true",
                    help="PRODUCTION CERTIFICATE MODE: record_kind becomes "
                         "'excision certificate', the record is written as "
                         "<segment>_excision_certificate.json, the mesh as "
                         "<segment>_excised.tifxyz, and the emitted mesh is "
                         "ALWAYS kept (it is the production artifact). No "
                         "measurement, threshold or policy changes.")
    ap.add_argument("--base-manifest", type=Path, default=None,
                    help="round-28 Q3 base manifest: take the input mesh and "
                         "the voxel scale from this segment's entry instead "
                         "of rediscovering them, and carry its hash-verified "
                         "provenance into the record")
    args = ap.parse_args(argv)

    t0 = time.time()
    timers = Timers()
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    base = None
    if args.base_manifest is not None:
        # The manifest pins EVERYTHING -- mesh, corpus, volume, voxel
        # scale. Nothing is rediscovered by globbing, which also makes
        # segments outside the pinned corpora roots (the expansion
        # manifest) runnable. Lookup is by exact segment name.
        base = manifest_entry(args.base_manifest, args.segment)
        seg_name = base["segment"]
        corpus = base.get("corpus")
        volume = base.get("volume")
        original = Path(base.get("original_mesh") or base["base_mesh"])
        voxel_um = (float(base["voxel_um"])
                    if base.get("voxel_um") is not None else None)
        chosen = {"path": Path(base["base_mesh"]),
                  "source": base.get("base_kind"),
                  "note": ("pinned by the round-28 base manifest "
                           f"{base['base_manifest']} "
                           f"(sha256 {base['base_manifest_sha256'][:12]}), "
                           "not rediscovered")}
    else:
        corpus, seg, volume = resolve_segment(args.segment)
        seg_name = seg.name
        original = sorted(seg.glob(f"mesh/*{volume}*.tifxyz"))[0]
        m = RES_UM.search(original.name)
        voxel_um = float(m.group(1)) if m else 7.91
        chosen = pick_input_mesh(seg_name, original)
    mesh = chosen["path"]
    kind = "certificate" if args.certificate else "shadow"
    keep_mesh = bool(args.keep_mesh or args.certificate)
    tag = hashlib.sha256((seg_name + kind).encode()).hexdigest()[:12]
    wdir = out_root / f"work_{tag}"
    wdir.mkdir(parents=True, exist_ok=True)
    rec_path = out_root / (f"{seg_name}_excision_certificate.json"
                           if args.certificate else f"{seg_name}_shadow.json")
    mesh_out = out_root / (f"{seg_name}_excised.tifxyz" if args.certificate
                           else f"{seg_name}_shadow.tifxyz")
    prov = code_provenance()

    R: dict = {
        "record_kind": (CERTIFICATE_RECORD_KIND if args.certificate
                        else SHADOW_RECORD_KIND),
        "purpose": ("round-26 Q2: measure what one segment-wide certified "
                    "excision costs, before promising a corpus"),
        "segment": seg_name, "corpus": corpus, "volume": volume,
        "operation_label": args.label,
        "voxel_um": voxel_um,
        "input_mesh": str(mesh), "input_mesh_source": chosen["source"],
        "input_mesh_note": chosen["note"],
        "original_mesh": str(original),
        "time_limit_s": args.time_limit,
        "strategy": args.strategy,
        "policy_version": (FROZEN_POLICY_VERSION
                           if args.strategy == "greedy_first" else "round27"),
        "policy_hash": (frozen_policy_hash()
                        if args.strategy == "greedy_first" else None),
        "frozen_policy": (FROZEN_POLICY if args.strategy == "greedy_first"
                          else None),
        "selection_policy": {
            "lp_time_limit_s": args.lp_time_limit,
            "lp_total_budget_s": args.lp_total_budget,
            "exact_max_constraints": args.exact_max_constraints,
            "exact_time_limit_s": args.exact_time_limit,
            "exact_total_budget_s": args.exact_total_budget,
            "stage2_max_constraints": args.stage2_max_constraints,
            "improve_budget_s": args.improve_budget,
            "scheduling_note": SCHEDULING_NOTE},
        "selection_status_rule": SELECTION_STATUS_RULE,
        "census_params": dict(CENSUS, diagonals=[0, 1]),
        "constraint_reduction_rule": REDUCTION_RULE,
        "invalidation": HYBRID_INVALIDATION,
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "code_provenance": prov,
        "source_tree_digest": prov["source_tree_digest"],
        "code_version": prov["code_version"],
        "terminal_disposition": None,
        "terminal_dispositions_defined": list(TERMINAL_DISPOSITIONS),
        "output_mesh": None,
        "output_mesh_hashes": None,
        "status": "started",
        "stage": "load",
    }
    if base is not None:
        R.update(base_provenance(base))

    def snap(**kw):
        """Persist after every stage: a `timeout` kill must still leave the
        numbers gathered so far, with the unfinished stage named."""
        R.update(kw)
        R["wall_seconds"] = round(time.time() - t0, 2)
        R["instrumentation"] = timers.report()
        rec_path.write_text(json.dumps(R, indent=1, default=str))

    def terminal(disposition, **kw):
        """Every exit of this script names its terminal disposition; there
        is no path that leaves the record without one (round-28 Q4(2))."""
        assert disposition in TERMINAL_DISPOSITIONS, disposition
        snap(terminal_disposition=disposition, **kw)

    def log(msg):
        print(msg, flush=True)

    log(f"{corpus} {seg_name}")
    log(f"  input {mesh}  ({chosen['source']})  voxel "
        f"{voxel_um if voxel_um is not None else 'unknown'} um")
    snap()

    try:
        surf = tifxyz.read(mesh)
        P_in = np.asarray(surf.points)
        P64 = np.asarray(P_in, np.float64)
        V_in = np.asarray(surf.valid, bool)
        nv, nu = surf.shape
        Q_in = retained_quads(P64, V_in, CENSUS["maxedge"])
        areas = quad_area_grids(P64)
        a_in_can = fsum(areas["canonical"][Q_in].tolist())
        snap(stage="census_before",
             grid_shape=[int(nv), int(nu)],
             n_valid=int(V_in.sum()),
             n_retained_quads=int(Q_in.sum()),
             input_area_canonical=a_in_can)
        log(f"  grid {nv}x{nu}  valid {int(V_in.sum()):,}  "
            f"retained quads {int(Q_in.sum()):,}")

        # ---- 0. triangle-empty / invalid input ---------------------------
        # Not a skip and not an error: an input with no valid vertices or no
        # retained quads has no triangles, so "transverse-clean" would be
        # vacuous. It gets its own terminal disposition, with the counts as
        # evidence (round-28 Q4(2)).
        if int(V_in.sum()) == 0 or int(Q_in.sum()) == 0:
            terminal("triangle_empty_invalid",
                     status="triangle_empty_invalid", stage="done",
                     evidence={"n_valid_vertices": int(V_in.sum()),
                               "n_retained_quads": int(Q_in.sum()),
                               "grid_shape": [int(nv), int(nu)],
                               "maxedge": CENSUS["maxedge"]},
                     note=("the input carries no triangles, so no census, no "
                           "cut and no cleanliness claim is defined on it"))
            log("  [TRIANGLE-EMPTY/INVALID] no triangles in the input")
            return 0

        # ---- 1. authoritative census -------------------------------------
        try:
            before = census_mesh(mesh, f"{tag}_before", wdir, nv, nu, timers)
        except Refusal as r:
            terminal("not_censusable", status="not_censusable",
                     stage="census_before", label=r.label, evidence=r.evidence)
            log(f"  [NOT CENSUSABLE] {r.label}")
            return 0
        cb = {f"d{d}": before["engine"][f"d{d}"] for d in (0, 1)}
        n_trans = sum(cb[f"d{d}"]["transverse"] for d in (0, 1))
        snap(stage="constraints", census_before=cb,
             input_transverse_total=n_trans)
        log(f"  census before: d0 transverse {cb['d0']['transverse']} | "
            f"d1 transverse {cb['d1']['transverse']}")

        if n_trans == 0:
            terminal("already_clean",
                     status="already_transverse_clean", stage="done",
                     note="no excision required; nothing to measure",
                     operational_retained_fraction=1.0,
                     headline_retained_fraction=1.0,
                     n_removed_quads=0,
                     core_gate_pass=True,
                     core_gate_note=("NO CUT WAS MADE, so every input "
                                     "component survives whole: R_main = 1 "
                                     "for the 99.9%-area core by "
                                     "construction, not by measurement"),
                     emptiness_guard={"clean_by_emptiness": False,
                                      "n_retained_quads": int(Q_in.sum()),
                                      "rule": ("nothing was removed, so the "
                                               "output cannot be empty")},
                     claimed_clean=True,
                     geometry_status=GEOMETRY_STATUS_CLEAN)
            log("  [ALREADY CLEAN] no transverse contacts to excise")
            return 0

        # ---- 2. segment-wide constraints ---------------------------------
        t = time.time()
        cons = coverage_from_rows(before["rows"], Q_in)
        timers.add("build_constraints", time.time() - t)
        t = time.time()
        red = reduce_constraints(cons, Q_in)
        timers.add("reduce_constraints", time.time() - t)
        red_stats = {k: v for k, v in red.items()
                     if k not in ("components", "kept", "rule")}
        red_stats["component_sizes"] = red_stats["component_sizes"][:100]
        red_stats["component_sizes_truncated_to"] = 100
        snap(stage="solve", constraint_reduction=red_stats,
             n_candidate_vertices=len({p for c in cons for p in c["coverage"]}))
        log(f"  constraints: {red['n_raw']} raw -> {red['n_after_dedup']} "
            f"deduped -> {red['n_after_dominance']} after dominance; "
            f"{red['n_components']} independent components "
            f"(largest {red['largest_component']})")

        # ---- 3. ONE segment-wide selection --------------------------------
        selection = None
        if args.strategy in ("greedy_first", "lp_round"):
            t = time.time()
            if args.strategy == "greedy_first":
                sel = select_global_frozen(cons, P64, Q_in, (),
                                           area_grid=areas["canonical"])
            else:
                sel = select_global(
                    cons, P64, Q_in, (),
                    lp_time_limit=args.lp_time_limit,
                    lp_total_budget=args.lp_total_budget,
                    exact_max_constraints=args.exact_max_constraints,
                    exact_time_limit=args.exact_time_limit,
                    exact_total_budget=args.exact_total_budget,
                    stage2_max_constraints=args.stage2_max_constraints,
                    improve_budget=args.improve_budget,
                    area_grid=areas["canonical"])
            sel_s = time.time() - t
            timers.add("selection", sel_s)
            comps = sel["components"]
            per_method: dict = {}
            for c in comps:
                b = per_method.setdefault(c["method"],
                                          {"n": 0, "constraints": 0,
                                           "achieved_area": 0.0,
                                           "lower_bound": 0.0,
                                           "bound_complete": True})
                b["n"] += 1
                b["constraints"] += c["n_constraints"]
                b["achieved_area"] += c.get("achieved_area", 0.0)
                if c["lower_bound"] is None:
                    b["bound_complete"] = False
                else:
                    b["lower_bound"] += c["lower_bound"]
            ratios = [c["ratio_achieved_over_bound"] for c in comps
                      if c.get("ratio_achieved_over_bound") is not None]
            phase = {k: round(v, 3) for k, v in sel["timings"].items()}
            selection = {
                "strategy": args.strategy,
                "policy_version": sel.get("policy_version"),
                "policy_hash": sel.get("policy_hash"),
                "scheduling_rule": sel.get("scheduling_rule"),
                "greedy_incumbent_area": sel.get("greedy_incumbent_area"),
                "greedy_construction_seconds":
                    sel.get("greedy_construction_seconds"),
                "improvement_over_greedy": sel.get("improvement_over_greedy"),
                "minimum_area_claim_admissible":
                    sel.get("minimum_area_claim_admissible"),
                "minimum_area_claim_rule": sel.get("minimum_area_claim_rule"),
                "n_components_lp_attempted":
                    sel.get("n_components_lp_attempted"),
                "n_components_lp_skipped": sel.get("n_components_lp_skipped"),
                "n_components_lp_improved":
                    sel.get("n_components_lp_improved"),
                "n_components_exact_attempted":
                    sel.get("n_components_exact_attempted"),
                "n_components_exact_skipped":
                    sel.get("n_components_exact_skipped"),
                "status": sel["status"],
                "selection_status": sel["selection_status"],
                "method_mix": sel["method_mix"],
                "per_method": per_method,
                "achieved_area_canonical": sel["achieved_area"],
                "achieved_area_bounded_subset":
                    sel["achieved_area_bounded_subset"],
                "achieved_area_unbounded_subset":
                    sel["achieved_area_unbounded_subset"],
                "n_components_without_bound":
                    sel["n_components_without_bound"],
                "combined_lower_bound": sel["combined_lower_bound"],
                "combined_lower_bound_complete":
                    sel["combined_lower_bound_complete"],
                "combined_lower_bound_rule": sel["combined_lower_bound_rule"],
                "ratio_achieved_over_bound": sel["ratio_achieved_over_bound"],
                "ratio_rule": sel["ratio_rule"],
                "ratio_covers_area_fraction":
                    sel["ratio_covers_area_fraction"],
                "ratio_min_component": (min(ratios) if ratios else None),
                "ratio_max_component": (max(ratios) if ratios else None),
                "k_max": sel["k_max"],
                "phase_timings_s": phase,
                "wall_seconds": sel_s,
                "n_components": len(comps),
                "components_sorted_by_area": sorted(
                    comps, key=lambda c: -c.get("achieved_area", 0.0)
                )[:REPORT_COMPONENTS],
                "components_truncated_to": REPORT_COMPONENTS,
                "n_exact_records": len(sel["milp_records"]),
                "rules": sel["rules"], "policy": sel["policy"],
                "scipy_version": sel["scipy_version"]}
            snap(stage="apply", selection=selection)
            log(f"  selection {sel['selection_status']} in {sel_s:.1f}s: "
                f"mix {sel['method_mix']} area {sel['achieved_area']:.4f} "
                f"bound {sel['combined_lower_bound']:.4f} ratio "
                f"{sel['ratio_achieved_over_bound']} "
                f"(phases {phase})")
            if sel["status"] != "ok":
                terminal("error", status=f"selection_{sel['status']}",
                         stage="done",
                         note=("no feasible mask: NO output, NO clean claim. "
                               "This is a recorded result, not an error."))
                log(f"  [SELECTION {sel['status'].upper()}] no feasible mask")
                return 0
            sol = {"chosen": sel["chosen"]}
        else:
            t = time.time()
            sol = solve_global(cons, P64, Q_in, set(), args.time_limit)
            milp_s = time.time() - t
            timers.add("milp", milp_s)
            agg = solve_totals(sol["records"], 1)
            solver = {"status": sol["status"],
                      "lexicographic": sol["lexicographic"],
                      "objective_achieved_area": agg["objective"],
                      "dual_bound": agg["dual_bound"],
                      "mip_gap": agg["mip_gap"],
                      "component_statuses": agg["statuses"],
                      "n_component_solves": agg["n_solves"],
                      "milp_seconds": milp_s,
                      "scipy_version": sol["scipy_version"],
                      "solves": sol["records"][:200],
                      "solves_truncated_to": 200,
                      "n_solve_records": len(sol["records"])}
            snap(stage="apply", solver=solver)
            log(f"  MILP {sol['status']} in {milp_s:.1f}s: area "
                f"{agg['objective']} bound {agg['dual_bound']} gap "
                f"{agg['mip_gap']}")
            if sol["status"] not in ("optimal", "best_found"):
                terminal("error", status=f"solver_{sol['status']}",
                         stage="done",
                         note=("no usable incumbent: NO output, NO clean "
                               "claim. A recorded result, not an error."))
                log(f"  [SOLVER {sol['status'].upper()}] no usable incumbent")
                return 0

        # ---- 4. apply the mask -------------------------------------------
        invalidated = sorted(sol["chosen"])
        X = np.zeros((nv, nu), bool)
        for v, u in invalidated:
            X[v, u] = True
        removed = (X[:-1, :-1] | X[1:, :-1] | X[:-1, 1:] | X[1:, 1:]) & Q_in
        kept = Q_in & ~removed
        V_out = V_in & ~X
        excised_cells = V_in & X
        Q_out = retained_quads(P64, V_out, CENSUS["maxedge"])
        assert np.array_equal(Q_out, kept), \
            "masked retained-quad set disagrees with the mask"

        area = {}
        for k, name in ((0, "d0"), (1, "d1"), ("canonical", "canonical")):
            A = areas[k]
            ai = fsum(A[Q_in].tolist())
            ac = fsum(A[kept].tolist())
            ae = fsum(A[removed].tolist())
            area[name] = {"A_input": ai, "A_retained": ac, "A_excised": ae,
                          "retained_fraction": (ac / ai) if ai else None,
                          "excised_fraction": (ae / ai) if ai else None,
                          "identity_residual": ac + ae - ai}
        area["denominator"] = (
            "OPERATIONAL: every fraction above is relative to the "
            "PRE-EXCISION BASE (the input mesh actually cut, which may be "
            "displacement-repaired). See `headline_area` for the second, "
            "original-coordinate denominator.")
        t = time.time()
        headline = (headline_denominator(original, removed, nv, nu,
                                         CENSUS["maxedge"])
                    if Path(original) != Path(mesh) else
                    dict(measurable=True,
                         definition=("the input IS the original published "
                                     "mesh, so the headline and operational "
                                     "denominators coincide by construction"),
                         original_mesh=str(original),
                         A_original_canonical=area["canonical"]["A_input"],
                         A_removed_priced_on_original=
                             area["canonical"]["A_excised"],
                         headline_retained_fraction=
                             area["canonical"]["retained_fraction"],
                         headline_excised_fraction=
                             area["canonical"]["excised_fraction"],
                         n_removed_quads=int(removed.sum()),
                         n_removed_quads_in_original=int(removed.sum()),
                         n_removed_quads_absent_from_original=0,
                         n_original_retained_quads=int(Q_in.sum()),
                         input_is_original=True))
        timers.add("headline_area", time.time() - t)
        n_edges, blen_vx = cut_boundary_arrays(P64, removed, kept)
        rec_comp = component_recovery(Q_in, Q_out, areas["canonical"])
        excision = {
            "n_invalidated_vertices": len(invalidated),
            "n_removed_quads": int(removed.sum()),
            "n_retained_quads": int(kept.sum()),
            "cut_boundary_edges": n_edges,
            "cut_boundary_length_vx": blen_vx,
            "cut_boundary_length_um": (blen_vx * voxel_um
                                       if voxel_um is not None else None)}
        snap(stage="emit", area=area, headline_area=headline,
             excision=excision, component_recovery=rec_comp)
        log(f"  headline retained "
            f"{headline.get('headline_retained_fraction')}"
            f"  (core gate pass "
            f"{rec_comp['core_gate']['core_gate_pass']}, core comps "
            f"{rec_comp['core_gate']['n_core_components']}, min R_main core "
            f"{rec_comp['core_gate']['min_R_main_core']})")
        log(f"  cut: {len(invalidated)} vertices, {int(removed.sum()):,} "
            f"quads; retained canonical "
            f"{area['canonical']['retained_fraction']:.6f}; components "
            f"{rec_comp['components_before']} -> "
            f"{rec_comp['components_after']}; min R_main "
            f"{rec_comp['min_R_main_any_component']}")

        # ---- 5. emit, reload, RECENSUS -----------------------------------
        emit = timers.timed("emit", emit_excised_tifxyz, mesh, mesh_out,
                            V_out, excised_cells)
        surf_out = tifxyz.read(mesh_out)
        Pout = np.asarray(surf_out.points)
        changed = np.any(Pout != P_in, axis=-1)
        reload_checks = {
            "grid_shape_equal": list(surf_out.shape) == [nv, nu],
            "dtype_equal": str(surf_out.points.dtype) == str(P_in.dtype),
            "retained_coordinate_bit_identity": bool(
                np.ascontiguousarray(Pout[V_out]).tobytes()
                == np.ascontiguousarray(P_in[V_out]).tobytes()),
            "coordinates_changed_only_at_excised_cells": bool(
                np.array_equal(changed, excised_cells)),
            "excised_cells_stamped_missing": bool(
                np.all(Pout[excised_cells] == MISSING)),
            "valid_mask_equals_intended": bool(np.array_equal(
                np.asarray(surf_out.valid, bool), V_out)),
            "changes_only_valid_to_invalid": bool(
                not (np.asarray(surf_out.valid, bool) & ~V_in).any())}
        snap(stage="census_after", emission={k: v for k, v in emit.items()
                                             if k != "files"},
             reload_checks=reload_checks,
             output_mesh=str(mesh_out),
             output_mesh_hashes=mesh_plane_hashes(mesh_out))

        try:
            after = census_mesh(mesh_out, f"{tag}_after", wdir, nv, nu,
                                timers)
            ca = {f"d{d}": after["engine"][f"d{d}"] for d in (0, 1)}
            after_trans = sum(ca[f"d{d}"]["transverse"] for d in (0, 1))
        except Refusal as r:
            terminal("not_censusable", status="output_not_censusable",
                     stage="done", label=r.label, evidence=r.evidence)
            log(f"  [OUTPUT NOT CENSUSABLE] {r.label}")
            return 0
        log(f"  recensus: d0 transverse {ca['d0']['transverse']} | "
            f"d1 transverse {ca['d1']['transverse']}")

        emptiness = {
            "retained_fraction_canonical":
                area["canonical"]["retained_fraction"],
            "n_retained_quads": int(kept.sum()),
            "clean_by_emptiness": bool(
                kept.sum() == 0
                or (area["canonical"]["retained_fraction"] or 0) < 0.5),
            "rule": ("a segment that satisfies the headline by emitting an "
                     "empty or near-empty surface does NOT count; flagged "
                     "when nothing is retained or under half the canonical "
                     "area survives")}
        ok = (after_trans == 0
              and all(v for k, v in reload_checks.items()
                      if isinstance(v, bool))
              and not emptiness["clean_by_emptiness"])
        geometry_status = (GEOMETRY_STATUS_CLEAN if after_trans == 0
                           else "residual_transverse")
        terminal("transformed" if after_trans == 0
                 else "residual_transverse",
                 status=("transverse_clean" if after_trans == 0
                         else "residual_transverse"),
                 stage="done", census_after=ca,
                 output_transverse_total=after_trans,
                 geometry_status=geometry_status,
                 selection_status=(selection or {}).get("selection_status"),
                 status_independence=SELECTION_STATUS_RULE,
                 emptiness_guard=emptiness,
                 core_gate_pass=bool(rec_comp["core_gate"]["core_gate_pass"]),
                 operational_retained_fraction=(
                     area["canonical"]["retained_fraction"]),
                 headline_retained_fraction=(
                     headline.get("headline_retained_fraction")),
                 claimed_clean=bool(ok))
        log(f"  [{R['status'].upper()}] retained "
            f"{area['canonical']['retained_fraction']:.6f} of canonical "
            f"area; wall {time.time() - t0:.0f}s")
    except Exception as exc:                     # a crash is also a RESULT
        terminal("error", status="error", error=repr(exc),
                 traceback=traceback.format_exc())
        log(f"  [ERROR] {exc!r}")
        return 1
    finally:
        # --certificate emits a PRODUCTION ARTIFACT: the mesh is always kept.
        if not keep_mesh and mesh_out.exists():
            shutil.rmtree(mesh_out)
            R["output_mesh"] = None
        R["output_mesh_kept"] = bool(keep_mesh)
        R["wall_seconds"] = round(time.time() - t0, 2)
        R["instrumentation"] = timers.report()
        rec_path.write_text(json.dumps(R, indent=1, default=str))
        print(f"  record {rec_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
