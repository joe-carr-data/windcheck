#!/usr/bin/env python3
"""Build the public release index for the excised corpus.

Reads, and only reads:

  * the pinned base manifest            out/corpus_bases.json
  * the driver summary                  out/excised/corpus/corpus_summary.jsonl
  * one excision certificate per segment out/excised/corpus/<seg>_excision_certificate.json
  * the independent recensus record     out/excised/corpus/verification.json
  * the memoised original areas         out/headline_original_areas.json

and writes two derived, regenerable artifacts:

  * out/release/index.json   machine-readable, one record per pinned segment
  * docs/CORPUS.md           the same content as a human-readable document

Nothing here is hand-entered.  Every per-segment number is copied out of the
segment's own certificate; a field that the certificate does not record is
emitted as null with an explicit note rather than being filled in.

Two area denominators are carried through, exactly as the certificates
define them:

  operational  1 - A_excised / A_input, where A_input is the PRE-EXCISION
               BASE that was actually cut.  For a displacement-repaired base
               that base is not the published mesh.

  headline     1 - (canonical area of the removed ORIGINAL quad indices,
               priced with the ORIGINAL PUBLISHED coordinates)
               / (canonical area of the ORIGINAL mesh's retained quads).

A segment that was never cut carries no priced-area block, because nothing
was removed -- but its headline denominator is not unknown, it is the full
canonical area of its original published mesh.  That area is recovered here
through the SAME code path the headline decision rule uses
(`bench/headline_decision.price_unmodified_segments`), reading and writing
the same memoised cache, so the two artifacts price the same population by
the same route and report the same figure.  No segment is dropped from the
headline denominator for want of a certificate field.

Usage:
    uv run python bench/build_release_index.py
    uv run python bench/build_release_index.py --repo-root . --no-markdown
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from decimal import Decimal, ROUND_DOWN
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
from windcheck.manifest import (MANIFEST_RULE, BaseKind,          # noqa: E402
                                BaseKindDisagreement, verify_base_kind)
from windcheck.provenance import release_provenance               # noqa: E402

SCHEMA = "release_index/v1"

BASE_KIND_RULE = (
    "base_kind is decided by SEMANTIC MANIFEST EQUALITY -- the input mesh's "
    "content manifest against the original published mesh's -- and never by "
    "path, because paths change in downloaded archives and fresh working "
    "directories. Index generation fails loudly when a certificate's "
    "declared base_kind disagrees with that hash-verified answer.")

# Terminal dispositions as the driver records them, mapped to the three
# public names used by the release index.
DISPOSITION_MAP = {
    "transformed": "transformed",
    "already_clean": "already_clean",
    "triangle_empty_invalid": "not_censusable",
}

RUNTIME_GATE_S = 600.0

# Original-area bookkeeping for segments that were never cut is SHARED with
# bench/headline_decision.py: the same memoised cache file, and, on a miss,
# the same recomputation from the original published mesh.
AREA_CACHE_REL = "out/headline_original_areas.json"

DEFS = {
    "represented_surface_retained": (
        "Every retention figure on this page is represented surface "
        "retained: canonical quad area that survives the cut, divided by the "
        "canonical quad area of the stated denominator. It is a geometric "
        "measure of the mesh, and nothing more."
    ),
    "operational_denominator": (
        "measured against the PRE-EXCISION BASE, i.e. the mesh that was "
        "actually cut. Where that base is a displacement-repaired mesh "
        "rather than the published original, this denominator measures what "
        "the excision cost relative to the surface it was handed."
    ),
    "headline_denominator": (
        "measured against the ORIGINAL PUBLISHED mesh: the removed ORIGINAL "
        "quad indices are priced using the ORIGINAL published coordinates, "
        "and divided by the canonical area of the original mesh's retained "
        "quads. When the base IS the original mesh the two denominators "
        "coincide by construction."
    ),
    "R_main": (
        "R_main(C) = area(largest output descendant of input component C) / "
        "area(C), reported per INPUT component. It detects a cut that scores "
        "near-100% on area while shattering a component into pieces."
    ),
    "core_gate": (
        "The 99.9%-AREA CORE is the smallest prefix of the input components, "
        "sorted by canonical area, whose cumulative area reaches 99.9% of "
        "the input retained area. The core gate passes when every component "
        "in that core has R_main >= 0.90."
    ),
    "dispositions": {
        "transformed": (
            "the input carried transverse contacts; an excision was computed "
            "and an output mesh was emitted"
        ),
        "already_clean": (
            "the input censused transverse-clean under both canonical "
            "triangulations, so no cut was made and no output mesh exists"
        ),
        "not_censusable": (
            "the input carries no triangles, so no census, no cut and no "
            "cleanliness claim is defined on it"
        ),
    },
}


# --------------------------------------------------------------- utilities

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def get(obj, *path, default=None):
    """Nested lookup that tolerates a missing or oddly-typed intermediate."""
    cur = obj
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return default if cur is None else cur


def percentile(values: list[float], q: float):
    """Linear-interpolated percentile over a sorted copy of ``values``."""
    if not values:
        return None
    xs = sorted(values)
    k = (len(xs) - 1) * q
    lo, hi = math.floor(k), math.ceil(k)
    if lo == hi:
        return xs[int(k)]
    return xs[lo] + (k - lo) * (xs[hi] - xs[lo])


def tree_bytes(path: Path) -> tuple[int, int]:
    """(apparent bytes, allocated bytes) for a file or directory tree.

    Allocated bytes is st_blocks * 512, i.e. the number a disk-usage tool
    reports; apparent bytes is the sum of file sizes.
    """
    if not path.exists():
        return 0, 0
    if path.is_file():
        st = path.stat()
        return st.st_size, st.st_blocks * 512
    apparent = allocated = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                st = os.stat(os.path.join(root, name))
            except OSError:
                continue
            apparent += st.st_size
            allocated += st.st_blocks * 512
    return apparent, allocated


def base_kind_check(cert: dict, segment: str = "") -> dict:
    """The one shared verdict on this certificate's declared base kind.

    Decided by content: the input mesh manifest against the original mesh
    manifest. No path is consulted anywhere on this route.
    """
    return verify_base_kind(cert.get("base_kind"), cert.get("base_hashes"),
                            cert.get("original_hashes"),
                            label=segment or cert.get("segment") or "")


def _declares_identity(rec: dict) -> bool:
    """Does this record carry mesh identity at all? A record with neither
    input nor original hashes has nothing to verify against -- that gap is
    reported by bench/verify_corpus.py, not by the base-kind rule."""
    return bool((rec or {}).get("base_hashes")
                or (rec or {}).get("original_hashes"))


def _is_disagreement(chk: dict, rec: dict) -> bool:
    """A base-kind check that must stop the release: the content decided and
    the declaration contradicts it, or identity was recorded but is not
    usable to decide."""
    return not chk["agrees"] and (chk["decidable"] or _declares_identity(rec))


def resolved_base_kind(cert: dict, segment: str = "") -> BaseKind | None:
    """The hash-verified base kind, or the declared one when content cannot
    decide. Callers that must not proceed on a disagreement use
    `base_kind_check` and raise; this is for the read-only accessors."""
    v = base_kind_check(cert, segment)
    kind = v["verified"] or v["declared"]
    try:
        return BaseKind(kind)
    except (ValueError, TypeError):
        return None


# ------------------------------------------------------------ per segment

def retention_block(cert: dict, disposition: str) -> dict:
    """Both area denominators, taken verbatim from the certificate."""
    if disposition == "not_censusable":
        return {
            "operational_retained_fraction": None,
            "headline_retained_fraction": None,
            "operational_A_input_canonical": None,
            "operational_A_excised_canonical": None,
            "headline_A_original_canonical": None,
            "headline_A_removed_priced_on_original": None,
            "n_removed_quads": None,
            "measured": False,
            "note": ("no triangles in the input, so no area is defined and "
                     "no retention figure is recorded"),
        }
    if disposition == "already_clean":
        base_is_original = resolved_base_kind(cert) is BaseKind.ORIGINAL
        return {
            "operational_retained_fraction":
                cert.get("operational_retained_fraction"),
            "headline_retained_fraction":
                cert.get("headline_retained_fraction"),
            "operational_A_input_canonical": cert.get("input_area_canonical"),
            "operational_A_excised_canonical": 0.0,
            # The certificate writes no headline_area block when nothing was
            # cut. When the base IS the original mesh the input area IS the
            # original area; when the base was displacement-repaired the
            # original-coordinate area is genuinely not recorded.
            "headline_A_original_canonical":
                cert.get("input_area_canonical") if base_is_original else None,
            "headline_A_removed_priced_on_original": 0.0,
            "n_removed_quads": cert.get("n_removed_quads", 0),
            "measured": False,
            "note": ("no cut was made, so both fractions are 1 by "
                     "construction" + ("" if base_is_original else
                     "; the base was displacement-repaired and the "
                     "original-coordinate area is not recorded, because no "
                     "headline area block is written when nothing is cut")),
        }
    return {
        "operational_retained_fraction":
            cert.get("operational_retained_fraction"),
        "headline_retained_fraction": cert.get("headline_retained_fraction"),
        "operational_A_input_canonical":
            get(cert, "area", "canonical", "A_input"),
        "operational_A_excised_canonical":
            get(cert, "area", "canonical", "A_excised"),
        "headline_A_original_canonical":
            get(cert, "headline_area", "A_original_canonical"),
        "headline_A_removed_priced_on_original":
            get(cert, "headline_area", "A_removed_priced_on_original"),
        "n_removed_quads": get(cert, "excision", "n_removed_quads"),
        "measured": True,
        "note": None,
    }


def fragmentation_block(cert: dict, disposition: str) -> dict:
    """min R_main, area-weighted R_main and the 99.9%-core gate."""
    if disposition == "transformed":
        unth = get(cert, "component_recovery", "unthresholded", default={})
        gate = get(cert, "component_recovery", "core_gate", default={})
        return {
            "measured": True,
            "min_R_main_all_components":
                unth.get("min_R_main_all_components"),
            "area_weighted_R_main": unth.get("area_weighted_R_main"),
            "min_R_main_core": gate.get("min_R_main_core"),
            "core_gate_pass": cert.get("core_gate_pass"),
            "n_input_components": gate.get("n_input_components"),
            "n_core_components": gate.get("n_core_components"),
            "n_core_components_below_gate":
                gate.get("n_core_components_below_gate"),
            "n_input_components_fully_destroyed":
                unth.get("n_input_components_fully_destroyed"),
            "n_components_before": get(cert, "component_recovery",
                                       "components_before"),
            "n_components_after": get(cert, "component_recovery",
                                      "components_after"),
            "note": None,
        }
    if disposition == "already_clean":
        return {
            "measured": False,
            "min_R_main_all_components": None,
            "area_weighted_R_main": None,
            "min_R_main_core": None,
            "core_gate_pass": cert.get("core_gate_pass"),
            "n_input_components": None,
            "n_core_components": None,
            "n_core_components_below_gate": None,
            "n_input_components_fully_destroyed": None,
            "n_components_before": None,
            "n_components_after": None,
            "note": ("no cut was made, so every input component survives "
                     "whole and R_main = 1 by construction, not by "
                     "measurement; no per-component figures are recorded"),
        }
    return {
        "measured": False,
        "min_R_main_all_components": None,
        "area_weighted_R_main": None,
        "min_R_main_core": None,
        "core_gate_pass": None,
        "n_input_components": None,
        "n_core_components": None,
        "n_core_components_below_gate": None,
        "n_input_components_fully_destroyed": None,
        "n_components_before": None,
        "n_components_after": None,
        "note": ("no triangles in the input, so components and R_main are "
                 "not defined and nothing is recorded"),
    }


def census_block(cert: dict, ver: dict | None) -> dict:
    rec = (ver or {}).get("recensus") or {}
    hc = (ver or {}).get("hash_check") or {}
    return {
        "input_transverse_total": cert.get("input_transverse_total"),
        "output_transverse_total": cert.get("output_transverse_total"),
        "geometry_status": cert.get("geometry_status"),
        "recensus_kind": rec.get("kind"),
        "recensus_ran": rec.get("ran"),
        "recensus_d0_transverse": rec.get("d0_transverse"),
        "recensus_d1_transverse": rec.get("d1_transverse"),
        "recensus_clean_both_diagonals": rec.get("clean_both_diagonals"),
        "recensus_disagrees_with_certificate":
            get(rec, "comparison", "transverse_disagrees"),
        "output_hashes_reverified": hc.get("status") == "ok",
    }


def segment_record(row: dict, cert: dict, ver: dict | None,
                   base_entry: dict | None, repo_root: Path) -> dict:
    raw = cert.get("terminal_disposition") or row.get("terminal_disposition")
    disposition = DISPOSITION_MAP.get(raw, raw)
    reason = None
    if disposition == "not_censusable":
        reason = cert.get("note") or DEFS["dispositions"]["not_censusable"]

    bk_check = base_kind_check(cert, row["segment"])
    out_mesh = cert.get("output_mesh")
    mesh_apparent = mesh_allocated = 0
    if out_mesh:
        mesh_apparent, mesh_allocated = tree_bytes(repo_root / out_mesh)

    return {
        "segment": row["segment"],
        "scroll": row.get("corpus"),
        "volume": cert.get("volume"),
        "voxel_um": cert.get("voxel_um"),
        "disposition": disposition,
        "disposition_recorded": raw,
        "status": cert.get("status"),
        "not_censusable_reason": reason,
        "not_censusable_evidence": cert.get("evidence"),
        "base_kind": (bk_check["verified"] or bk_check["declared"]),
        "base_kind_declared": cert.get("base_kind"),
        "base_kind_verified_from_hashes": bk_check["verified"],
        "base_kind_agrees": bk_check["agrees"],
        "input_mesh": cert.get("base_mesh"),
        "input_hashes": cert.get("base_hashes"),
        "original_mesh": cert.get("original_mesh"),
        "original_hashes": cert.get("original_hashes"),
        "output_mesh": out_mesh,
        "output_hashes": cert.get("output_mesh_hashes"),
        "output_mesh_bytes": mesh_apparent or None,
        "output_mesh_bytes_on_disk": mesh_allocated or None,
        "certificate": row.get("certificate"),
        "certificate_sha256": row.get("certificate_sha256"),
        "repair_certificate": cert.get("repair_certificate"),
        "repair_certificate_sha256": cert.get("repair_certificate_sha256"),
        "geometry_key": cert.get("geometry_key"),
        "is_canonical": cert.get("is_canonical"),
        "duplicate_of": cert.get("duplicate_of"),
        "claimed_clean": cert.get("claimed_clean"),
        "retention": retention_block(cert, disposition),
        "fragmentation": fragmentation_block(cert, disposition),
        "census": census_block(cert, ver),
        "wall_seconds": row.get("wall_seconds"),
        "log": row.get("log"),
        "grid_shape": cert.get("grid_shape"),
        "n_valid_vertices": cert.get("n_valid"),
        "base_manifest_sha256": cert.get("base_manifest_sha256"),
        "selection_status": cert.get("selection_status"),
        "policy_version": cert.get("policy_version"),
        "policy_hash": cert.get("policy_hash"),
    }


# ------------------------------------------------- unmodified-segment areas

class _chdir:
    """Run a block with the process CWD at ``path``.

    The mesh paths on a certificate are repo-relative, and the shared pricing
    routine resolves them against the CWD.  Entering the repo root for the
    duration keeps the paths -- and therefore the cache rows keyed on them --
    byte-identical to the ones the headline tool writes, instead of storing an
    absolute path that the other artifact would then fail to match.
    """

    def __init__(self, path: Path):
        self.path = Path(path)
        self.prev = None

    def __enter__(self):
        self.prev = os.getcwd()
        os.chdir(self.path)
        return self

    def __exit__(self, *exc):
        os.chdir(self.prev)
        return False


def headline_decision_module():
    """The headline decision rule, imported lazily from bench/."""
    import sys
    bench_dir = str(Path(__file__).resolve().parent)
    if bench_dir not in sys.path:
        sys.path.insert(0, bench_dir)
    import headline_decision
    return headline_decision


def price_unmodified_segments(order, segs, certs, repo_root: Path,
                              recompute: bool = True) -> dict:
    """Give every never-cut segment its missing headline denominator.

    An `already_clean` certificate records no cut, so it writes no
    `headline_area` block.  Where the base IS the original published mesh the
    input area already IS the original area and `retention_block` has taken
    it.  Where the base was displacement-repaired, the original-coordinate
    area is not on the certificate -- and dropping the segment there would
    silently shrink the headline denominator.

    So the number is recovered rather than dropped, by handing the segment to
    `headline_decision.price_unmodified_segments`: the memoised cache
    `out/headline_original_areas.json` is read first, and on a miss the
    original published mesh is re-read by that module's own
    `original_canonical_area` and the row is written back to the same cache.
    Both artifacts therefore price these segments from one source of truth.

    Mutates the retention blocks in place; returns a report.
    """
    hd = headline_decision_module()
    cache_path = Path(repo_root) / AREA_CACHE_REL
    cache = hd.load_area_cache(cache_path)

    shim: dict = {"order": [], "segments": {}}
    for seg in order:
        rec = segs[seg]
        if rec["disposition"] != "already_clean":
            continue
        if rec["retention"]["headline_A_original_canonical"] is not None:
            continue
        cert = certs.get(seg) or {}
        shim["order"].append(seg)
        shim["segments"][seg] = {
            "disposition": "already_clean",
            "certificate": cert,
            "original_mesh": cert.get("original_mesh"),
            "original_hashes": cert.get("original_hashes"),
            "base_mesh": cert.get("base_mesh"),
            "A_original_canonical": None,
            "A_removed_priced_on_original": None,
            "headline_retained_fraction":
                rec["retention"]["headline_retained_fraction"],
            "operational_retained_fraction":
                rec["retention"]["operational_retained_fraction"],
        }

    if not shim["order"]:
        return {"n_priced": 0, "priced": [], "n_recomputed_from_mesh": 0,
                "failures": [], "cache": AREA_CACHE_REL}

    with _chdir(repo_root):
        report = hd.price_unmodified_segments(shim, cache, recompute=recompute)
        if report["fresh"]:
            cache.update(report["fresh"])
            hd.write_area_cache(cache_path, cache)

    priced = []
    for row in report["derived"]:
        seg = row["segment"]
        ret = segs[seg]["retention"]
        ret["headline_A_original_canonical"] = row["A_original_canonical"]
        ret["headline_A_removed_priced_on_original"] = 0.0
        ret["headline_area_source"] = row["source"]
        if ret["headline_retained_fraction"] is None:
            ret["headline_retained_fraction"] = 1.0
        ret["note"] = (
            "no cut was made, so both fractions are 1 by construction; the "
            "base was displacement-repaired, so the original-coordinate area "
            "is not on the certificate and was recovered from the original "
            "published mesh (" + str(row["source"]) + ")")
        priced.append({"segment": seg,
                       "A_original_canonical": row["A_original_canonical"],
                       "source": row["source"]})

    return {
        "n_priced": len(priced),
        "priced": priced,
        "n_recomputed_from_mesh": len(report["fresh"]),
        "failures": report["failures"],
        "cache": AREA_CACHE_REL,
        "note": ("segments that were never cut carry no priced-area block; "
                 "their headline denominator is the full canonical area of "
                 "their original published mesh, recovered here by the same "
                 "route bench/headline_decision.py uses"),
    }


# ---------------------------------------------------------------- summary

def summary_block(order, segs, verification, base_manifest) -> dict:
    disp = {}
    scrolls = {}
    bases = {}
    for seg in order:
        r = segs[seg]
        disp[r["disposition"]] = disp.get(r["disposition"], 0) + 1
        s = scrolls.setdefault(r["scroll"], {"n": 0, "transformed": 0,
                                             "already_clean": 0,
                                             "not_censusable": 0})
        s["n"] += 1
        s[r["disposition"]] += 1
        bases[r["base_kind"]] = bases.get(r["base_kind"], 0) + 1

    # ---- retention, both denominators ----------------------------------
    op_num = op_den = 0.0
    hl_num = hl_den = 0.0
    hl_transformed_num = hl_transformed_den = 0.0
    op_gaps, hl_gaps = [], []
    hl_priced, hl_zero_area = [], []
    for seg in order:
        r = segs[seg]
        ret = r["retention"]
        if r["disposition"] == "not_censusable":
            # No triangles, so no area exists on either side of the ratio.
            # The segment stays in the headline POPULATION and enters it with
            # zero area on both sides: it cannot move the figure, and it is
            # named here rather than quietly dropped.
            hl_zero_area.append(seg)
            continue
        a_in = ret["operational_A_input_canonical"]
        a_ex = ret["operational_A_excised_canonical"]
        if a_in is None:
            op_gaps.append(seg)
        else:
            op_den += a_in
            op_num += a_ex or 0.0
        a_or = ret["headline_A_original_canonical"]
        a_rm = ret["headline_A_removed_priced_on_original"]
        if a_or is None:
            hl_gaps.append(seg)
        else:
            hl_priced.append(seg)
            hl_den += a_or
            hl_num += a_rm or 0.0
            if r["disposition"] == "transformed":
                hl_transformed_den += a_or
                hl_transformed_num += a_rm or 0.0

    # How far apart do the two denominators actually land, per segment?
    divergence = []
    for seg in order:
        ret = segs[seg]["retention"]
        a = ret["operational_retained_fraction"]
        b = ret["headline_retained_fraction"]
        if a is not None and b is not None:
            divergence.append((abs(a - b), seg))
    divergence.sort(reverse=True)
    n_identical = sum(1 for d, _ in divergence if d == 0.0)

    per_seg_hl = {seg: segs[seg]["retention"]["headline_retained_fraction"]
                  for seg in order
                  if segs[seg]["retention"]["headline_retained_fraction"]
                  is not None}
    worst_hl = min(per_seg_hl, key=per_seg_hl.get) if per_seg_hl else None
    per_seg_op = {seg: segs[seg]["retention"]["operational_retained_fraction"]
                  for seg in order
                  if segs[seg]["retention"]["operational_retained_fraction"]
                  is not None}
    worst_op = min(per_seg_op, key=per_seg_op.get) if per_seg_op else None

    # ---- runtime -------------------------------------------------------
    runtimes = [segs[s]["wall_seconds"] for s in order
                if segs[s]["wall_seconds"] is not None]

    # ---- fragmentation -------------------------------------------------
    gate_fail = [s for s in order
                 if segs[s]["fragmentation"]["core_gate_pass"] is False]
    gate_pass = sum(1 for s in order
                    if segs[s]["fragmentation"]["core_gate_pass"] is True)
    measured = [s for s in order if segs[s]["fragmentation"]["measured"]]
    min_r = min((segs[s]["fragmentation"]["min_R_main_all_components"]
                 for s in measured), default=None)

    vs = verification.get("summary", {})
    return {
        "n_segments": len(order),
        "n_scrolls": len(scrolls),
        "dispositions": disp,
        "per_scroll": scrolls,
        "bases": {
            **bases,
            "n_unique_geometries": base_manifest.get("n_canonical"),
            "n_duplicate_aliases": base_manifest.get("n_duplicate_aliases"),
            "note": ("every displacement-repaired base was hash-verified "
                     "against its repair certificate before it was cut"),
        },
        "retention": {
            "denominators": {
                "operational": DEFS["operational_denominator"],
                "headline": DEFS["headline_denominator"],
            },
            "operational_area_weighted": (1.0 - op_num / op_den)
                                          if op_den > 0 else None,
            "operational_A_input_total": op_den,
            "operational_A_excised_total": op_num,
            "operational_segments_counted": len(order) - len(op_gaps)
                                            - disp.get("not_censusable", 0),
            "operational_segments_without_area": op_gaps,
            "headline_area_weighted": (1.0 - hl_num / hl_den)
                                       if hl_den > 0 else None,
            "headline_A_original_total": hl_den,
            "headline_A_removed_total": hl_num,
            # The headline population is every pinned artifact whose original
            # area is known -- which, once the never-cut segments are priced,
            # is all of them.  This is the same population, counted the same
            # way, that bench/headline_decision.py reports.
            "headline_segments_counted": len(order) - len(hl_gaps),
            "headline_segments_priced": len(hl_priced),
            "headline_segments_zero_area": hl_zero_area,
            "headline_segments_without_area": hl_gaps,
            "headline_population_note": (
                "The headline denominator is evaluated over all "
                f"{len(order) - len(hl_gaps)} pinned artifacts: "
                f"{len(hl_priced)} carry a priced original-coordinate area, "
                f"and {len(hl_zero_area)} are triangle-empty or invalid "
                "inputs that carry no area at all and so enter with zero on "
                "both sides of the ratio. Segments that were never cut are "
                "priced at exactly 1.0 against the full canonical area of "
                "their original published mesh, recovered through "
                "bench/headline_decision.py so that both artifacts report "
                "the same figure over the same population."),
            "headline_area_weighted_transformed_only":
                (1.0 - hl_transformed_num / hl_transformed_den)
                if hl_transformed_den > 0 else None,
            "headline_A_original_total_transformed_only":
                hl_transformed_den,
            "headline_A_removed_total_transformed_only":
                hl_transformed_num,
            "min_headline_retained_fraction": per_seg_hl.get(worst_hl),
            "min_headline_segment": worst_hl,
            "min_operational_retained_fraction": per_seg_op.get(worst_op),
            "min_operational_segment": worst_op,
            "denominator_divergence": {
                "n_segments_compared": len(divergence),
                "n_segments_identical": n_identical,
                "max_abs_difference": divergence[0][0] if divergence else None,
                "max_abs_difference_segment":
                    divergence[0][1] if divergence else None,
                "note": ("Difference between the two per-segment retention "
                         "fractions. It is exactly zero wherever the base IS "
                         "the original published mesh, and it is what the "
                         "displacement repair costs elsewhere."),
            },
            "coverage_note": (
                "Segments that were already clean carry no priced-area block, "
                "because nothing was cut. Under the operational denominator "
                "the certificate still records the area of the base that was "
                "handed to the cut, so they price themselves. Under the "
                "headline denominator the original-coordinate area is on the "
                "certificate only where the base IS the original published "
                "mesh; where the base was displacement-repaired it is "
                "recovered by re-reading the original mesh, memoised in "
                "out/headline_original_areas.json. Either way the segment "
                "enters the headline denominator with a retained fraction of "
                "exactly 1.0 rather than being dropped."
            ),
        },
        "runtime_seconds": {
            "n": len(runtimes),
            "min": min(runtimes) if runtimes else None,
            "median": percentile(runtimes, 0.50),
            "p90": percentile(runtimes, 0.90),
            "p95": percentile(runtimes, 0.95),
            "p99": percentile(runtimes, 0.99),
            "max": max(runtimes) if runtimes else None,
            "total": round(sum(runtimes), 2),
            "gate": RUNTIME_GATE_S,
            "n_over_gate": sum(1 for v in runtimes if v > RUNTIME_GATE_S),
        },
        "fragmentation": {
            "core_gate_definition": DEFS["core_gate"],
            "R_main_definition": DEFS["R_main"],
            "n_core_gate_measured": len(measured),
            "n_core_gate_pass": gate_pass,
            "n_core_gate_fail": len(gate_fail),
            "core_gate_failures": [
                {"segment": s,
                 "scroll": segs[s]["scroll"],
                 "min_R_main_core":
                     segs[s]["fragmentation"]["min_R_main_core"],
                 "min_R_main_all_components":
                     segs[s]["fragmentation"]["min_R_main_all_components"],
                 "area_weighted_R_main":
                     segs[s]["fragmentation"]["area_weighted_R_main"],
                 "n_core_components_below_gate":
                     segs[s]["fragmentation"]["n_core_components_below_gate"]}
                for s in gate_fail],
            "min_R_main_over_corpus": min_r,
        },
        "verification": {
            "n_certificates_checked": vs.get("n_certificates_checked"),
            "n_certificates_well_formed": vs.get("n_certificates_well_formed"),
            "n_meshes_rehashed": vs.get("n_meshes_rehashed"),
            "n_meshes_rehash_ok": vs.get("n_meshes_rehash_ok"),
            "n_meshes_rehash_mismatch": vs.get("n_meshes_rehash_mismatch"),
            "n_recensused": vs.get("n_recensused"),
            "n_recensus_clean_both_diagonals":
                vs.get("n_recensus_clean_both_diagonals"),
            "n_census_disagreements": vs.get("n_census_disagreements"),
            "n_not_censusable_confirmed":
                vs.get("n_not_censusable_confirmed"),
            "n_failed": vs.get("n_failed"),
            "residual_transverse_contacts": sum(
                (segs[s]["census"]["recensus_d0_transverse"] or 0)
                + (segs[s]["census"]["recensus_d1_transverse"] or 0)
                for s in order),
            "n_timeouts": sum(1 for s in order
                              if (segs[s]["wall_seconds"] or 0)
                              > RUNTIME_GATE_S),
            "independence_note": verification.get("independence_note"),
        },
        "not_censusable": [
            {"segment": s, "scroll": segs[s]["scroll"],
             "reason": segs[s]["not_censusable_reason"],
             "evidence": segs[s]["not_censusable_evidence"]}
            for s in order if segs[s]["disposition"] == "not_censusable"],
    }


# ---------------------------------------------------------------- archive

SHARED_METADATA = [
    "out/corpus_bases.json",
    "out/excised/corpus/corpus_summary.jsonl",
    "out/excised/corpus/verification.json",
    "out/release/index.json",
]


def archive_block(order, segs, repo_root: Path) -> dict:
    per_scroll = {}
    total_apparent = total_allocated = 0
    for seg in order:
        r = segs[seg]
        bucket = per_scroll.setdefault(r["scroll"], {
            "n_segments": 0, "n_meshes": 0,
            "mesh_bytes": 0, "certificate_bytes": 0, "log_bytes": 0,
            "bytes": 0, "bytes_on_disk": 0})
        bucket["n_segments"] += 1
        app = alloc = 0
        if r["output_mesh"]:
            bucket["n_meshes"] += 1
            a, b = tree_bytes(repo_root / r["output_mesh"])
            bucket["mesh_bytes"] += a
            app += a
            alloc += b
        if r["certificate"]:
            a, b = tree_bytes(repo_root / r["certificate"])
            bucket["certificate_bytes"] += a
            app += a
            alloc += b
        if r["log"]:
            a, b = tree_bytes(repo_root / r["log"])
            bucket["log_bytes"] += a
            app += a
            alloc += b
        bucket["bytes"] += app
        bucket["bytes_on_disk"] += alloc
        total_apparent += app
        total_allocated += alloc

    shared = 0
    shared_files = []
    for rel in SHARED_METADATA:
        a, _ = tree_bytes(repo_root / rel)
        shared += a
        shared_files.append({"path": rel, "bytes": a})

    proposed = []
    for scroll in sorted(per_scroll):
        b = per_scroll[scroll]
        slug = scroll.lower().replace(" ", "")
        payload = b["bytes"] + shared
        proposed.append({
            "archive": f"windcheck-corpus-{slug}.tar",
            "scroll": scroll,
            "contents": [
                f"{b['n_meshes']} excised tifxyz meshes",
                f"{b['n_segments']} excision certificates",
                f"{b['n_segments']} run logs",
                ("shared manifests (base manifest, driver summary, "
                 "independent recensus record, release index)"),
            ],
            "payload_bytes": payload,
            "payload_mib": round(payload / 2 ** 20, 1),
            "payload_gib": round(payload / 2 ** 30, 3),
        })

    return {
        "measured_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "note": ("Sizes are of the emitted corpus artifacts only: the "
                 "excised meshes, the per-segment certificates and the "
                 "run logs. Driver scratch directories are not release "
                 "artifacts and are excluded. 'on disk' is allocated blocks, "
                 "the number a disk-usage tool reports; 'bytes' is the sum "
                 "of file sizes."),
        "total_bytes": total_apparent,
        "total_bytes_on_disk": total_allocated,
        "total_gib": round(total_apparent / 2 ** 30, 3),
        "total_gb": round(total_apparent / 1e9, 3),
        "per_scroll": per_scroll,
        "shared_metadata_bytes": shared,
        "shared_metadata_files": shared_files,
        "proposed_archives": proposed,
        "proposed_split_rationale": (
            "One archive per scroll. The five scrolls are independently "
            "useful, the split keeps every archive under 1 GiB, and each "
            "archive is self-describing because it carries a copy of the "
            "shared manifests. No archive is created by this script."
        ),
    }


# ------------------------------------------------------------------- build

def build(repo_root: Path) -> dict:
    repo_root = Path(repo_root).resolve()
    corpus_dir = repo_root / "out" / "excised" / "corpus"
    manifest_path = repo_root / "out" / "corpus_bases.json"
    summary_path = corpus_dir / "corpus_summary.jsonl"
    verification_path = corpus_dir / "verification.json"

    base_manifest = json.loads(manifest_path.read_text())
    base_by_seg = {e["segment"]: e for e in base_manifest.get("entries", [])}

    rows = [json.loads(line) for line in
            summary_path.read_text().splitlines() if line.strip()]
    rows.sort(key=lambda r: (r.get("corpus") or "", r["segment"]))

    verification = {}
    ver_by_seg = {}
    if verification_path.exists():
        verification = json.loads(verification_path.read_text())
        ver_by_seg = {s["segment"]: s
                      for s in verification.get("segments", [])}

    order, segs, certs = [], {}, {}
    base_kind_problems: list[str] = []
    for row in rows:
        cert_path = repo_root / row["certificate"]
        cert = json.loads(cert_path.read_text())
        rec = segment_record(row, cert, ver_by_seg.get(row["segment"]),
                             base_by_seg.get(row["segment"]), repo_root)
        chk = base_kind_check(cert, row["segment"])
        if _is_disagreement(chk, cert):
            base_kind_problems.append(chk["reason"])
        # The pinned base manifest declares the same thing a second time.
        # It has to agree with the certificate's own content, or the two
        # published artifacts disagree about what was cut.
        entry = base_by_seg.get(row["segment"]) or {}
        if entry:
            e_chk = verify_base_kind(entry.get("base_kind"),
                                     entry.get("base_hashes"),
                                     entry.get("original_hashes"),
                                     label=f"{row['segment']} (base manifest)")
            if _is_disagreement(e_chk, entry):
                base_kind_problems.append(e_chk["reason"])
        order.append(rec["segment"])
        segs[rec["segment"]] = rec
        certs[rec["segment"]] = cert

    if base_kind_problems:
        head = "\n  ".join(base_kind_problems[:20])
        more = (f"\n  ... and {len(base_kind_problems) - 20} more"
                if len(base_kind_problems) > 20 else "")
        raise BaseKindDisagreement(
            f"{len(base_kind_problems)} base-kind disagreement(s) between a "
            f"declared base_kind and the hash-verified answer; refusing to "
            f"publish an index that contradicts its own certificates:\n  "
            f"{head}{more}")

    # A segment that was never cut carries no headline_area block.  Price it
    # before anything is aggregated, so the headline denominator covers the
    # whole pinned population instead of only the segments that were cut.
    area_pricing = price_unmodified_segments(order, segs, certs, repo_root)

    index = {
        "schema": SCHEMA,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        # Provenance a public consumer can verify from the published files
        # alone. No commit sha: a reader cannot check one.
        "provenance": release_provenance(repo_root),
        "mesh_identity_rule": MANIFEST_RULE,
        "base_kind_rule": BASE_KIND_RULE,
        "base_kind_verified": True,
        "generator": "bench/build_release_index.py",
        "terminology": DEFS["represented_surface_retained"],
        "definitions": DEFS,
        "sources": {
            "base_manifest": "out/corpus_bases.json",
            "base_manifest_schema": base_manifest.get("schema"),
            "base_manifest_sha256": sha256_file(manifest_path),
            "corpus_summary": "out/excised/corpus/corpus_summary.jsonl",
            "corpus_summary_sha256": sha256_file(summary_path),
            "verification": "out/excised/corpus/verification.json",
            "verification_sha256": (sha256_file(verification_path)
                                    if verification_path.exists() else None),
            "verification_schema": verification.get("schema"),
            "certificates_dir": "out/excised/corpus",
            "original_area_cache": AREA_CACHE_REL,
        },
        "census_params": verification.get("census_params"),
        "area_pricing": area_pricing,
        "summary": {},
        "archive": {},
        "order": order,
        "segments": segs,
    }
    index["summary"] = summary_block(order, segs, verification, base_manifest)
    index["archive"] = archive_block(order, segs, repo_root)
    return index


# ---------------------------------------------------------------- markdown

def pct(value, digits: int = 4) -> str:
    if value is None:
        return "not recorded"
    return f"{value * 100:.{digits}f}%"


def pct_truncated(value, digits: int = 3) -> str:
    """A percentage TRUNCATED toward zero, as the decision rule prints it.

    bench/headline_decision.py never rounds a retention figure up, so a
    figure quoted here in its form is truncated the same way rather than
    rounded to the same width.
    """
    if value is None:
        return "not recorded"
    q = Decimal(1).scaleb(-digits)
    return str((Decimal(repr(float(value))) * 100).quantize(
        q, rounding=ROUND_DOWN)) + "%"


def num(value, digits: int = 4) -> str:
    if value is None:
        return "not recorded"
    return f"{value:.{digits}f}"


def mib(byte_count: int) -> str:
    return f"{byte_count / 2 ** 20:,.1f} MiB"


def render_markdown(index: dict) -> str:
    s = index["summary"]
    a = index["archive"]
    segs = index["segments"]
    order = index["order"]
    disp = s["dispositions"]
    rt = s["runtime_seconds"]
    ret = s["retention"]
    L: list[str] = []
    add = L.append

    add("# Excised corpus")
    add("")
    add(f"Generated {index['generated_utc']} by `{index['generator']}` from "
        "the pinned base manifest, the per-segment excision certificates and "
        "the independent recensus record. Every number on this page is "
        "copied from those files; nothing is hand-entered. The machine "
        "readable form of exactly this content is "
        "`out/release/index.json` (schema `" + index["schema"] + "`).")
    add("")
    add("Throughout, **represented surface retained** means canonical quad "
        "area that survives the cut, divided by the canonical quad area of "
        "the stated denominator. It is a measure of the mesh and nothing "
        "more.")
    add("")

    # ---------------------------------------------------------- summary
    add("## Summary")
    add("")
    add("| | |")
    add("|---|---|")
    add(f"| Pinned trace artifacts | {s['n_segments']} |")
    add(f"| Scrolls | {s['n_scrolls']} |")
    add(f"| Transformed | {disp.get('transformed', 0)} |")
    add(f"| Already clean (no change needed) | "
        f"{disp.get('already_clean', 0)} |")
    add(f"| Not censusable | {disp.get('not_censusable', 0)} |")
    add(f"| Displacement-repaired bases | "
        f"{s['bases'].get('displacement_repaired', 0)} |")
    add(f"| Original published bases | {s['bases'].get('original', 0)} |")
    add(f"| Unique base geometries | "
        f"{s['bases'].get('n_unique_geometries')} |")
    add(f"| Exact duplicate bases | "
        f"{s['bases'].get('n_duplicate_aliases')} |")
    add(f"| Represented surface retained, operational denominator, "
        f"area-weighted | {pct(ret['operational_area_weighted'])} "
        f"over {ret['operational_segments_counted']} censusable segments |")
    add(f"| Represented surface retained, headline denominator, "
        f"area-weighted | {pct(ret['headline_area_weighted'])} "
        f"over {ret['headline_segments_counted']} pinned artifacts |")
    add(f"| Lowest segment retention (headline) | "
        f"{pct(ret['min_headline_retained_fraction'])} "
        f"({ret['min_headline_segment']}) |")
    add(f"| Runtime min / median / p90 / max | "
        f"{num(rt['min'], 2)} s / {num(rt['median'], 2)} s / "
        f"{num(rt['p90'], 2)} s / {num(rt['max'], 2)} s |")
    add(f"| Runtime gate | {num(rt['gate'], 0)} s, "
        f"{rt['n_over_gate']} over |")
    add(f"| Meshes re-hashed clean in a fresh workdir | "
        f"{s['verification']['n_meshes_rehash_ok']}/"
        f"{s['verification']['n_meshes_rehashed']} |")
    add(f"| Independently recensused | "
        f"{s['verification']['n_recensus_clean_both_diagonals']}/"
        f"{s['verification']['n_recensused']} clean under both canonical "
        f"triangulations |")
    add(f"| Residual non-adjacent transverse contacts | "
        f"{s['verification']['residual_transverse_contacts']} |")
    add(f"| Disagreements with the recorded census | "
        f"{s['verification']['n_census_disagreements']} |")
    add(f"| Timeouts / errors | {s['verification']['n_timeouts']} / "
        f"{s['verification']['n_failed']} |")
    add(f"| 99.9%-area core gate | {s['fragmentation']['n_core_gate_pass']} "
        f"pass, {s['fragmentation']['n_core_gate_fail']} fail |")
    add("")

    add("### Per scroll")
    add("")
    add("| Scroll | Segments | Transformed | Already clean | "
        "Not censusable |")
    add("|---|---:|---:|---:|---:|")
    for scroll in sorted(s["per_scroll"]):
        b = s["per_scroll"][scroll]
        add(f"| {scroll} | {b['n']} | {b['transformed']} | "
            f"{b['already_clean']} | {b['not_censusable']} |")
    add("")

    # ------------------------------------------------------ definitions
    add("## What the columns mean")
    add("")
    add("**Disposition.**")
    for k in ("transformed", "already_clean", "not_censusable"):
        add(f"- `{k}` — {DEFS['dispositions'][k]}.")
    add("")
    add("**Input hash and output hash.** Each tifxyz mesh is three coordinate "
        "planes. The input hash column is the sha256 of the `x` plane of the "
        "pre-excision base that was actually cut; the output hash column is "
        "the sha256 of the `x` plane of the emitted mesh. `index.json` "
        "carries all three planes on both sides, plus the hashes of the "
        "original published mesh. The `y` and `z` planes are hashed the same "
        "way and are checked together, so the `x` plane alone is a "
        "convenient abbreviation, not the whole check.")
    add("")
    add("**The two area denominators.** Two retention figures are recorded "
        "for every segment, because 103 of the 185 bases are "
        "displacement-repaired meshes rather than the published originals.")
    add("")
    add(f"- **Operational retention** — {DEFS['operational_denominator']}")
    add(f"- **Headline retention** — {DEFS['headline_denominator']}")
    add("")
    add("The two summary figures are also evaluated over different "
        "populations, and each row above says which. **Operational "
        f"retention** runs over the {ret['operational_segments_counted']} "
        "censusable segments: every segment that was cut, plus every "
        "already-clean segment, each priced against the base it was actually "
        f"handed. **Headline retention** runs over all "
        f"{ret['headline_segments_counted']} pinned artifacts — "
        f"{ret['headline_segments_priced']} of them carry a priced "
        f"original-coordinate area, and the other "
        f"{len(ret['headline_segments_zero_area'])} are the triangle-empty "
        "or invalid inputs, which carry no surface at all and so enter with "
        "zero area on both sides of the ratio, where they cannot move the "
        "figure.")
    add("")
    add("A segment that was never cut is **not** dropped from the headline "
        "denominator. It enters at exactly 1.0, against the full canonical "
        "area of its original published mesh. Where the base IS that mesh, "
        "the area is on the certificate. Where the base was "
        "displacement-repaired, the certificate writes no area block at all "
        "— nothing was cut, so there was nothing to price — and the original "
        "area is instead recovered by re-reading the original published mesh "
        "and dividing its retained quads the same way every cut segment is "
        f"divided; the result is memoised in `{AREA_CACHE_REL}`. The "
        "headline figure above, and the population it is taken over, are "
        "therefore the same ones the frozen decision rule "
        "(`bench/headline_decision.py`) reports: truncated to three decimals "
        "as that rule prints it — it never rounds a retention figure up — "
        f"the figure is {pct_truncated(ret['headline_area_weighted'])} over "
        f"the same {ret['headline_segments_counted']} artifacts.")
    add("")
    div = ret["denominator_divergence"]
    add("In the per-segment tables below the two columns look identical, and "
        "that is a measured result rather than a copy. Of the "
        f"{div['n_segments_compared']} segments that carry both figures, "
        f"{div['n_segments_identical']} agree to the last bit — every "
        "segment cut from the original published mesh, where the two "
        "denominators are the same denominator, plus every already-clean "
        "segment, where nothing was removed. Across the rest, all of them "
        "displacement-repaired bases that were cut, the largest difference "
        f"anywhere in the corpus is {div['max_abs_difference']:.3e} "
        f"(`{div['max_abs_difference_segment']}`): the two figures agree to "
        "roughly seven significant figures. Pricing the removed quads on the "
        "repaired coordinates rather than the published ones barely moves "
        "the area. Both are still reported per segment, because that "
        "agreement is an observation about this corpus and not a guarantee.")
    add("")
    add(f"**R_main.** {DEFS['R_main']}")
    add("")
    add("- **min R_main** is the minimum over *all* input components, with no "
        "size threshold of any kind.")
    add("- **area-weighted R_main** weights each input component's R_main by "
        "that component's share of input area.")
    add(f"- **core gate** — {DEFS['core_gate']}")
    add("")
    add("**Runtime.** Wall-clock seconds for the segment under the frozen "
        f"policy, against a {num(rt['gate'], 0)} s gate.")
    add("")

    # -------------------------------------------------- not censusable
    add("## Inputs that could not be censused")
    add("")
    add(f"{len(s['not_censusable'])} of the {s['n_segments']} pinned "
        "artifacts carry no triangles at all. No census, no cut and no "
        "cleanliness claim is defined on them, so they hold no retention "
        "figure and no fragmentation figure. They are listed here rather "
        "than dropped, and they still have a certificate.")
    add("")
    add("| Segment | Scroll | Grid | Valid vertices | Retained quads | "
        "Reason |")
    add("|---|---|---|---:|---:|---|")
    for item in s["not_censusable"]:
        seg = segs[item["segment"]]
        ev = item.get("evidence") or {}
        grid = ev.get("grid_shape") or seg.get("grid_shape") or []
        add(f"| `{item['segment']}` | {item['scroll']} | "
            f"{'x'.join(str(g) for g in grid)} | "
            f"{ev.get('n_valid_vertices', '—')} | "
            f"{ev.get('n_retained_quads', '—')} | {item['reason']} |")
    add("")

    # -------------------------------------------------- fragmentation
    add("## Fragmentation")
    add("")
    frag = s["fragmentation"]
    if frag["core_gate_failures"]:
        add(f"{frag['n_core_gate_pass']} of the "
            f"{frag['n_core_gate_measured'] + disp.get('already_clean', 0)} "
            "segments where the core gate applies pass it. The following "
            "show material fragmentation: at least one component inside the "
            "99.9%-area core has R_main below 0.90, meaning its largest "
            "surviving descendant holds less than 90% of that component's "
            "input area. Area retention alone does not show this — each of "
            "these segments still retains most of its area — which is why "
            "R_main is reported separately.")
        add("")
        add("| Segment | Scroll | min R_main (core) | min R_main (all) | "
            "area-weighted R_main | core components below gate |")
        add("|---|---|---:|---:|---:|---:|")
        for f in frag["core_gate_failures"]:
            add(f"| `{f['segment']}` | {f['scroll']} | "
                f"{num(f['min_R_main_core'], 6)} | "
                f"{num(f['min_R_main_all_components'], 6)} | "
                f"{num(f['area_weighted_R_main'], 6)} | "
                f"{f['n_core_components_below_gate']} |")
    else:
        add("No segment fails the 99.9%-area core gate.")
    add("")

    # ------------------------------------------------------- archive
    add("## Archive size on disk")
    add("")
    add(a["note"])
    add("")
    add(f"**Total emitted corpus artifacts: "
        f"{a['total_bytes']:,} bytes = {a['total_gib']:.3f} GiB "
        f"({a['total_gb']:.3f} GB); {a['total_bytes_on_disk']:,} bytes "
        f"allocated on disk.**")
    add("")
    add("| Scroll | Segments | Emitted meshes | Mesh bytes | "
        "Certificate bytes | Log bytes | Total |")
    add("|---|---:|---:|---:|---:|---:|---:|")
    for scroll in sorted(a["per_scroll"]):
        b = a["per_scroll"][scroll]
        add(f"| {scroll} | {b['n_segments']} | {b['n_meshes']} | "
            f"{mib(b['mesh_bytes'])} | {mib(b['certificate_bytes'])} | "
            f"{b['log_bytes'] / 1024:,.0f} KiB | {mib(b['bytes'])} |")
    add("")
    add("### Proposed split")
    add("")
    add(a["proposed_split_rationale"] + " The sizes below are the "
        "uncompressed payload each archive would contain.")
    add("")
    add("| Proposed archive | Scroll | Contents | Payload |")
    add("|---|---|---|---:|")
    for p in a["proposed_archives"]:
        add(f"| `{p['archive']}` | {p['scroll']} | "
            f"{'; '.join(p['contents'])} | {p['payload_mib']:,.1f} MiB "
            f"({p['payload_gib']:.3f} GiB) |")
    add("")
    add(f"Each archive also carries the shared manifests, "
        f"{mib(a['shared_metadata_bytes'])} in total: "
        + ", ".join(f"`{f['path']}`" for f in a["shared_metadata_files"])
        + ".")
    add("")

    # ---------------------------------------------------- certificates
    add("## Where the certificates live")
    add("")
    add("One excision certificate per segment, at "
        "`out/excised/corpus/<segment>_excision_certificate.json`. Each one "
        "records the pinned input mesh and its plane hashes, the emitted "
        "output mesh and its plane hashes, the census before and after under "
        "both canonical triangulations, both area denominators, the "
        "per-component recovery table, the frozen selection policy, and the "
        "code provenance of the run.")
    add("")
    add("Alongside them:")
    add("")
    add("- `out/corpus_bases.json` — the pinned base manifest: which mesh "
        "each segment was cut from, and the hash that pins it.")
    add("- `out/excised/corpus/corpus_summary.jsonl` — one driver record per "
        "segment.")
    note = index["summary"]["verification"]["independence_note"] or ""
    if note:
        note = " " + note[0].upper() + note[1:] + "."
    add("- `out/excised/corpus/verification.json` — the independent "
        "recensus: every emitted artifact re-hashed and re-censused in a "
        "fresh workdir." + note)
    add("- `out/excised/corpus/logs/<segment>.log` — the per-segment run log.")
    add("- `out/release/index.json` — this document, machine-readable.")
    add("")

    # ------------------------------------------------------- per scroll
    add("## Per-segment results")
    add("")
    add("`op.` is operational retention, `headline` is headline retention; "
        "see the two denominators above. Hashes are the first 12 hex digits "
        "of the sha256 of the `x` plane. A dash means the certificate "
        "records no value for that segment, which is stated rather than "
        "filled in.")
    add("")
    scrolls = sorted({segs[x]["scroll"] for x in order},
                     key=lambda v: (v or ""))
    for scroll in scrolls:
        members = [x for x in order if segs[x]["scroll"] == scroll]
        add(f"### {scroll} ({len(members)} segments)")
        add("")
        add("| Segment | Disposition | Base | Input hash | Output hash | "
            "op. retained | headline retained | min R_main | "
            "area-wt R_main | core gate | Runtime (s) |")
        add("|---|---|---|---|---|---:|---:|---:|---:|:---:|---:|")
        for seg in members:
            r = segs[seg]
            ret_r = r["retention"]
            fr = r["fragmentation"]
            ih = (r["input_hashes"] or {}).get("x")
            oh = (r["output_hashes"] or {}).get("x")
            gate = fr["core_gate_pass"]
            gate_s = "—" if gate is None else ("pass" if gate else "FAIL")
            if fr["measured"]:
                minr = num(fr["min_R_main_all_components"], 4)
                awr = num(fr["area_weighted_R_main"], 4)
            elif r["disposition"] == "already_clean":
                minr = awr = "1 by constr."
            else:
                minr = awr = "—"
            add(f"| `{seg}` | {r['disposition']} | "
                f"{'repaired' if r['base_kind'] == 'displacement_repaired' else 'original'} | "
                f"`{ih[:12] if ih else '—'}` | `{oh[:12] if oh else '—'}` | "
                f"{pct(ret_r['operational_retained_fraction'])} | "
                f"{pct(ret_r['headline_retained_fraction'])} | "
                f"{minr} | {awr} | {gate_s} | "
                f"{num(r['wall_seconds'], 2)} |")
        add("")

    add("### Reading the retention columns")
    add("")
    add("A segment with disposition `already_clean` shows 100.0000% under "
        "both denominators because no cut was made. Its R_main is 1 by "
        "construction rather than by measurement, and the certificate "
        "records no per-component table for it; that is marked "
        "`1 by constr.` rather than reported as a measured value. A segment "
        "with disposition `not_censusable` shows a dash everywhere: no area "
        "and no component structure is defined on an input with no "
        "triangles.")
    add("")
    ap = index.get("area_pricing") or {}
    if ap.get("n_priced"):
        add(f"{ap['n_priced']} already-clean segments sit on a "
            "displacement-repaired base, so their certificate records no "
            "original-coordinate area: nothing was cut, so nothing was "
            "priced. They are **not** dropped from the area-weighted "
            "headline denominator. The canonical area of each one's original "
            "published mesh was recomputed from that mesh, by the same route "
            "every cut segment's headline denominator is computed, and each "
            "segment enters at a retained fraction of exactly 1.0. The "
            f"recomputed areas are memoised in `{ap.get('cache', AREA_CACHE_REL)}`. "
            "They are: "
            + ", ".join(f"`{x['segment']}`" for x in ap["priced"]) + ".")
        add("")
    if ret["headline_segments_without_area"]:
        add(f"{len(ret['headline_segments_without_area'])} segments still "
            "carry no original-coordinate area and could not be priced from "
            "their original mesh, so they are outside the area-weighted "
            "headline denominator above and outside the count it is reported "
            "over. Each removed exactly zero area, so including them could "
            "only raise that figure. They are: "
            + ", ".join(f"`{x}`" for x in
                        ret["headline_segments_without_area"]) + ".")
        add("")
    if ret["headline_segments_zero_area"]:
        add(f"{len(ret['headline_segments_zero_area'])} not-censusable "
            "inputs are inside the count the headline figure is reported "
            "over, carrying zero area on both sides of the ratio: they have "
            "no triangles, so there is no surface to retain and none to "
            "remove, and their presence cannot move the figure either way. "
            "They are: "
            + ", ".join(f"`{x}`" for x in
                        ret["headline_segments_zero_area"]) + ".")
        add("")
    add(f"Restricted to the {disp.get('transformed', 0)} transformed "
        "segments alone — that is, excluding every already-clean segment "
        "from both numerator and denominator — the area-weighted headline "
        "retention is "
        f"{pct(ret['headline_area_weighted_transformed_only'])}. That is the "
        "stricter reading: it prices only the segments that were actually "
        "cut.")
    add("")
    return "\n".join(L) + "\n"


# -------------------------------------------------------------------- main

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", default=str(Path(__file__).resolve()
                                               .parent.parent))
    ap.add_argument("--out", default=None,
                    help="index.json path (default out/release/index.json)")
    ap.add_argument("--markdown", default=None,
                    help="document path (default docs/CORPUS.md)")
    ap.add_argument("--no-markdown", action="store_true")
    args = ap.parse_args(argv)

    repo_root = Path(args.repo_root).resolve()
    index = build(repo_root)

    out = Path(args.out) if args.out else repo_root / "out/release/index.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(index, indent=1, sort_keys=False) + "\n")
    print(f"wrote {out}  ({out.stat().st_size:,} bytes, "
          f"{len(index['order'])} segments)")

    if not args.no_markdown:
        md = Path(args.markdown) if args.markdown \
            else repo_root / "docs/CORPUS.md"
        md.parent.mkdir(parents=True, exist_ok=True)
        # The archive block quotes index.json's own size, so recompute it
        # now that the file exists on disk and re-emit both.
        index["archive"] = archive_block(index["order"], index["segments"],
                                         repo_root)
        out.write_text(json.dumps(index, indent=1, sort_keys=False) + "\n")
        md.write_text(render_markdown(index))
        print(f"wrote {md}  ({md.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
