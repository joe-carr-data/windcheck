#!/usr/bin/env python3
"""Build the release index for the EXPANSION corpus, plus combined totals.

Additive companion to bench/build_release_index.py, exactly as
out/expand_bases.json is the additive companion to out/corpus_bases.json:
the pinned index (out/release/index.json, 185 segments) is READ, never
restated, and this program writes

  * out/release/expansion_index.json   one record per expansion segment
  * docs/CORPUS-EXPANSION.md           the same content, human-readable

Every per-segment number is copied out of the segment's own excision
certificate or the independent verification record; nothing is
hand-entered. The combined block recomputes corpus-wide retention over
BOTH inventories from per-segment areas -- the pinned 99.505% figure is
never carried forward by copy.

Denominator reconciliation (stated here once, verbatim in the output):
  185 pinned segments (179 censusable) + 99 expansion traces
  (95 censusable) = 284 indexed artifacts, 274 censusable. The July
  full-corpus census reported "278 traces" = 179 pinned censusable
  + 99 expansion traces processed; the 4 expansion refusals below the
  census validity floor and the 6 pinned not-censusable inputs are why
  278 != 284 and 274 != 278.

For expansion segments base_kind is "original" throughout (no
displacement-repair campaign ran), so the operational and headline
denominators coincide by construction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
for _p in (REPO_ROOT / "src",):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

SCHEMA = "expansion_index/v1"

RECONCILIATION = (
    "185 pinned segments (179 censusable) + 99 expansion traces (95 "
    "censusable) = 284 indexed artifacts, 274 censusable. The July "
    "full-corpus census reported 278 traces = 179 pinned censusable + 99 "
    "expansion traces processed; the 4 expansion refusals below the census "
    "validity floor and the 6 pinned not-censusable inputs account for "
    "278 != 284 and 274 != 278.")


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build(repo_root: Path) -> dict:
    out_dir = repo_root / "out" / "excised" / "expand"
    bases = json.loads((repo_root / "out" / "expand_bases.json").read_text())
    ver = json.loads((out_dir / "verification.json").read_text())
    ncc = json.loads((out_dir / "not_censusable_confirmed.json").read_text())
    ver_by_seg = {s["segment"]: s for s in ver["segments"]}

    if ver["summary"]["n_failed"] != 0:
        raise SystemExit("verification records failures; refusing to index")

    segments = {}
    order = []
    for e in sorted(bases["entries"], key=lambda e: e["segment"]):
        seg = e["segment"]
        cpath = out_dir / f"{seg}_excision_certificate.json"
        cert = json.loads(cpath.read_text())
        disp = cert["terminal_disposition"]
        v = ver_by_seg[seg]
        rec = {
            "segment": seg,
            "sample": e["corpus"],
            "volume": e.get("volume"),
            "voxel_um": e.get("voxel_um"),
            "voxel_um_note": e.get("voxel_um_note"),
            "s3_prefix": e.get("s3_prefix"),
            "disposition": disp,
            "base_kind": "original",
            "input_mesh": e["base_mesh"],
            "input_hashes": e["base_hashes"],
            "output_mesh": cert.get("output_mesh")
            if disp == "transformed" else None,
            "output_hashes": cert.get("output_mesh_hashes")
            if disp == "transformed" else None,
            "certificate": str(cpath.relative_to(repo_root)),
            "certificate_sha256": sha256_file(cpath),
            "geometry_key": e["geometry_key"],
            "is_canonical": e["is_canonical"],
            "duplicate_of": e.get("duplicate_of"),
            "grid_shape": cert.get("grid_shape"),
            "n_valid_vertices": cert.get("n_valid"),
            "wall_seconds": cert.get("wall_seconds"),
            "policy_version": cert.get("policy_version"),
            "policy_hash": cert.get("policy_hash"),
            "verification_ok": v["ok"],
        }
        if disp == "transformed":
            cr = cert.get("component_recovery") or {}
            rec["retention"] = {
                "retained_fraction": cert.get("headline_retained_fraction"),
                "A_original_canonical": cert.get("input_area_canonical"),
                "note": ("base is the original published mesh, so the "
                         "operational and headline denominators coincide"),
            }
            rec["fragmentation"] = {
                "components_before": cr.get("components_before"),
                "components_after": cr.get("components_after"),
                "core_gate_pass": cert.get("core_gate_pass"),
                "min_R_main_any_component":
                    cr.get("min_R_main_any_component"),
            }
            rec["census"] = {
                "before_transverse_total":
                    cert.get("input_transverse_total"),
                "after_transverse_total":
                    cert.get("output_transverse_total"),
                "independent_recensus_clean_both_diagonals":
                    v["recensus"].get("clean_both_diagonals"),
            }
        elif disp == "already_clean":
            rec["retention"] = {
                "retained_fraction": 1.0,
                "A_original_canonical": cert.get("input_area_canonical"),
                "note": "no cut was made",
            }
            rec["census"] = {
                "before_transverse_total": 0,
                "independent_recensus_clean_both_diagonals":
                    v["recensus"].get("clean_both_diagonals"),
            }
        else:   # not_censusable
            rec["retention"] = None
            rec["census"] = None
            rec["not_censusable"] = {
                "reason": ("below the census validity floor (5000 valid "
                           "cells); no audit or excision is defined"),
                "n_valid": cert.get("n_valid"),
                "decline_confirmed":
                    ncc["segments"][seg]["declined"],
            }
        segments[seg] = rec
        order.append(seg)

    # ---- expansion-only retention (area-weighted, per-segment areas) ----
    a_tot = a_ret = 0.0
    counted = 0
    for seg in order:
        r = segments[seg]
        if r["retention"] is None or not r["retention"].get(
                "A_original_canonical"):
            continue
        a = float(r["retention"]["A_original_canonical"])
        a_tot += a
        a_ret += a * float(r["retention"]["retained_fraction"])
        counted += 1
    exp_retention = {
        "area_weighted_retained_fraction": a_ret / a_tot,
        "A_original_total_canonical": a_tot,
        "segments_counted": counted,
        "worst_segment": min(
            (s for s in order if segments[s]["retention"]),
            key=lambda s: segments[s]["retention"]["retained_fraction"]),
    }
    ws = exp_retention["worst_segment"]
    exp_retention["worst_segment_retained"] = (
        segments[ws]["retention"]["retained_fraction"])
    exp_retention["worst_segment_core_gate_pass"] = (
        (segments[ws].get("fragmentation") or {}).get("core_gate_pass"))

    # ---- combined totals over BOTH inventories (recomputed, not copied) --
    pinned = json.loads(
        (repo_root / "out" / "release" / "index.json").read_text())
    ps = pinned["summary"]["retention"]
    comb_A = float(ps["headline_A_original_total"]) + a_tot
    comb_removed = float(ps["headline_A_removed_total"]) + (a_tot - a_ret)
    combined = {
        "what": ("area-weighted retention over the pinned AND expansion "
                 "inventories, recomputed from per-segment areas"),
        "pinned_index_sha256": sha256_file(
            repo_root / "out" / "release" / "index.json"),
        "n_indexed": 185 + len(order),
        "n_censusable": 179 + counted,
        "dispositions": {
            "pinned": pinned["summary"]["dispositions"],
            "expansion": {
                "transformed": sum(1 for s in order
                                   if segments[s]["disposition"]
                                   == "transformed"),
                "already_clean": sum(1 for s in order
                                     if segments[s]["disposition"]
                                     == "already_clean"),
                "not_censusable": sum(1 for s in order
                                      if segments[s]["disposition"]
                                      == "not_censusable"),
            },
        },
        "combined_area_weighted_retained_fraction":
            1.0 - comb_removed / comb_A,
        "combined_A_original_total": comb_A,
        "combined_A_removed_total": comb_removed,
        "area_unit_note": (
            "areas are in canonical grid vx^2 of each trace's own volume; "
            "voxel sizes differ across scrolls, so the combined figure is a "
            "vx^2-weighted summary, not a physical-area claim"),
        "reconciliation": RECONCILIATION,
    }

    return {
        "schema": SCHEMA,
        "generated_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "sources": {
            "base_manifest": "out/expand_bases.json",
            "base_manifest_sha256": sha256_file(
                repo_root / "out" / "expand_bases.json"),
            "verification": "out/excised/expand/verification.json",
            "verification_sha256": sha256_file(out_dir /
                                               "verification.json"),
            "census_records": "out/expand/expand.jsonl",
        },
        "census_params": ver["census_params"],
        "reconciliation": RECONCILIATION,
        "expansion_retention": exp_retention,
        "combined": combined,
        "order": order,
        "segments": segments,
    }


def render_markdown(idx: dict) -> str:
    L = []
    L.append("# The expansion corpus: 99 further traces, dispositioned "
             "and verified\n")
    L.append("Additive companion to `docs/CORPUS.md` (the 185 pinned "
             "segments). Same frozen transform policy, same census "
             "parameters, same independent verification. Nothing in the "
             "pinned corpus is restated.\n")
    L.append("## Reconciliation of the denominators\n")
    L.append(idx["reconciliation"] + "\n")
    er = idx["expansion_retention"]
    cb = idx["combined"]
    L.append("## Retention\n")
    L.append(f"- Expansion inventory: **{er['area_weighted_retained_fraction']:.6f}** "
             f"area-weighted over {er['segments_counted']} censusable "
             "traces.")
    L.append(f"- Worst single trace: `{er['worst_segment']}` retains "
             f"**{er['worst_segment_retained']:.4f}** and its 99.9%-core "
             f"gate is **{'PASS' if er['worst_segment_core_gate_pass'] else 'FAIL'}** "
             "-- reported, not averaged away.")
    L.append(f"- Combined over both inventories (recomputed from "
             f"per-segment areas, {cb['n_censusable']} censusable of "
             f"{cb['n_indexed']} indexed): "
             f"**{cb['combined_area_weighted_retained_fraction']:.6f}**.")
    L.append(f"- {cb['area_unit_note']}\n")
    L.append("## Per-trace records\n")
    L.append("| trace | sample | disposition | retained | core gate | "
             "verified |")
    L.append("|---|---|---|---|---|---|")
    for seg in idx["order"]:
        r = idx["segments"][seg]
        ret = (f"{r['retention']['retained_fraction']:.6f}"
               if r["retention"] else "--")
        core = ""
        if r.get("fragmentation"):
            core = "PASS" if r["fragmentation"]["core_gate_pass"] else "FAIL"
        L.append(f"| `{seg}` | {r['sample']} | {r['disposition']} | {ret} "
                 f"| {core} | {'yes' if r['verification_ok'] else 'NO'} |")
    L.append("\nEvery number above is copied from the trace's own "
             "excision certificate or the independent verification record "
             "(`verification.json`); the four census refusals carry an "
             "explicit decline confirmation "
             "(`not_censusable_confirmed.json`).\n")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    ap.add_argument("--no-markdown", action="store_true")
    args = ap.parse_args(argv)
    idx = build(args.repo_root)
    out = args.repo_root / "out" / "release" / "expansion_index.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(idx, indent=1))
    print(f"{out}: {len(idx['order'])} segments; expansion retention "
          f"{idx['expansion_retention']['area_weighted_retained_fraction']:.6f}; "
          f"combined {idx['combined']['combined_area_weighted_retained_fraction']:.6f} "
          f"over {idx['combined']['n_censusable']} censusable")
    if not args.no_markdown:
        md = args.repo_root / "docs" / "CORPUS-EXPANSION.md"
        md.write_text(render_markdown(idx))
        print(f"wrote {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
