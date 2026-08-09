"""PREREG B7 acceptance test: the evaluator must agree with the certificates.

The benchmark is not published unless, for every packaged trace, scoring
the REFERENCE DERIVATIVE as if it were a candidate reproduces what the
segment's own certificate already says: the same clean verdict, the same
retained fraction, the same removed-cell count, the same fragmentation,
and bit-identical retained coordinates.

That is the whole point of the test. The evaluator and the certificates
were written months apart by different code paths; if they disagree, one
of them is wrong and neither number can be published.

Two modes:

  --emit-manifest   write the TSV that bench/census_batch.sh consumes, so
                    every reference surface is censused in ONE container.
                    Per-surface containers wedged the Docker daemon at
                    3 of 154 last time this was done the naive way.

  --reports DIR     run the comparison, taking the census from that sweep.

The census is by the OFFICIAL validator, not ours. Our own engine agreeing
with our own certificates would prove nothing about interoperability --
and the one defect this project shipped was invisible to our engine and
obvious to the official one.

TOLERANCES, fixed here and not tuned afterwards:

  retained_fraction   1e-9 relative. The evaluator recomputes canonical
                      quad area in float64 from the same arrays, so this
                      is float-summation order, nothing else.
  removed cells       exact. An off-by-one is a real disagreement.
  fragmentation       1e-9 absolute on R_main.
  clean               exact.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from topology_eval import evaluate  # noqa: E402

REL_TOL = 1e-9
ABS_TOL = 1e-9


def surfaces(bench_index: Path, indices: list[Path]) -> list[dict]:
    idx = json.loads(bench_index.read_text())
    entries = {}
    for p in indices:
        d = json.loads(p.read_text())
        entries.update(d["segments"])
    rows = []
    for t in idx["traces"]:
        if not t.get("scored"):
            continue
        e = entries[t["segment"]]
        ref = e.get("output_mesh") or e["input_mesh"]
        rows.append({"segment": t["segment"], "input": e["input_mesh"],
                     "reference": ref,
                     "certificate": e["certificate"],
                     "provenance": str(bench_index.parent / "traces"
                                       / t["segment"] / "provenance.json"),
                     "is_identity": e.get("output_mesh") is None})
    return rows


def close(a, b, rel=REL_TOL) -> bool:
    if a is None or b is None:
        return a is b
    return abs(a - b) <= rel * max(1.0, abs(a), abs(b))


def agree(out: dict, key: str, ev, ce, *, tol: float | None = None) -> None:
    """Record one comparison, and FAIL when only one side has a value.

    The earlier version skipped whenever either side was None. A broken
    evaluator that emitted `None` for a field therefore passed against a
    certificate that had a real value -- the check disappeared exactly
    when it was needed.
    """
    out["evaluator"][key] = ev
    out["certificate"][key] = ce
    if ce is None:
        if ev is not None:
            out["not_compared"].append(
                f"{key}: the certificate carries no value")
        return
    if ev is None:
        out["disagreements"].append(
            f"{key}: the evaluator produced no value while the certificate "
            f"says {ce!r}")
        return
    if tol is None:
        if bool(ev) != bool(ce) if isinstance(ce, bool) else ev != ce:
            out["disagreements"].append(
                f"{key}: evaluator {ev!r}, certificate {ce!r}")
    elif abs(ev - ce) > tol * max(1.0, abs(ev), abs(ce)):
        out["disagreements"].append(
            f"{key}: evaluator {ev!r}, certificate {ce!r}")


def compare(row: dict, report: Path, binding: Path | None) -> dict:
    cert = json.loads(Path(row["certificate"]).read_text())
    r = evaluate(Path(row["input"]), Path(row["reference"]), [],
                 census_report=report, census_binding=binding,
                 provenance=(Path(row["provenance"])
                             if row.get("provenance") else None))
    out = {"segment": row["segment"], "disagreements": [],
           "not_compared": [], "evaluator": {}, "certificate": {}}

    if not r.get("comparable") or r.get("refused"):
        out["disagreements"].append(
            f"the evaluator refused the reference: "
            f"{r.get('refused') or 'not comparable'}")
        return out

    # ---- clean -----------------------------------------------------------
    cert_clean = bool(cert.get("claimed_clean"))
    out["evaluator"]["clean"] = r.get("clean")
    out["certificate"]["clean"] = cert_clean
    out["census_source"] = r["census"].get("source")
    out["census_transverse"] = r["census"].get("transverse")
    out["census_triangles"] = r["census"].get("triangles")
    if bool(r.get("clean")) != cert_clean:
        out["disagreements"].append(
            f"clean: evaluator {r.get('clean')}, certificate {cert_clean} "
            f"({r['census'].get('transverse')})")
    # A clean verdict over zero triangles is how the mask defect hid.
    tri = r["census"].get("triangles") or {}
    if r.get("clean") and all((v or 0) == 0 for v in tri.values()):
        out["disagreements"].append(
            "clean over ZERO triangles: the official loader sees no "
            "surface here, so the verdict is vacuous")
    if not r.get("scored"):
        return out

    # ---- retention -------------------------------------------------------
    agree(out, "retained_fraction", r.get("retained_fraction"),
          cert.get("operational_retained_fraction"), tol=REL_TOL)

    # ---- removed cells ---------------------------------------------------
    cert_removed = (cert.get("excision") or {}).get("n_invalidated_vertices")
    if cert_removed is None and row["is_identity"]:
        cert_removed = 0
    agree(out, "removed_cells", r.get("removed_cells"), cert_removed)
    agree(out, "added_cells", r.get("added_cells"),
          0 if cert.get("terminal_disposition") else None)

    # ---- fragmentation ---------------------------------------------------
    # A trace that needed no cut has no component_recovery block, so its
    # expected values are stated HERE rather than left absent: nothing was
    # removed, so every input component survives whole.
    cr = cert.get("component_recovery") or {}
    cg = cr.get("core_gate") or {}
    un = cr.get("unthresholded") or {}
    if row["is_identity"] and not cr:
        # NOT derived from the evaluator's own output -- comparing a tool
        # to itself is a check that cannot fail. The certificate gives no
        # component COUNT for an uncut trace, so that stays uncompared and
        # is reported as such; what an uncut trace does pin exactly is
        # that every component survives whole.
        exp = {"components_input": None,
               "components_candidate": None,
               "R_main_area_weighted": 1.0, "R_main_min": 1.0,
               "n_core_components": None, "min_R_main_core": 1.0,
               "core_gate_pass": bool(cert.get("core_gate_pass", True))}
        out["identity_expectation"] = (
            "no cut was made, so R_main = 1 for every input component and "
            "the core gate passes by construction. The certificate records "
            "no component COUNT for such a trace, so the counts are "
            "reported as not compared rather than checked against a number "
            "this tool produced itself.")
    else:
        exp = {"components_input": cg.get("n_input_components"),
               "components_candidate": None,
               "R_main_area_weighted": un.get("area_weighted_R_main"),
               "R_main_min": un.get("min_R_main_all_components"),
               "n_core_components": cg.get("n_core_components"),
               "min_R_main_core": cg.get("min_R_main_core"),
               "core_gate_pass": cg.get("core_gate_pass")}
    agree(out, "components_input", r.get("components_input"),
          exp["components_input"])
    agree(out, "components_candidate", r.get("components_candidate"),
          exp["components_candidate"])
    agree(out, "R_main_area_weighted", r.get("R_main_area_weighted"),
          exp["R_main_area_weighted"], tol=ABS_TOL)
    agree(out, "R_main_min", r.get("R_main_min"), exp["R_main_min"],
          tol=ABS_TOL)
    core = r.get("core_gate") or {}
    agree(out, "n_core_components", core.get("n_core_components"),
          exp["n_core_components"])
    agree(out, "min_R_main_core", core.get("min_R_main_core"),
          exp["min_R_main_core"], tol=ABS_TOL)
    agree(out, "core_gate_pass", core.get("core_gate_pass"),
          exp["core_gate_pass"])

    # ---- coordinate fidelity --------------------------------------------
    cert_fid = (cert.get("reload_checks") or {}).get(
        "retained_coordinate_bit_identity")
    if cert_fid is None and row["is_identity"]:
        cert_fid = True
    agree(out, "coordinate_fidelity", r.get("coordinate_fidelity"), cert_fid)

    # ---- the input was verified BY HASH, not by shape --------------------
    ii = r.get("input_identity") or {}
    out["input_identity_verified"] = bool(ii.get("verified"))
    if not ii.get("verified"):
        out["disagreements"].append(
            "the input was not verified by hash against provenance.json")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark-index", default="out/benchmark/index.json")
    ap.add_argument("--indices", nargs="+",
                    default=["out/release/index.json",
                             "out/release/expansion_index.json"])
    ap.add_argument("--emit-manifest", default=None,
                    help="write the census-batch TSV and stop")
    ap.add_argument("--reports", default=None,
                    help="directory of <segment>.json census reports")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    rows = surfaces(Path(a.benchmark_index), [Path(p) for p in a.indices])

    if a.emit_manifest:
        p = Path(a.emit_manifest)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w") as fh:
            for r in rows:
                fh.write(f"{r['segment']}\t"
                         f"{Path(r['reference']).resolve()}\n")
        print(f"{len(rows)} reference surfaces -> {p}")
        return 0

    if not a.reports:
        print("REFUSING: give --reports (or --emit-manifest first)")
        return 2
    reports = Path(a.reports)

    results, missing = [], []
    for i, row in enumerate(rows, 1):
        rep = reports / f"{row['segment']}.json"
        bnd = reports / f"{row['segment']}.binding.json"
        if not bnd.is_file():
            results.append({"segment": row["segment"], "not_compared": [],
                            "disagreements": [f"no census binding at {bnd}"]})
            print(f"[{i:3d}] {row['segment']:52s} NO BINDING")
            continue
        if not rep.is_file():
            missing.append(row["segment"])
            results.append({"segment": row["segment"],
                            "disagreements": [f"no census report at {rep}"]})
            print(f"[{i:3d}] {row['segment']:52s} NO REPORT")
            continue
        res = compare(row, rep, bnd)
        results.append(res)
        mark = "agrees" if not res["disagreements"] else \
            f"DISAGREES ({len(res['disagreements'])})"
        print(f"[{i:3d}] {row['segment']:52s} {mark}", flush=True)

    bad = [r for r in results if r["disagreements"]]
    verdict = "PASS" if not bad else "FAIL"

    # How many traces actually EXERCISED each comparison. A test that
    # reports "274/274 agree" while a field was absent on 90 of them is
    # overstating itself; the certificates of traces that needed no cut
    # carry no fragmentation block, so there is nothing to compare.
    fields = ("clean", "retained_fraction", "removed_cells",
              "components_input", "R_main_area_weighted",
              "coordinate_fidelity")
    coverage = {}
    for k in fields:
        n = sum(1 for r in results
                if (r.get("evaluator") or {}).get(k) is not None
                and (r.get("certificate") or {}).get(k) is not None)
        coverage[k] = {"compared": n, "not_compared": len(results) - n}
    out = {
        "test": "PREREG-TOPOLOGY-BENCHMARK B7",
        "rule": ("the benchmark is not published unless the evaluator "
                 "reproduces every packaged trace's certificate"),
        "census": ("the OFFICIAL vc_tifxyz_selfcross, swept in one "
                   "container; our own engine agreeing with our own "
                   "certificates would prove nothing about "
                   "interoperability"),
        "tolerances": {"retained_fraction_relative": REL_TOL,
                       "R_main_absolute": ABS_TOL,
                       "removed_cells": "exact",
                       "clean": "exact"},
        "n_traces": len(results),
        "coverage": coverage,
        "coverage_note": (
            "fragmentation is compared only on traces that were actually "
            "cut: a certificate for a trace that needed no excision "
            "carries no component_recovery block, so there is nothing to "
            "disagree with. This is stated rather than left to look like "
            "full coverage."),
        "n_agree": len(results) - len(bad),
        "n_disagree": len(bad),
        "n_missing_report": len(missing),
        "verdict": verdict,
        "disagreements": bad[:50],
        "results": results,
    }
    if a.out:
        Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"\nB7: {out['n_agree']}/{out['n_traces']} agree with their "
          f"certificate -- {verdict}")
    for r in bad[:10]:
        print(f"  {r['segment']}: {r['disagreements'][0]}")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
