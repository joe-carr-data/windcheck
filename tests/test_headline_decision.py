"""Tests for the round-28 Q4 headline decision rule (bench/headline_decision.py).

Every corpus here is FABRICATED under tmp_path: nothing depends on the real
185-segment run, so these pass on a checkout with no out/ at all.  The point
of the file is the arithmetic that decides the headline -- unique-geometry
weighting, the eight gates, and above all the rule that a failing 98.997%
must never be displayed as 99%.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "bench"))

import headline_decision as hd                                  # noqa: E402


# ------------------------------------------------------------------ helpers

def cert_body(segment: str, *, disposition: str = "transformed",
              a_orig: float = 100000.0, a_removed: float = 0.0,
              d0: int = 0, d1: int = 0, wall: float = 120.0,
              core_gate: bool = True, empty: bool = False,
              reload_ok: bool = True, area_blocks: bool = True,
              retained_fraction: bool = True,
              input_mesh: str | None = None,
              original_mesh: str | None = None,
              input_area_canonical: float | None = None,
              maxedge: float | None = None) -> dict:
    """One certificate.

    `area_blocks=False` reproduces the shape of a REAL already-clean record:
    no cut was made, so the certificate carries neither `headline_area` nor
    the operational `area` block -- only `input_area_canonical`, the canonical
    area of the input mesh's retained quads.
    """
    retained = 1.0 - a_removed / a_orig if a_orig else None
    body = {
        "record_kind": "excision certificate",
        "segment": segment,
        "terminal_disposition": disposition,
        "status": "transverse_clean" if d0 == d1 == 0 else "residual_transverse",
        "policy_hash": "deadbeef" * 8,
        "policy_version": "round28-frozen",
        "code_commit": "0" * 40,
        "census_after": {"d0": {"transverse": d0, "triangles": 10},
                         "d1": {"transverse": d1, "triangles": 10}},
        "output_transverse_total": d0 + d1,
        "reload_checks": {"grid_shape_equal": True, "dtype_equal": True,
                          "valid_mask_equals_intended": reload_ok,
                          "changes_only_valid_to_invalid": True},
        "operational_retained_fraction": retained,
        "headline_retained_fraction": retained,
        "core_gate_pass": core_gate,
        "component_recovery": {"core_gate": {"core_gate_pass": core_gate,
                                             "n_core_components": 1,
                                             "min_R_main_core":
                                                 0.99 if core_gate else 0.42}},
        "emptiness_guard": {"clean_by_emptiness": empty},
        "wall_seconds": wall,
        "headline_area": {"A_original_canonical": a_orig,
                          "A_removed_priced_on_original": a_removed},
        "area": {"canonical": {"A_input": a_orig,
                               "A_excised": a_removed,
                               "retained_fraction": retained}},
        "output_mesh": f"out/excised/corpus/{segment}.tifxyz",
        "output_mesh_hashes": {"x.tif": "a" * 64},
    }
    if not retained_fraction:
        body.pop("headline_retained_fraction")
        body.pop("operational_retained_fraction")
    if not area_blocks:
        body.pop("headline_area")
        body.pop("area")
    if input_area_canonical is not None:
        body["input_area_canonical"] = input_area_canonical
    if input_mesh is not None:
        body["input_mesh"] = input_mesh
    if original_mesh is not None:
        body["original_mesh"] = original_mesh
    if maxedge is not None:
        body["census_params"] = {"cell": 40.0, "exclude": 1,
                                 "maxedge": maxedge, "threads": 0,
                                 "touch_tol": 1e-3, "diagonals": [0, 1]}
    return body


def build_corpus(tmp_path: Path, specs: list[dict], *, name: str = "corpus",
                 record_hashes: bool = True) -> dict:
    """Write a base manifest, certificates and a summary index.

    Each spec is {"segment": s, "duplicate_of": t | None, **cert kwargs}.
    Returns the CLI paths plus the on-disk certificate paths.
    """
    root = tmp_path / name
    certs = root / "certs"
    certs.mkdir(parents=True, exist_ok=True)

    entries, rows, cert_paths = [], [], {}
    for spec in specs:
        spec = dict(spec)
        seg = spec.pop("segment")
        dup = spec.pop("duplicate_of", None)
        manifest_extra = spec.pop("manifest_extra", None) or {}
        entries.append({
            "segment": seg,
            "is_canonical": dup is None,
            "duplicate_of": dup,
            "geometry_key": hashlib.sha256((dup or seg).encode()).hexdigest(),
            "base_kind": "original",
            "original_mesh": f"data/synthetic/{seg}.tifxyz",
            **manifest_extra,
        })
        body = cert_body(seg, **spec)
        p = certs / f"{seg}{hd.CERT_SUFFIX}"
        p.write_text(json.dumps(body, indent=1))
        cert_paths[seg] = p
        row = {"segment": seg, "certificate": str(p),
               "terminal_disposition": body["terminal_disposition"]}
        if record_hashes:
            row["certificate_sha256"] = hd.sha256_file(p)
        rows.append(row)

    manifest = root / "corpus_bases.json"
    manifest.write_text(json.dumps({"schema": "corpus_bases/v1",
                                    "n_segments": len(entries),
                                    "entries": entries}, indent=1))
    summary = root / "corpus_summary.jsonl"
    summary.write_text("".join(json.dumps(r) + "\n" for r in rows))
    return {"manifest": manifest, "summary": summary, "certs": certs,
            "cert_paths": cert_paths, "root": root}


def run(corpus: dict, capsys, *, total: int | None = None,
        extra: list[str] | None = None) -> tuple[int, str]:
    argv = ["--base-manifest", str(corpus["manifest"]),
            "--summary", str(corpus["summary"]),
            "--certificates-dir", str(corpus["certs"]),
            # Never let a test touch the repository's own area cache.
            "--area-cache", str(corpus["root"] / "original_areas.json"),
            "--expected-total",
            str(total if total is not None
                else len(json.loads(corpus["manifest"].read_text())["entries"]))]
    argv += extra or []
    code = hd.main(argv)
    return code, capsys.readouterr().out


def conds(out: str) -> dict[int, str]:
    """{condition number: 'PASS'|'FAIL'} scraped from the printed report."""
    got = {}
    for line in out.splitlines():
        if line.startswith("COND "):
            parts = line.split()
            got[int(parts[1])] = parts[2].strip("[]")
    return got


def offender_block(out: str, n: int) -> str:
    """Everything printed under COND n, up to the next condition."""
    lines, keep, buf = out.splitlines(), False, []
    for line in lines:
        if line.startswith(f"COND {n} "):
            keep = True
        elif line.startswith("COND "):
            keep = False
        if keep:
            buf.append(line)
    return "\n".join(buf)


# A clean three-segment corpus: 300 of 300000 removed -> 99.900% overall,
# worst segment 99.700%.
CLEAN_SPECS = [
    {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 300.0},
    {"segment": "20240101000002", "a_orig": 100000.0, "a_removed": 0.0,
     "disposition": "already_clean"},
    {"segment": "20240101000003", "a_orig": 100000.0, "a_removed": 0.0},
]


# ------------------------------------------------------- all eight passing

def test_all_conditions_pass_gives_strong_sentence(tmp_path, capsys):
    corpus = build_corpus(tmp_path, CLEAN_SPECS)
    code, out = run(corpus, capsys)
    assert conds(out) == {n: "PASS" for n in range(1, 9)}
    assert code == 0
    assert "VERDICT: STRONG" in out
    assert (
        "windcheck produced a reload-verified, "
        "transverse-self-intersection-free tifxyz version of all 3 pinned "
        "trace artifacts from five scrolls, retaining 99.900% of original "
        "represented surface area overall and at least 99.700% per trace."
    ) in out
    # The fallback wording must not leak into a strong result.
    assert "were audited" not in out


def test_strong_run_emits_json_verdict(tmp_path, capsys):
    corpus = build_corpus(tmp_path, CLEAN_SPECS)
    out_json = tmp_path / "verdict" / "headline.json"
    code, _ = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    assert code == 0
    assert doc["verdict"] == "STRONG"
    assert doc["n_conditions_passed"] == 8
    assert [c["n"] for c in doc["conditions"]] == list(range(1, 9))
    assert all(c["text"] for c in doc["conditions"])
    assert doc["retention"]["X_percent"] == "99.900"
    assert doc["retention"]["Y_percent"] == "99.700"
    assert doc["counts"]["transverse_clean_artifacts_N"] == 3
    assert doc["sentence"] == doc["sentence"].strip()


def test_verbatim_condition_text_matches_decisions_entry(tmp_path, capsys):
    corpus = build_corpus(tmp_path, CLEAN_SPECS)
    _, out = run(corpus, capsys)
    for fragment in (
            "every unique censusable mesh emits a reload-verified aggregate "
            "with transverse 0/0 both diagonals",
            "have a terminal disposition (transformed | already-clean | "
            "duplicate alias | explicitly triangle-empty/invalid)",
            "unique-geometry-weighted corpus retention >= 99.000%",
            "every claimed-clean segment >= 95.000% of ORIGINAL canonical area",
            "every component in every segment's input 99.9%-area core has "
            "R_main >= 0.90",
            "no segment clean by empty/near-empty output",
            "every production artifact completes within ten minutes under the "
            "frozen policy",
            "no missing certificate, failed reload, unresolved transverse "
            "contact, timeout or unaccounted hash mismatch"):
        assert fragment in out


# ------------------------------- THE REQUIRED FAILING CASE: 98.997%

def test_98997_percent_fails_condition_three_and_is_never_shown_as_99(
        tmp_path, capsys):
    """A corpus retaining exactly 98.997% must FAIL, and must never print 99%."""
    corpus = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 50000.0, "a_removed": 1003.0},
        {"segment": "20240101000002", "a_orig": 50000.0, "a_removed": 0.0},
    ])
    code, out = run(corpus, capsys)
    assert conds(out)[3] == "FAIL"
    assert code == 1
    assert "VERDICT: FALLBACK" in out
    assert "98.997" in out
    # "99.000%" occurs exactly twice as the GATE (the frozen condition text and
    # the evidence line that quotes it).  Outside those, no measured figure may
    # be shown as 99% -- that is the whole point of the rule.
    measured = out.replace(
        "unique-geometry-weighted corpus retention >= 99.000%", "<COND3>")
    measured = measured.replace("gate 99.000%", "<GATE>")
    assert "99.0" not in measured
    assert "99%" not in measured
    assert ("Of 2 pinned trace artifacts, 2 were censusable. All 2 now have "
            "reload-verified tifxyz outputs with zero non-adjacent transverse "
            "contacts under both canonical triangulations: 2 were transformed "
            "and 0 required no change. Zero triangle-empty or invalid inputs "
            "have explicit terminal records.") in out
    # The shortfall belongs to the AREA QUALIFICATION, printed separately --
    # it is never folded into the sentence as a clause.
    sentence = out.split("SENTENCE:\n", 1)[1].split("\n", 1)[0]
    assert "98.997" not in sentence
    assert "98.997" in out.split("QUALIFICATIONS", 1)[1]
    # The other seven gates are untouched by a shortfall in area.
    assert [conds(out)[n] for n in (1, 2, 4, 5, 6, 7, 8)] == ["PASS"] * 7


def test_98997_percent_json_carries_the_untouched_float(tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 50000.0, "a_removed": 1003.0},
        {"segment": "20240101000002", "a_orig": 50000.0, "a_removed": 0.0},
    ])
    out_json = tmp_path / "headline.json"
    code, _ = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    assert code == 1
    assert doc["verdict"] == "FALLBACK"
    assert doc["retention"]["X_percent"] == "98.997"
    assert doc["retention"]["unique_geometry_weighted"] < 0.99
    # The retention figure lives in the qualifications, not the sentence.
    assert "98.997" not in doc["sentence"]
    area = next(q for q in doc["qualifications"]
                if q["kind"] == "area_retention")
    assert "98.997" in area["text"]
    assert "99.0" not in area["text"].replace("gate 99.000%", "<GATE>")


# ------------------------------------------------ duplicate aliases counted once

def test_duplicate_alias_does_not_change_unique_retention(tmp_path, capsys):
    base_specs = [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 1000.0},
        {"segment": "20240101000002", "a_orig": 100000.0, "a_removed": 0.0},
    ]
    alias_specs = base_specs + [
        {"segment": "20240101000003", "duplicate_of": "20240101000001",
         "disposition": "duplicate_alias", "a_orig": 100000.0,
         "a_removed": 1000.0},
    ]
    j_plain = tmp_path / "plain.json"
    j_alias = tmp_path / "alias.json"
    run(build_corpus(tmp_path, base_specs, name="plain"), capsys,
        extra=["--json", str(j_plain)])
    run(build_corpus(tmp_path, alias_specs, name="alias"), capsys,
        extra=["--json", str(j_alias)])
    plain, alias = json.loads(j_plain.read_text()), json.loads(j_alias.read_text())

    assert plain["retention"]["unique_geometry_weighted"] == \
        alias["retention"]["unique_geometry_weighted"]
    assert alias["retention"]["X_percent"] == plain["retention"]["X_percent"] \
        == "99.500"
    # Denominators identical: the duplicated geometry is priced exactly once.
    assert alias["retention"]["unique_A_original_total"] == 200000.0
    assert alias["retention"]["unique_segments_counted"] == 2
    # ... while the artifact-count figure DOES see all three artifacts, so the
    # equality above is a real exclusion and not a coincidence.
    assert alias["retention"]["artifact_count_weighted"] < \
        alias["retention"]["unique_geometry_weighted"]
    assert alias["counts"]["duplicate_aliases"] == 1
    assert alias["counts"]["transverse_clean_artifacts_N"] == 3
    assert alias["verdict"] == "STRONG"


# ------------------------------------------------------- individual gate trips

def test_segment_below_95_percent_fails_condition_four_by_name(tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 5100.0},
        {"segment": "20240101000002", "a_orig": 100000.0, "a_removed": 0.0},
    ])
    code, out = run(corpus, capsys)
    assert conds(out)[4] == "FAIL"
    assert code == 1
    block = offender_block(out, 4)
    assert "20240101000001" in block
    assert "94.900" in block
    assert "20240101000002" not in block


def test_wall_seconds_601_fails_condition_seven(tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        dict(CLEAN_SPECS[0]), dict(CLEAN_SPECS[1]),
        {"segment": "20240101000004", "a_orig": 100000.0, "a_removed": 0.0,
         "wall": 601.0},
    ])
    code, out = run(corpus, capsys)
    assert conds(out)[7] == "FAIL"
    assert conds(out)[8] == "FAIL"          # a timeout is also an accounting failure
    assert code == 1
    assert "20240101000004" in offender_block(out, 7)


def test_wall_seconds_exactly_600_still_passes(tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 0.0,
         "wall": 600.0}])
    code, out = run(corpus, capsys)
    assert conds(out)[7] == "PASS"
    assert code == 0


def test_missing_certificate_fails_conditions_two_and_eight(tmp_path, capsys):
    corpus = build_corpus(tmp_path, CLEAN_SPECS)
    corpus["cert_paths"]["20240101000002"].unlink()
    code, out = run(corpus, capsys)
    got = conds(out)
    assert got[2] == "FAIL" and got[8] == "FAIL"
    assert code == 1
    assert "20240101000002" in offender_block(out, 2)
    assert "missing certificate" in offender_block(out, 8)


def test_certificate_hash_mismatch_fails_condition_eight(tmp_path, capsys):
    corpus = build_corpus(tmp_path, CLEAN_SPECS)
    p = corpus["cert_paths"]["20240101000003"]
    body = json.loads(p.read_text())
    body["note"] = "edited after the summary was written"
    p.write_text(json.dumps(body, indent=1))
    code, out = run(corpus, capsys)
    assert conds(out)[8] == "FAIL"
    assert code == 1
    assert "hash mismatch" in offender_block(out, 8)
    # ... and the check is exactly what --no-strict-artifacts turns off.
    code2, out2 = run(corpus, capsys, extra=["--no-strict-artifacts"])
    assert conds(out2)[8] == "PASS" and code2 == 0


def test_core_gate_false_fails_condition_five(tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        dict(CLEAN_SPECS[0]),
        {"segment": "20240101000005", "a_orig": 100000.0, "a_removed": 0.0,
         "core_gate": False},
    ])
    code, out = run(corpus, capsys)
    assert conds(out)[5] == "FAIL"
    assert code == 1
    assert "20240101000005" in offender_block(out, 5)


def test_clean_by_emptiness_fails_condition_six(tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        dict(CLEAN_SPECS[0]),
        {"segment": "20240101000006", "a_orig": 100000.0, "a_removed": 0.0,
         "empty": True},
    ])
    code, out = run(corpus, capsys)
    assert conds(out)[6] == "FAIL"
    assert code == 1
    assert "20240101000006" in offender_block(out, 6)


def test_residual_transverse_on_d1_fails_condition_one(tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        dict(CLEAN_SPECS[0]),
        {"segment": "20240101000007", "a_orig": 100000.0, "a_removed": 0.0,
         "d1": 1},
    ])
    code, out = run(corpus, capsys)
    got = conds(out)
    assert got[1] == "FAIL"
    assert got[8] == "FAIL"                 # unresolved transverse contact
    assert code == 1
    block = offender_block(out, 1)
    assert "20240101000007" in block and "d1 transverse 1" in block


def test_failed_reload_check_fails_conditions_one_and_eight(tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        dict(CLEAN_SPECS[0]),
        {"segment": "20240101000008", "a_orig": 100000.0, "a_removed": 0.0,
         "reload_ok": False},
    ])
    code, out = run(corpus, capsys)
    got = conds(out)
    assert got[1] == "FAIL" and got[8] == "FAIL"
    assert code == 1
    assert "valid_mask_equals_intended" in offender_block(out, 1)


def test_non_terminal_dispositions_fail_condition_two(tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        dict(CLEAN_SPECS[0]),
        {"segment": "20240101000009", "disposition": "error"},
        {"segment": "20240101000010", "disposition": "residual_transverse",
         "d0": 3},
        {"segment": "20240101000011", "disposition": "not_censusable"},
    ])
    code, out = run(corpus, capsys)
    assert conds(out)[2] == "FAIL"
    assert code == 1
    block = offender_block(out, 2)
    for seg in ("20240101000009", "20240101000010", "20240101000011"):
        assert seg in block
    assert "20240101000001" not in block


def test_triangle_empty_invalid_is_terminal_and_not_an_area_gap(
        tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        dict(CLEAN_SPECS[0]),
        {"segment": "20240101000012", "disposition": "triangle_empty_invalid",
         "a_orig": 0.0, "a_removed": 0.0},
    ])
    out_json = tmp_path / "v.json"
    code, out = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    assert conds(out)[2] == "PASS"
    assert conds(out)[3] == "PASS"
    assert doc["retention"]["unique_missing_area"] == []
    assert doc["counts"]["transverse_clean_artifacts_N"] == 1   # not the empty one
    # An explicitly triangle-empty/invalid record IS terminal under Q4, so all
    # eight conditions still pass and the verdict stays STRONG.
    assert code == 0
    assert doc["verdict"] == "STRONG"


def test_roster_shorter_than_expected_total_fails_condition_two(
        tmp_path, capsys):
    corpus = build_corpus(tmp_path, CLEAN_SPECS)
    code, out = run(corpus, capsys, total=185)
    assert conds(out)[2] == "FAIL"
    assert code == 1
    assert "roster has 3 entries, expected 185" in offender_block(out, 2)
    assert "all 185 entries have a terminal disposition" in out


def test_offender_list_is_capped_at_twenty_with_a_count(tmp_path, capsys):
    specs = [{"segment": f"202401010{i:05d}", "disposition": "error"}
             for i in range(25)]
    code, out = run(build_corpus(tmp_path, specs), capsys)
    block = offender_block(out, 2)
    assert block.count("      - ") == hd.MAX_LISTED
    assert "... and 5 more" in block
    assert code == 1


# ------------------------------------------------------------- rounding rules

def test_truncation_never_rounds_up():
    assert hd.truncate_pct(0.98999) == "98.999"
    assert hd.truncate_pct(0.98997) == "98.997"
    assert hd.truncate_pct(0.9899999) == "98.999"
    assert hd.truncate_pct(0.99) == "99.000"
    assert hd.truncate_pct(1.0) == "100.000"
    assert hd.truncate_pct(None) == "n/a"


def test_failing_retention_is_floored_below_99(tmp_path):
    # A value a hair under the gate whose shortest repr would print as 99.000
    # is still a FAILING value, so it must not be displayed as 99%.
    hair = math.nextafter(0.99, 0.0)
    assert hair < hd.RETENTION_GATE
    assert hd.fmt_retention(hair, passed=False) == "98.999"
    assert hd.fmt_retention(0.999, passed=True) == "99.900"


def test_condition_three_compares_the_exact_float(tmp_path, capsys):
    # 0.99 exactly passes; one ulp below it fails.
    passing = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 1000.0}],
        name="edge_pass")
    code, out = run(passing, capsys)
    assert conds(out)[3] == "PASS" and code == 0
    assert "99.000%" in out

    failing = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 1000.01}],
        name="edge_fail")
    code2, out2 = run(failing, capsys)
    assert conds(out2)[3] == "FAIL" and code2 == 1
    assert "98.999" in out2


# ------------------------------------------------------------- input handling

def test_missing_base_manifest_exits_two(tmp_path, capsys):
    code = hd.main(["--base-manifest", str(tmp_path / "nope.json"),
                    "--summary", str(tmp_path / "nope.jsonl"),
                    "--certificates-dir", str(tmp_path)])
    assert code == 2


def test_works_without_a_summary_file(tmp_path, capsys):
    """Certificates are found at the default path when no summary exists."""
    corpus = build_corpus(tmp_path, CLEAN_SPECS)
    corpus["summary"].unlink()
    code, out = run(corpus, capsys)
    assert conds(out) == {n: "PASS" for n in range(1, 9)}
    assert code == 0


def test_summary_row_without_a_base_entry_fails_two_and_eight(
        tmp_path, capsys):
    corpus = build_corpus(tmp_path, CLEAN_SPECS)
    with open(corpus["summary"], "a") as fh:
        fh.write(json.dumps({"segment": "29999999999999",
                             "terminal_disposition": "transformed"}) + "\n")
    code, out = run(corpus, capsys)
    got = conds(out)
    assert got[2] == "FAIL" and got[8] == "FAIL"
    assert code == 1
    assert "29999999999999" in offender_block(out, 8)


def test_both_summaries_are_reported(tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 1000.0},
        {"segment": "20240101000002", "a_orig": 100000.0, "a_removed": 0.0},
        {"segment": "20240101000003", "duplicate_of": "20240101000001",
         "disposition": "duplicate_alias", "a_orig": 100000.0,
         "a_removed": 1000.0},
    ])
    _, out = run(corpus, capsys)
    assert "unique-geometry-weighted retention : 99.500%" in out
    assert "artifact-count-weighted retention  : 99.333%" in out


# ------------------------------------ pricing the segments nobody had to cut
#
# An already-clean certificate records that NOTHING WAS DONE to the segment:
# no cut, no output mesh, and so no `headline_area` block.  Before this was
# fixed those segments fell out of the retention denominator entirely, and
# condition 3 -- which gates on the corpus-wide ratio -- was being evaluated
# over the cut segments alone.  That is a bookkeeping defect, not an area
# result: an unmodified segment retains exactly 1.0 of its original area.

def unpriced_clean_spec(segment: str, area: float) -> dict:
    """An already-clean spec shaped like the real record: no area blocks."""
    return {"segment": segment, "disposition": "already_clean",
            "area_blocks": False, "input_area_canonical": area,
            "input_mesh": f"data/synthetic/{segment}.tifxyz",
            "original_mesh": f"data/synthetic/{segment}.tifxyz",
            "maxedge": 60.0}


def test_unmodified_segment_is_priced_at_exactly_one(tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 1000.0},
        unpriced_clean_spec("20240101000002", 100000.0),
    ])
    out_json = tmp_path / "v.json"
    code, out = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    r = doc["retention"]

    # It is in the denominator, at its full original area, with nothing removed.
    assert r["unique_missing_area"] == []
    assert r["unique_segments_counted"] == 2
    assert r["unique_A_original_total"] == 200000.0
    assert r["unique_A_removed_total"] == 1000.0
    assert r["X_percent"] == "99.500"
    assert conds(out)[3] == "PASS" and code == 0

    priced = doc["area_pricing"]
    assert priced["n_derived"] == 1
    assert priced["failures"] == []
    assert priced["derived"][0]["segment"] == "20240101000002"
    assert priced["derived"][0]["A_original_canonical"] == 100000.0
    # No mesh had to be re-read: the input IS the original published mesh.
    assert priced["n_recomputed_from_mesh"] == 0
    assert "input_area_canonical" in priced["derived"][0]["source"]


def test_unpriced_clean_segments_no_longer_sink_condition_three(
        tmp_path, capsys):
    """The exact defect: the gate failing on bookkeeping, not on area.

    The cut segment alone retains 98.500%.  Counting the untouched segment --
    which lost nothing -- the corpus retains 99.250%, and 99.250% is what the
    gate is entitled to see.
    """
    specs = [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 1500.0},
        unpriced_clean_spec("20240101000002", 100000.0),
    ]
    fixed_code, fixed_out = run(build_corpus(tmp_path, specs, name="fixed"),
                                capsys)
    out_json = tmp_path / "fixed.json"
    run(build_corpus(tmp_path, specs, name="fixed2"), capsys,
        extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    assert conds(fixed_out)[3] == "PASS" and fixed_code == 0
    assert doc["retention"]["X_percent"] == "99.250"
    assert doc["retention"]["unique_segments_counted"] == 2

    # The cut segment on its own is the 98.500% that used to be reported as
    # the whole corpus, and it is still exactly what its own record says.
    assert doc["retention"]["min_headline_retained_fraction"] == 0.985


def test_condition_four_counts_the_unmodified_segments_at_one(
        tmp_path, capsys):
    """A certificate that states no retained fraction still means 1.0."""
    spec = unpriced_clean_spec("20240101000002", 100000.0)
    spec["retained_fraction"] = False
    corpus = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 1000.0},
        spec,
    ])
    out_json = tmp_path / "v.json"
    code, out = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    assert conds(out)[4] == "PASS" and code == 0
    # Two claimed-clean artifacts are checked, not one, and the untouched
    # segment is the best of them rather than an absentee.
    assert "2/2 claimed-clean artifacts" in offender_block(out, 4)
    assert doc["retention"]["min_headline_segment"] == "20240101000001"
    assert doc["retention"]["min_headline_retained_fraction"] == 0.99


# --------------------------------------- deriving an area from a real mesh

def write_planar_tifxyz(path: Path, nv: int, nu: int, step: float) -> float:
    """A flat (nv x nu) grid of squares; returns its exact canonical area."""
    import numpy as np
    import tifffile

    path.mkdir(parents=True, exist_ok=True)
    vv, uu = np.meshgrid(np.arange(nv), np.arange(nu), indexing="ij")
    planes = {"x": uu * step, "y": vv * step, "z": np.zeros((nv, nu))}
    for axis, plane in planes.items():
        tifffile.imwrite(path / f"{axis}.tif", plane.astype("float32"))
    return (nv - 1) * (nu - 1) * step * step


def test_area_is_recomputed_from_the_original_mesh_and_then_cached(
        tmp_path, capsys):
    """A displacement-repaired base: the original mesh must be re-read.

    The area is the one `bench/excise_shadow.headline_denominator` would
    build -- the canonical area of the ORIGINAL mesh's retained quads -- so a
    flat grid of 10um squares must come back as exactly its square count.
    """
    mesh = tmp_path / "meshes" / "original.tifxyz"
    exact = write_planar_tifxyz(mesh, 4, 5, 10.0)
    assert exact == 1200.0

    spec = unpriced_clean_spec("20240101000002", 999.0)   # repaired-base input
    spec["input_mesh"] = str(tmp_path / "meshes" / "repaired.tifxyz")
    spec["original_mesh"] = str(mesh)
    # The base is displacement-repaired, and that is established by CONTENT
    # -- the base manifest's hashes differ from the original's -- never by
    # the two meshes sitting at different paths.
    spec["manifest_extra"] = {"base_kind": "displacement_repaired",
                              "base_hashes": {"x": "d" * 64},
                              "original_hashes": {"x": "a" * 64}}
    corpus = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 100.0},
        spec,
    ])
    cache = corpus["root"] / "original_areas.json"
    out_json = tmp_path / "v.json"
    code, out = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())

    row = doc["area_pricing"]["derived"][0]
    assert row["segment"] == "20240101000002"
    assert row["A_original_canonical"] == exact
    assert row["n_original_retained_quads"] == 12
    assert row["grid_shape"] == [4, 5]
    assert doc["area_pricing"]["n_recomputed_from_mesh"] == 1
    # The input mesh's own area is NOT what got used: the headline denominator
    # is the ORIGINAL published mesh.
    assert row["A_original_canonical"] != 999.0
    assert doc["retention"]["unique_A_original_total"] == 100000.0 + exact
    assert code == 0

    # Cached, so a second run costs no mesh reads and reaches the same number.
    assert json.loads(cache.read_text())["areas"]["20240101000002"][
        "A_original_canonical"] == exact
    second_json = tmp_path / "v2.json"
    run(corpus, capsys, extra=["--json", str(second_json)])
    again = json.loads(second_json.read_text())
    assert again["area_pricing"]["n_recomputed_from_mesh"] == 0
    assert again["area_pricing"]["derived"][0]["A_original_canonical"] == exact
    assert again["retention"]["unique_geometry_weighted"] == \
        doc["retention"]["unique_geometry_weighted"]


def test_base_kind_is_decided_by_content_not_by_path(tmp_path, capsys):
    """Two DIFFERENT paths holding the SAME bytes is one mesh.

    The old rule compared `input_mesh` to `original_mesh` as paths, so a
    downloaded archive or a fresh workdir -- where the same mesh sits
    somewhere else -- read as a displacement-repaired base and sent the
    segment down the recompute route. The manifests decide now, so the
    shortcut fires and no mesh is re-read.
    """
    spec = unpriced_clean_spec("20240101000002", 4242.0)
    spec["input_mesh"] = str(tmp_path / "fresh_workdir" / "seg.tifxyz")
    spec["original_mesh"] = str(tmp_path / "download" / "seg.tifxyz")
    same = {"x": "a" * 64, "y": "b" * 64, "z": "c" * 64}
    spec["manifest_extra"] = {"base_kind": "original",
                              "base_hashes": dict(same),
                              "original_hashes": dict(same)}
    corpus = build_corpus(tmp_path, [spec])
    out_json = tmp_path / "v.json"
    run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    row = doc["area_pricing"]["derived"][0]
    assert row["A_original_canonical"] == 4242.0
    assert doc["area_pricing"]["n_recomputed_from_mesh"] == 0
    assert "the input IS the original published mesh" in row["source"]


def test_a_stale_cache_row_is_not_trusted(tmp_path, capsys):
    """A cached area belongs to one mesh with one set of hashes."""
    mesh = tmp_path / "meshes" / "original.tifxyz"
    exact = write_planar_tifxyz(mesh, 4, 5, 10.0)
    spec = unpriced_clean_spec("20240101000002", 999.0)
    spec["input_mesh"] = str(tmp_path / "meshes" / "repaired.tifxyz")
    spec["original_mesh"] = str(mesh)
    spec["manifest_extra"] = {"base_kind": "displacement_repaired",
                              "base_hashes": {"x": "d" * 64},
                              "original_hashes": {"x": "a" * 64}}
    corpus = build_corpus(tmp_path, [spec])
    cache = corpus["root"] / "original_areas.json"
    cache.write_text(json.dumps({"schema": hd.AREA_CACHE_SCHEMA, "areas": {
        "20240101000002": {"A_original_canonical": 1.0,
                           "original_mesh": str(mesh),
                           "original_hashes": {"x": "b" * 64},
                           "maxedge": 60.0}}}))
    out_json = tmp_path / "v.json"
    run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    assert doc["area_pricing"]["derived"][0]["A_original_canonical"] == exact
    assert doc["area_pricing"]["n_recomputed_from_mesh"] == 1
    assert json.loads(cache.read_text())["areas"]["20240101000002"][
        "original_hashes"] == {"x": "a" * 64}


def test_an_underivable_area_is_reported_as_a_gap_not_guessed(
        tmp_path, capsys):
    """No mesh, no cache, no shortcut: the segment stays visibly unpriced."""
    spec = unpriced_clean_spec("20240101000002", 100000.0)
    spec["input_mesh"] = "data/synthetic/repaired.tifxyz"
    spec["original_mesh"] = str(tmp_path / "absent.tifxyz")
    spec["manifest_extra"] = {"base_kind": "displacement_repaired",
                              "base_hashes": {"x": "d" * 64},
                              "original_hashes": {"x": "a" * 64}}
    corpus = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 100.0},
        spec,
    ])
    out_json = tmp_path / "v.json"
    code, out = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    assert conds(out)[3] == "FAIL" and code == 1
    assert doc["retention"]["unique_missing_area"] == ["20240101000002"]
    assert doc["area_pricing"]["derived"] == []
    assert doc["area_pricing"]["failures"]
    assert "20240101000002" in offender_block(out, 3)


# ------------------------------------------- the sentence and what qualifies it

def corpus_shaped_like_the_real_pass(tmp_path) -> dict:
    """154 transformed, 25 unmodified, 6 triangle-empty/invalid: 185 records.

    One transformed segment keeps only 94.678% of its original area, so
    condition 4 fails and the run lands on the fallback sentence -- which is
    the case the wording has to be right for.
    """
    specs = [{"segment": f"2024010100{i:04d}", "a_orig": 100000.0,
              "a_removed": 10.0} for i in range(154)]
    specs[0]["a_removed"] = 5322.0
    specs += [unpriced_clean_spec(f"2025010100{i:04d}", 100000.0)
              for i in range(25)]
    specs += [{"segment": f"2026010100{i:04d}",
               "disposition": "triangle_empty_invalid",
               "a_orig": 0.0, "a_removed": 0.0} for i in range(6)]
    return build_corpus(tmp_path, specs, name="full")


def test_fallback_sentence_states_the_census_and_nothing_else(
        tmp_path, capsys):
    corpus = corpus_shaped_like_the_real_pass(tmp_path)
    out_json = tmp_path / "v.json"
    code, out = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    assert code == 1 and doc["verdict"] == "FALLBACK"
    assert doc["sentence"] == (
        "Of 185 pinned trace artifacts, 179 were censusable. All 179 now "
        "have reload-verified tifxyz outputs with zero non-adjacent "
        "transverse contacts under both canonical triangulations: 154 were "
        "transformed and 25 required no change. Six triangle-empty or "
        "invalid inputs have explicit terminal records.")
    assert doc["sentence"] in out
    # The six exceptions are triangle-empty/invalid NON-AUDITS.  They are not
    # residual dispositions and they are not refusals, and the sentence must
    # not call them either.
    for wrong in ("residual", "refusal", "were audited"):
        assert wrong not in doc["sentence"]
    # Every number in it is derived, not spelled into the template.
    k = doc["counts"]
    assert (k["censusable"], k["transformed"], k["already_clean_unchanged"],
            k["not_censusable_triangle_empty_invalid"]) == (179, 154, 25, 6)


def test_qualifications_are_printed_below_and_apart_from_the_sentence(
        tmp_path, capsys):
    corpus = corpus_shaped_like_the_real_pass(tmp_path)
    out_json = tmp_path / "v.json"
    _, out = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())

    kinds = [q["kind"] for q in doc["qualifications"]]
    assert kinds == ["area_retention", "fragmentation"]

    # Neither qualification is inside the sentence, as a clause or otherwise.
    assert "retention" not in doc["sentence"]
    assert "area" not in doc["sentence"]
    assert "R_main" not in doc["sentence"]
    assert doc["retention"]["X_percent"] not in doc["sentence"]

    # ... and they are printed UNDER it, under their own heading.
    head, tail = out.split("SENTENCE:\n", 1)
    body, quals = tail.split(
        "QUALIFICATIONS (separate statements, NOT part of the sentence "
        "above):\n", 1)
    assert body.strip() == doc["sentence"]
    for q in doc["qualifications"]:
        assert q["text"] in quals
    assert doc["retention"]["X_percent"] in quals


def test_fragmentation_qualification_names_the_worst_core_recovery(
        tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        dict(CLEAN_SPECS[0]),
        {"segment": "20240101000005", "a_orig": 100000.0, "a_removed": 0.0,
         "core_gate": False},
    ])
    out_json = tmp_path / "v.json"
    _, out = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    frag = next(q for q in doc["qualifications"] if q["kind"] == "fragmentation")
    assert "1 segment does not meet the R_main >= 0.90 recovery gate" in \
        frag["text"]
    assert "0.42" in frag["text"] and "20240101000005" in frag["text"]
    assert frag["text"] in out


def test_a_clean_corpus_says_so_in_both_qualifications(tmp_path, capsys):
    corpus = build_corpus(tmp_path, CLEAN_SPECS)
    out_json = tmp_path / "v.json"
    code, _ = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    assert code == 0
    area, frag = doc["qualifications"]
    assert "99.900" in area["text"] and "99.700" in area["text"]
    assert "below the gate" not in area["text"]
    assert "every component in every segment's input 99.9%-area core meets" \
        in frag["text"]


def test_singular_wording_when_exactly_one_input_is_not_censusable(
        tmp_path, capsys):
    corpus = build_corpus(tmp_path, [
        {"segment": "20240101000001", "a_orig": 100000.0, "a_removed": 5100.0},
        {"segment": "20240101000012", "disposition": "triangle_empty_invalid",
         "a_orig": 0.0, "a_removed": 0.0},
    ])
    out_json = tmp_path / "v.json"
    _, _ = run(corpus, capsys, extra=["--json", str(out_json)])
    doc = json.loads(out_json.read_text())
    assert doc["sentence"].endswith(
        "One triangle-empty or invalid input has an explicit terminal record.")
