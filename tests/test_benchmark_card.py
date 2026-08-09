"""Regressions for the benchmark card's refusals.

The card previously printed five different kinds of figure it had not
established, under a heading promising it never does. Each refusal that
replaced them was demonstrated once, by hand, against a mutated copy of
the real records — which means each survives only as long as someone
remembers to redo it. These are the same mutations, run every time.

Every case builds a complete copy of the card's inputs, breaks exactly
one thing, and requires a refusal naming that thing. The control proves
the fixture itself still produces a card.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CARD = ROOT / "bench" / "benchmark_card.py"
BENCH = ROOT / "out" / "benchmark"
RECORDS = ("index.json", "B7-acceptance.json", "selfcheck.json",
           "witness-binding.json", "negative-controls.json", "locators.json")
INDICES = (ROOT / "out/release/index.json",
           ROOT / "out/release/expansion_index.json")
RECON = ROOT / "out" / "reconciliation.json"

pytestmark = pytest.mark.skipif(
    not (all((BENCH / r).is_file() for r in RECORDS)
         and all(p.is_file() for p in INDICES) and RECON.is_file()),
    reason="benchmark records are not present in this checkout")


@pytest.fixture
def package(tmp_path):
    d = tmp_path / "pkg"
    (d / "release").mkdir(parents=True)
    for r in RECORDS:
        shutil.copy(BENCH / r, d / r)
    for p in INDICES:
        shutil.copy(p, d / "release" / p.name)
    shutil.copy(RECON, d / "reconciliation.json")
    return d


def edit(d: Path, name: str, fn):
    p = d / name
    j = json.loads(p.read_text())
    fn(j)
    p.write_text(json.dumps(j))


def run(d: Path):
    r = subprocess.run(
        [sys.executable, str(CARD), "--benchmark", str(d),
         "--reconciliation", str(d / "reconciliation.json"),
         "--indices", str(d / "release/index.json"),
         str(d / "release/expansion_index.json"),
         "--out", str(d / "CARD.md")],
        capture_output=True, text=True, cwd=ROOT)
    return r.returncode, (r.stdout + r.stderr)


def refusal(d: Path) -> str:
    rc, out = run(d)
    assert rc != 0, f"the card did NOT refuse; it wrote:\n{out}"
    assert "REFUSING" in out, out
    return out


# --- the control ---------------------------------------------------------

def test_the_real_records_produce_a_card(package):
    rc, out = run(package)
    assert rc == 0, out
    card = (package / "CARD.md").read_text()
    assert "274 scorable cases" in card
    assert "7/7 negative rejected; 1/1 positive accepted" in card
    # The claim that was published wrong. It is derived now, so this is a
    # regression against the specific sentence, not just the machinery.
    assert "**6 of 274**" in card
    assert "20251217234605-w2_20251217234605189" in card


# --- a missing verification record --------------------------------------

@pytest.mark.parametrize("record", RECORDS[1:])
def test_a_missing_record_is_not_an_unknown_number(package, record):
    (package / record).unlink()
    assert "REFUSING" in refusal(package)


def test_a_missing_reconciliation_refuses(package):
    (package / "reconciliation.json").unlink()
    assert "Population figures come from the reconciliation" in refusal(package)


# --- classification the card must not guess at --------------------------

def test_unknown_base_kind(package):
    edit(package, "release/index.json",
         lambda j: j["segments"][next(iter(j["segments"]))]
         .update(base_kind="repaired"))
    assert "will not classify an unknown base as repaired" in refusal(package)


def test_unknown_disposition(package):
    def m(j):
        next(t for t in j["traces"] if t.get("scored"))["disposition"] = "cut"
    edit(package, "index.json", m)
    assert "not one of ['transformed', 'already_clean']" in refusal(package)


# --- the core gate, which was asserted rather than read -----------------

def test_missing_core_gate_verdict(package):
    edit(package, "B7-acceptance.json",
         lambda j: j["results"][7]["evaluator"].pop("core_gate_pass"))
    assert "not a boolean" in refusal(package)


def test_core_gate_verdict_as_a_string(package):
    edit(package, "B7-acceptance.json",
         lambda j: j["results"][3]["evaluator"]
         .__setitem__("core_gate_pass", "false"))
    assert "not a boolean" in refusal(package)


def test_evaluator_and_certificate_disagree(package):
    edit(package, "B7-acceptance.json",
         lambda j: j["results"][3]["certificate"]
         .__setitem__("core_gate_pass", False))
    assert "core-gate verdict disagrees" in refusal(package)


def test_card_naming_different_failures_from_the_reconciliation(package):
    def m(j):
        j["claims"]["core_gate_failures"] = \
            j["claims"]["core_gate_failures"][:-1] + ["another-segment"]
    edit(package, "reconciliation.json", m)
    assert "name different core-gate failures" in refusal(package)


# --- a partial verification -----------------------------------------------

def test_acceptance_record_covering_part_of_the_population(package):
    edit(package, "B7-acceptance.json",
         lambda j: j.__setitem__("results", j["results"][:200]))
    assert "does not cover the packaged population exactly" in refusal(package)


def test_duplicate_acceptance_row(package):
    edit(package, "B7-acceptance.json",
         lambda j: j["results"].append(j["results"][0]))
    assert "lists a segment twice" in refusal(package)


@pytest.mark.parametrize("record,key,value", [
    ("B7-acceptance.json", "n_agree", 273),
    ("selfcheck.json", "n_ok", 12),
    ("witness-binding.json", "n_identical", 0),
])
def test_a_partial_result_is_not_a_verification(package, record, key, value):
    edit(package, record, lambda j: j.__setitem__(key, value))
    assert "reports" in refusal(package)


def test_locator_coverage_short_of_the_population(package):
    edit(package, "locators.json",
         lambda j: j["coverage"].__setitem__("inputs_located", 200))
    assert "covers 200 traces" in refusal(package)


def test_acceptance_disagreements_must_be_zero(package):
    edit(package, "B7-acceptance.json", lambda j: j.__setitem__("n_disagree", 2))
    assert "n_disagree" in refusal(package)


def test_a_failing_verdict_is_not_narrated(package):
    edit(package, "negative-controls.json",
         lambda j: j.__setitem__("verdict", "FAIL"))
    assert "does not narrate a failing one" in refusal(package)


def test_an_empty_control_suite_is_not_a_passing_one(package):
    """A forged record declaring 0 of each control and a PASS verdict
    satisfies numerator == denominator."""
    def m(j):
        j.update(n_negative=0, n_negative_rejected=0,
                 n_positive=0, n_positive_accepted=0, verdict="PASS")
    edit(package, "negative-controls.json", m)
    assert "the suite has" in refusal(package)


def test_one_control_short(package):
    edit(package, "negative-controls.json",
         lambda j: j.__setitem__("n_negative_rejected", 6))
    assert "negative controls rejected is 6/7" in refusal(package)


# --- counts and fractions -------------------------------------------------

def test_negative_contact_count(package):
    def m(j):
        next(t for t in j["traces"] if t.get("scored"))["n_contacts"] = -5
    edit(package, "index.json", m)
    assert "not a non-negative integer" in refusal(package)


def test_retention_outside_the_unit_interval(package):
    def m(j):
        next(t for t in j["traces"]
             if t.get("scored"))["retained_fraction_certificate"] = 1.7
    edit(package, "index.json", m)
    assert "not a finite fraction in [0, 1]" in refusal(package)


def test_summary_total_disagreeing_with_the_rows(package):
    edit(package, "index.json",
         lambda j: j["summary"].__setitem__("n_contacts_total", 999))
    assert "summed per-case contacts" in refusal(package)


def test_population_disagreeing_with_the_reconciliation(package):
    edit(package, "reconciliation.json",
         lambda j: j["claims"].__setitem__("n_censusable", 273))
    assert "the reconciliation says 273" in refusal(package)


# --- exclusions must carry the confirmation the heading claims ----------

def test_expansion_exclusion_without_a_confirmed_decline(package):
    def m(j):
        t = next(t for t in j["traces"]
                 if t.get("excluded_reason") and t.get("inventory") == "expansion")
        t["excluded_evidence"]["decline_confirmed"] = False
    edit(package, "index.json", m)
    assert "no confirmed census decline" in refusal(package)


def test_pinned_exclusion_that_retains_quads(package):
    def m(j):
        t = next(t for t in j["traces"]
                 if t.get("excluded_reason") and t.get("inventory") == "pinned")
        t["excluded_evidence"]["n_retained_quads"] = 5
    edit(package, "index.json", m)
    assert "retains 5 quads" in refusal(package)


def test_exclusion_with_no_valid_cell_count(package):
    def m(j):
        t = next(t for t in j["traces"] if t.get("excluded_reason"))
        t["excluded_evidence"] = {}
    edit(package, "index.json", m)
    assert "no valid-cell count" in refusal(package)
