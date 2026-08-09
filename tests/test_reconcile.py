"""Negative tests for the reconciliation gate.

The gate's whole value is that it stops the build when two records
disagree about a segment. A gate that never refuses is decoration, so
every refusal it claims is exercised here against a mutated copy of the
real sources, with an untouched control proving the fixture itself is
sound.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bench"))
import reconcile  # noqa: E402

REAL = {k: ROOT / v for k, v in reconcile.SOURCES.items()}

pytestmark = pytest.mark.skipif(
    not all(p.is_file() for p in REAL.values()),
    reason="reconciliation sources are not present in this checkout")


@pytest.fixture(scope="module")
def sources():
    return {k: json.loads(p.read_text()) for k, p in REAL.items()}


def run(tmp_path, sources, mutate=None):
    data = copy.deepcopy(sources)
    if mutate:
        mutate(data)
    paths = {}
    for k, v in data.items():
        p = tmp_path / f"{k}.json"
        p.write_text(json.dumps(v))
        paths[k] = p
    return reconcile.build(paths)


def refusal(tmp_path, sources, mutate) -> str:
    with pytest.raises(reconcile.Contradiction) as e:
        run(tmp_path, sources, mutate)
    return str(e.value)


# --- the control ---------------------------------------------------------

def test_the_real_corpus_reconciles(tmp_path, sources):
    """Not a tautology: this is the assertion that the shipped records
    agree with each other, segment by segment."""
    rec = run(tmp_path, sources)
    c = rec["claims"]
    assert c["n_segments"] == len(rec["rows"])
    assert c["n_pinned"] + c["n_expansion"] == c["n_segments"]
    assert c["n_censusable"] + c["n_excluded_not_censusable"] == c["n_segments"]
    assert c["n_transformed"] + c["n_already_clean"] == c["n_censusable"]
    assert c["n_defect_bearing"] <= c["n_intrinsic_processed"]
    assert c["defect_bearing_operand"] == "published original trace"


def test_every_row_is_complete(tmp_path, sources):
    rec = run(tmp_path, sources)
    for r in rec["rows"]:
        assert r["inventory"] in ("pinned", "expansion")
        assert r["disposition"] in reconcile.DISPOSITIONS
        assert r["base_kind"] in reconcile.BASE_KINDS
        assert isinstance(r["censusable"], bool)
        assert isinstance(r["defect_bearing"], bool)
        if r["censusable"]:
            assert isinstance(r["core_gate_pass"], bool)
            assert r["official_clean"] is True


def test_sources_are_hashed(tmp_path, sources):
    rec = run(tmp_path, sources)
    assert set(rec["sources"]) == set(reconcile.SOURCES)
    for s in rec["sources"].values():
        assert len(s["sha256"]) == 64


# --- a missing row -------------------------------------------------------

def test_segment_missing_from_the_benchmark_index(tmp_path, sources):
    def m(d):
        d["benchmark_index"]["traces"].pop()
    assert "absent from the benchmark index" in refusal(tmp_path, sources, m)


def test_censusable_segment_with_no_acceptance_row(tmp_path, sources):
    def m(d):
        d["acceptance"]["results"].pop()
    assert "no acceptance row" in refusal(tmp_path, sources, m)


# --- a duplicate row -----------------------------------------------------

def test_duplicate_benchmark_row(tmp_path, sources):
    def m(d):
        d["benchmark_index"]["traces"].append(d["benchmark_index"]["traces"][0])
    assert "appears twice in the benchmark index" in refusal(tmp_path, sources, m)


def test_duplicate_acceptance_row(tmp_path, sources):
    def m(d):
        d["acceptance"]["results"].append(d["acceptance"]["results"][0])
    assert "appears twice in the acceptance record" in refusal(tmp_path, sources, m)


def test_segment_in_both_inventories(tmp_path, sources):
    def m(d):
        seg, ent = next(iter(d["release_pinned"]["segments"].items()))
        d["release_expansion"]["segments"][seg] = ent
    assert "more than one release index" in refusal(tmp_path, sources, m)


# --- a scope swap: the failure class this file exists for ----------------

def test_inventory_swap_is_caught(tmp_path, sources):
    """A pinned segment relabelled as expansion. This is the shape of the
    error that produced 'one trace fails the core gate'."""
    def m(d):
        t = next(t for t in d["benchmark_index"]["traces"]
                 if t["inventory"] == "pinned")
        t["inventory"] = "expansion"
    assert "says inventory='expansion'" in refusal(tmp_path, sources, m)


def test_disposition_disagreement_is_caught(tmp_path, sources):
    def m(d):
        t = next(t for t in d["benchmark_index"]["traces"]
                 if t["disposition"] == "transformed")
        t["disposition"] = "already_clean"
    assert "disposition" in refusal(tmp_path, sources, m)


def test_scored_disagreeing_with_censusability(tmp_path, sources):
    def m(d):
        t = next(t for t in d["benchmark_index"]["traces"]
                 if t["disposition"] == "not_censusable")
        t["scored"] = True
    assert "but disposition says" in refusal(tmp_path, sources, m)


def test_events_without_a_repaired_base_or_a_cut(tmp_path, sources):
    """The distinction the gate surfaced: events with no excision are
    admissible ONLY where displacement repair removed them."""
    def m(d):
        pinned = {s: e for s, e in d["release_pinned"]["segments"].items()}
        seg = next(s for s, e in pinned.items()
                   if e["disposition"] == "already_clean"
                   and e["base_kind"] == "displacement_repaired")
        pinned[seg]["base_kind"] = "original"
        for t in d["benchmark_index"]["traces"]:
            if t["segment"] == seg:
                break
    msg = refusal(tmp_path, sources, m)
    assert "must be one whose base was repaired" in msg


# --- a null or malformed field ------------------------------------------

def test_null_disposition(tmp_path, sources):
    def m(d):
        next(iter(d["release_pinned"]["segments"].values()))["disposition"] = None
    assert "disposition None is not one of" in refusal(tmp_path, sources, m)


def test_unknown_base_kind(tmp_path, sources):
    def m(d):
        next(iter(d["release_pinned"]["segments"].values()))["base_kind"] = "repaired"
    assert "base_kind 'repaired' is not one of" in refusal(tmp_path, sources, m)


def test_voxel_scale_as_a_string_is_not_a_scale(tmp_path, sources):
    def m(d):
        next(iter(d["release_pinned"]["segments"].values()))["voxel_um"] = "7.91"
    assert "positive finite number nor null" in refusal(tmp_path, sources, m)


def test_censusable_segment_with_no_contact_count(tmp_path, sources):
    def m(d):
        t = next(t for t in d["benchmark_index"]["traces"] if t.get("scored"))
        t["n_contacts"] = None
    assert "n_contacts is None" in refusal(tmp_path, sources, m)


# --- a flipped or malformed core verdict ---------------------------------

def test_flipped_core_verdict(tmp_path, sources):
    def m(d):
        d["acceptance"]["results"][0]["certificate"]["core_gate_pass"] = \
            not d["acceptance"]["results"][0]["certificate"]["core_gate_pass"]
    assert "core-gate verdict disagrees" in refusal(tmp_path, sources, m)


def test_core_verdict_as_a_string(tmp_path, sources):
    def m(d):
        d["acceptance"]["results"][0]["evaluator"]["core_gate_pass"] = "false"
    assert "not a boolean" in refusal(tmp_path, sources, m)


def test_clean_over_zero_triangles_is_not_clean(tmp_path, sources):
    """The mask defect's exact shape: transverse 0, triangles 0."""
    def m(d):
        d["acceptance"]["results"][0]["census_triangles"] = {"d0": 0, "d1": 0}
    msg = refusal(tmp_path, sources, m)
    assert "officially clean over positive triangles = censusable" in msg


def test_negative_transverse_count(tmp_path, sources):
    def m(d):
        d["acceptance"]["results"][0]["census_transverse"]["d0"] = -1
    assert "census transverse d0 is -1" in refusal(tmp_path, sources, m)


# --- the intrinsic half --------------------------------------------------

def test_intrinsic_row_for_a_segment_in_no_index(tmp_path, sources):
    def m(d):
        e = copy.deepcopy(d["intrinsic_pinned_d0"][0])
        e["segment"] = "a-segment-that-does-not-exist"
        d["intrinsic_pinned_d0"].append(e)
    assert "in no index" in refusal(tmp_path, sources, m)


def test_pinned_intrinsic_row_for_an_expansion_segment(tmp_path, sources):
    def m(d):
        e = copy.deepcopy(d["intrinsic_pinned_d0"][0])
        e["segment"] = next(iter(d["release_expansion"]["segments"]))
        d["intrinsic_pinned_d0"].append(e)
    assert "its inventory is 'expansion'" in refusal(tmp_path, sources, m)


def test_expansion_intrinsic_trace_not_ok(tmp_path, sources):
    def m(d):
        d["intrinsic_expansion"]["traces"][0]["ok"] = False
    assert "records ok=False" in refusal(tmp_path, sources, m)


def test_expansion_intrinsic_missing_a_diagonal(tmp_path, sources):
    def m(d):
        d["intrinsic_expansion"]["traces"][0]["diagonals"].pop("d1")
    assert "require exactly d0 and d1" in refusal(tmp_path, sources, m)


# --- diagonal coverage and measurement completeness ----------------------

def test_missing_d1_file_leaves_half_a_measurement(tmp_path, sources):
    """Without this, a missing d1 still leaves every pinned trace marked
    processed from d0 alone, with event totals quietly halved."""
    def m(d):
        d["intrinsic_pinned_d1"] = []
    msg = refusal(tmp_path, sources, m)
    assert "does not cover the pinned inventory exactly" in msg


@pytest.mark.parametrize("key", ["intrinsic_pinned_d0", "intrinsic_pinned_d1"])
def test_missing_one_pinned_row(tmp_path, sources, key):
    def m(d):
        d[key].pop()
    assert "does not cover the pinned inventory exactly" in \
        refusal(tmp_path, sources, key and m)


@pytest.mark.parametrize("key", ["intrinsic_pinned_d0", "intrinsic_pinned_d1"])
def test_duplicate_intrinsic_row(tmp_path, sources, key):
    def m(d):
        d[key].append(copy.deepcopy(d[key][0]))
    assert "appears twice in" in refusal(tmp_path, sources, m)


def test_pinned_event_not_measured(tmp_path, sources):
    def m(d):
        e = next(x for x in d["intrinsic_pinned_d0"] if x["events_total"] > 0)
        e["events_measured"] = e["events_total"] - 1
    assert "events" in refusal(tmp_path, sources, m)


def test_expansion_event_not_measured(tmp_path, sources):
    def m(d):
        t = d["intrinsic_expansion"]["traces"][0]
        t["diagonals"]["d0"]["events_measured"] = \
            t["diagonals"]["d0"]["events_total"] - 1
    assert "measured" in refusal(tmp_path, sources, m)


def test_per_diagonal_counts_are_preserved(tmp_path, sources):
    rec = run(tmp_path, sources)
    for r in rec["rows"]:
        if r["intrinsic_processed"]:
            assert r["intrinsic_events_d0"] is not None
            assert r["intrinsic_events_d1"] is not None
            assert (r["intrinsic_events_d0"] + r["intrinsic_events_d1"]
                    == r["intrinsic_events"])


def test_ambiguous_aggregates_name_their_operand(tmp_path, sources):
    c = run(tmp_path, sources)["claims"]
    assert "input base" in c["contacts_operand"]
    assert "published original" in c["defect_bearing_operand"]
    assert "not a cross-diagonal union" in c["intrinsic_event_semantics"]
    assert "NOT the denominator" in c["voxel_scale_note"]


# --- permissive types ----------------------------------------------------

def test_scored_as_a_string_is_not_true(tmp_path, sources):
    def m(d):
        next(t for t in d["benchmark_index"]["traces"])["scored"] = "false"
    assert "not a boolean" in refusal(tmp_path, sources, m)


def test_voxel_scale_of_zero_is_refused(tmp_path, sources):
    def m(d):
        next(iter(d["release_pinned"]["segments"].values()))["voxel_um"] = 0
    assert "positive finite number" in refusal(tmp_path, sources, m)


def test_fractional_contact_count_is_refused(tmp_path, sources):
    def m(d):
        t = next(t for t in d["benchmark_index"]["traces"] if t.get("scored"))
        t["n_contacts"] = 3.5
    assert "non-negative integer" in refusal(tmp_path, sources, m)


def test_retention_above_one_is_refused(tmp_path, sources):
    def m(d):
        t = next(t for t in d["benchmark_index"]["traces"] if t.get("scored"))
        t["retained_fraction_certificate"] = 1.4
    assert "finite fraction in [0, 1]" in refusal(tmp_path, sources, m)


def test_duplicate_expansion_intrinsic_row(tmp_path, sources):
    """A duplicated expansion row left processed and defect-bearing
    unchanged, so nothing noticed it."""
    def m(d):
        d["intrinsic_expansion"]["traces"].append(
            copy.deepcopy(d["intrinsic_expansion"]["traces"][0]))
    assert "appears twice in intrinsic_expansion" in refusal(tmp_path, sources, m)


def test_missing_expansion_intrinsic_row(tmp_path, sources):
    def m(d):
        d["intrinsic_expansion"]["traces"].pop()
    assert "does not cover the transformed expansion traces exactly" in \
        refusal(tmp_path, sources, m)
