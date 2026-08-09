"""Regressions for the quickstart gate.

The gate exists because the quickstart used to treat any non-zero exit on
the ORIGINAL trace as success. Every test below is a report that would
have printed as a passing demonstration under that rule, or a binding
that a merely-clean surface would satisfy without being the packaged
reference. Feeding reports directly means these run without a validator,
a container, or a 300 MB trace.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bench"))
from quickstart_gate import check  # noqa: E402


def _write(p: Path, obj) -> Path:
    p.write_text(json.dumps(obj))
    return p


def dirty_report(inp: Path, cand: Path, **over) -> dict:
    d = {
        "input": str(inp), "candidate": str(cand),
        "comparable": True,
        "input_identity": {"verified": True},
        "census": {"ran": True, "authoritative": True, "parameters_ok": True,
                   "transverse": {"d0": 11, "d1": 14},
                   "triangles": {"d0": 500, "d1": 500}},
        "clean": False, "scored": False,
    }
    d.update(over)
    return d


def clean_report(inp: Path, cand: Path, **over) -> dict:
    d = {
        "input": str(inp), "candidate": str(cand),
        "comparable": True,
        "input_identity": {"verified": True},
        "census": {"ran": True, "authoritative": True, "parameters_ok": True,
                   "transverse": {"d0": 0, "d1": 0},
                   "triangles": {"d0": 500, "d1": 500}},
        "clean": True, "scored": True,
        "retained_fraction": 0.9999685874068042,
        "coordinate_fidelity": True,
        "added_cells": 0,
        "fragmentation_comparable": True,
    }
    d.update(over)
    return d


@pytest.fixture
def dirs(tmp_path):
    inp = tmp_path / "input.tifxyz"
    cand = tmp_path / "reference.tifxyz"
    inp.mkdir()
    cand.mkdir()
    return inp, cand


# --- the shapes that used to pass ---------------------------------------

def test_dirty_report_accepted(tmp_path, dirs):
    inp, cand = dirs
    r = _write(tmp_path / "r.json", dirty_report(inp, inp))
    assert check(r, "dirty", 3, inp, inp) == []


def test_missing_report_is_not_a_pass(tmp_path, dirs):
    inp, _ = dirs
    bad = check(tmp_path / "absent.json", "dirty", 3, inp, inp)
    assert bad and "wrote no report" in bad[0]


def test_dead_validator_is_not_a_dirty_demonstration(tmp_path, dirs):
    """The census never ran; exit was non-zero for the wrong reason."""
    inp, _ = dirs
    d = dirty_report(inp, inp)
    d["census"] = {"ran": False, "detail": "Cannot open: .../meta.json"}
    r = _write(tmp_path / "r.json", d)
    bad = check(r, "dirty", 2, inp, inp)
    assert bad == ["the census did not run: Cannot open: .../meta.json"]


def test_vacuous_census_over_zero_triangles_is_rejected(tmp_path, dirs):
    inp, _ = dirs
    d = dirty_report(inp, inp)
    d["census"]["triangles"] = {"d0": 0, "d1": 0}
    r = _write(tmp_path / "r.json", d)
    bad = check(r, "dirty", 3, inp, inp)
    assert any("vacuous" in b for b in bad)


def test_a_refusal_is_not_a_rejection(tmp_path, dirs):
    inp, _ = dirs
    d = dirty_report(inp, inp, refused="census parameters differ")
    r = _write(tmp_path / "r.json", d)
    bad = check(r, "dirty", 2, inp, inp)
    assert bad == ["the evaluator refused: census parameters differ"]


def test_precomputed_census_is_not_authoritative(tmp_path, dirs):
    inp, _ = dirs
    d = dirty_report(inp, inp)
    d["census"]["authoritative"] = False
    r = _write(tmp_path / "r.json", d)
    assert any("authoritative" in b for b in check(r, "dirty", 3, inp, inp))


def test_stale_report_from_another_operand_is_rejected(tmp_path, dirs):
    """A report left by a previous invocation must not be read as this
    run's result."""
    inp, cand = dirs
    r = _write(tmp_path / "r.json",
               dirty_report(tmp_path / "some-other-trace", inp))
    bad = check(r, "dirty", 3, inp, inp)
    assert any("the report is for input" in b for b in bad)


def test_unverified_input_is_rejected(tmp_path, dirs):
    inp, _ = dirs
    d = dirty_report(inp, inp)
    d["input_identity"] = {"verified": False}
    r = _write(tmp_path / "r.json", d)
    assert any("NOT verified by hash" in b for b in check(r, "dirty", 3, inp, inp))


def test_clean_original_does_not_demonstrate_rejection(tmp_path, dirs):
    inp, _ = dirs
    d = dirty_report(inp, inp)
    d["census"]["transverse"] = {"d0": 0, "d1": 0}
    r = _write(tmp_path / "r.json", d)
    assert any("does not demonstrate rejection" in b
               for b in check(r, "dirty", 3, inp, inp))


def test_wrong_exit_status_is_rejected(tmp_path, dirs):
    inp, _ = dirs
    r = _write(tmp_path / "r.json", dirty_report(inp, inp))
    assert any("exit status 2" in b for b in check(r, "dirty", 2, inp, inp))


def test_malformed_json_is_rejected(tmp_path, dirs):
    inp, _ = dirs
    p = tmp_path / "r.json"
    p.write_text("{not json")
    assert check(p, "dirty", 3, inp, inp)


def test_true_is_not_a_triangle_count(tmp_path, dirs):
    """`True == 1` in Python, so a boolean must not read as a count."""
    inp, _ = dirs
    d = dirty_report(inp, inp)
    d["census"]["triangles"] = {"d0": True, "d1": True}
    r = _write(tmp_path / "r.json", d)
    assert any("vacuous" in b for b in check(r, "dirty", 3, inp, inp))


# --- the clean side, and its binding to the packaged reference ----------

def _make_planes(d: Path, payload: bytes) -> dict:
    hashes = {}
    for a in "xyz":
        (d / f"{a}.tif").write_bytes(payload + a.encode())
        hashes[a] = hashlib.sha256(payload + a.encode()).hexdigest()
    return hashes


CANONICAL_RF = 0.9999685874068042


@pytest.fixture
def bound(tmp_path, dirs):
    """A candidate that IS the packaged reference, with the two records
    the clean gate requires."""
    inp, cand = dirs
    real = _make_planes(cand, b"the packaged reference")
    prov = _write(tmp_path / "provenance.json",
                  {"reference_derivative": {"hashes": real}})
    fields = _write(tmp_path / "fields.json",
                    {"retention": {"area": {"canonical": {
                        "retained_fraction": CANONICAL_RF}}}})
    return inp, cand, prov, fields


def test_clean_report_accepted(tmp_path, bound):
    inp, cand, prov, fields = bound
    r = _write(tmp_path / "r.json", clean_report(inp, cand))
    assert check(r, "clean", 0, inp, cand, provenance=prov,
                 fields=fields) == []


@pytest.mark.parametrize("missing", ["provenance", "fields"])
def test_clean_mode_requires_both_bindings(tmp_path, bound, missing):
    """The module's contract is that a clean pass means THE packaged
    reference. Without these it would be the weaker claim."""
    inp, cand, prov, fields = bound
    kw = {"provenance": prov, "fields": fields}
    kw[missing] = None
    bad = check(_write(tmp_path / "r.json", clean_report(inp, cand)),
                "clean", 0, inp, cand, **kw)
    assert any(f"no {'provenance' if missing == 'provenance' else 'fields.json'}"
               in b for b in bad)


@pytest.mark.parametrize("field,value,needle", [
    ("coordinate_fidelity", False, "MOVED geometry"),
    ("added_cells", 7, "added 7 cells"),
    ("fragmentation_comparable", False, "not comparable"),
    ("retained_fraction", None, "not a finite number"),
    ("retained_fraction", float("nan"), "not a finite number"),
    ("retained_fraction", 1.5, "outside [0, 1]"),
])
def test_clean_but_not_an_excision_is_rejected(tmp_path, bound, field, value,
                                               needle):
    inp, cand, prov, fields = bound
    r = _write(tmp_path / "r.json", clean_report(inp, cand, **{field: value}))
    assert any(needle in b for b in check(r, "clean", 0, inp, cand,
                                          provenance=prov, fields=fields))


def test_some_other_clean_surface_is_not_the_reference(tmp_path, bound):
    """The load-bearing binding: a corrupted masks.npz that removed far
    more surface would still be clean, still be an excision, and still
    score. Only the coordinate hashes rule it out."""
    inp, cand, prov, fields = bound
    _make_planes(cand, b"a different clean derivative")
    bad = check(_write(tmp_path / "r.json", clean_report(inp, cand)),
                "clean", 0, inp, cand, provenance=prov, fields=fields)
    assert sum("it is not THE reference" in b for b in bad) == 3


def test_retention_must_match_the_shipped_certificate_figure(tmp_path, bound):
    inp, cand, prov, fields = bound
    off = clean_report(inp, cand, retained_fraction=0.97)
    bad = check(_write(tmp_path / "r.json", off), "clean", 0, inp, cand,
                provenance=prov, fields=fields)
    assert any("does not match the certificate figure" in b for b in bad)


@pytest.mark.parametrize("hashes,needle", [
    ({}, "no valid x hash"),
    ("PARTIAL", "no valid y hash"),
    ({"x": "not-a-hash", "y": "b" * 64, "z": "c" * 64}, "no valid x hash"),
])
def test_incomplete_reference_hashes_are_refused(tmp_path, bound, hashes,
                                                 needle):
    """A provenance record carrying only `x` must not satisfy a binding
    this module describes as x/y/z."""
    inp, cand, _, fields = bound
    real = _make_planes(cand, b"the packaged reference")
    if hashes == "PARTIAL":
        hashes = {"x": real["x"]}
    prov = _write(tmp_path / "prov2.json",
                  {"reference_derivative": {"hashes": hashes}})
    bad = check(_write(tmp_path / "r.json", clean_report(inp, cand)),
                "clean", 0, inp, cand, provenance=prov, fields=fields)
    assert any(needle in b for b in bad)


def test_nonfinite_certificate_figure_is_refused(tmp_path, bound):
    inp, cand, prov, _ = bound
    fields = _write(tmp_path / "f2.json",
                    {"retention": {"area": {"canonical": {
                        "retained_fraction": None}}}})
    bad = check(_write(tmp_path / "r.json", clean_report(inp, cand)),
                "clean", 0, inp, cand, provenance=prov, fields=fields)
    assert any("no finite canonical retained fraction" in b for b in bad)
