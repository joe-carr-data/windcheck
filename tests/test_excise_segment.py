"""Driver-logic tests for bench/excise_segment.py -- NO real data, NO engine.

The real excision's authority is the C++ census, which these tests never
run. What they pin is the driver logic that stands between the engine and
the MILP, because that is where a real run can go quietly wrong:

1. The census-row -> coverage-set mapping (synthetic schema-v2 CSV rows):
   verdict filtering, round-24 eight-corner coverage, the ORIGINAL
   (v, u, tri) identities surviving the mapping, per-diagonal triangle
   corners, legacy-schema refusal, out-of-range accounting, and the two
   hard errors (retention divergence, adjacency-exclusion violation).
2. The shared-support (round-23 third class) REFUSAL path, both as the
   pure detector and end-to-end through `main()` with a mocked census:
   exit 0, a labelled refusal record, no output mesh, no clean claim.
3. The excision-fraction PRE-REGISTRATION gate: a fabricated over-budget
   solution must be refused before any mesh is emitted.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
import tifffile

sys.path.insert(0, "bench")
import excise_segment as es                                     # noqa: E402
from windcheck.check import PAIR_DTYPE                          # noqa: E402
from windcheck.intrinsic import oriented_events                 # noqa: E402

VOL = "20990101000000"
HEADER = es.SCHEMA_V2_HEADER


def write_csv(path: Path, rows) -> Path:
    """A synthetic AUTHORITATIVE-schema census CSV (engine sort order)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    body = [HEADER]
    for v1, u1, v2, u2, verdict, pen, ang, t1, t2 in rows:
        body.append(f"{v1},{u1},{v2},{u2},{verdict},{pen},{ang},{t1},{t2}")
    path.write_text("\n".join(body) + "\n")
    return path


# ------------------------------------------- 1. census row -> coverage sets
def test_parse_census_csv_keys_and_multiset(tmp_path):
    csv = write_csv(tmp_path / "c_d0.csv", [
        (2, 3, 7, 9, "transverse", 1.5, 30.0, 0, 1),
        (2, 3, 7, 9, "transverse", 0.5, 31.0, 1, 1),
        (4, 4, 8, 8, "grazing", 0.0, 2.0, 1, 0),
        (4, 4, 9, 9, "coplanar", 0.0, 0.0, 0, 0),
    ])
    out = es.parse_census_csv(csv, 0, 20, 20)
    assert out["header"] == HEADER
    assert out["n_lines"] == 4
    assert out["out_of_range"] == 0
    # keys keep the engine's ORIGINAL (v, u, tri) identities, per diagonal
    assert [r["key"] for r in out["rows"]] == [
        (0, 2, 3, 0, 7, 9, 1, "transverse"),
        (0, 2, 3, 1, 7, 9, 1, "transverse"),
        (0, 4, 4, 1, 8, 8, 0, "grazing"),
        (0, 4, 4, 0, 9, 9, 0, "coplanar")]
    assert sum(out["multiset"].values()) == 4
    assert out["multiset"][(0, 2, 3, 0, 7, 9, 1, "transverse")] == 1
    # the diagonal tag comes from the caller (one CSV per diagonal)
    assert all(k[0] == 1 for k in es.parse_census_csv(csv, 1, 20, 20)
               ["multiset"])


def test_parse_census_csv_out_of_range_counted_not_dropped(tmp_path):
    # nv = 6 -> quad rows only exist for v <= 4; v = 5 is out of range
    csv = write_csv(tmp_path / "c_d0.csv", [
        (0, 0, 3, 3, "transverse", 1.0, 20.0, 0, 0),
        (5, 0, 3, 3, "transverse", 1.0, 20.0, 0, 0),
        (0, 0, 3, 9, "transverse", 1.0, 20.0, 0, 0),
    ])
    out = es.parse_census_csv(csv, 0, 6, 10)
    assert len(out["rows"]) == 1
    assert out["out_of_range"] == 2


def test_parse_census_csv_refuses_legacy_and_unknown_schema(tmp_path):
    legacy = tmp_path / "legacy.csv"
    legacy.write_text("v1,u1,v2,u2,verdict,penetration,angle_deg\n"
                      "0,0,3,3,transverse,1,20\n")
    with pytest.raises(ValueError, match="legacy census CSV"):
        es.parse_census_csv(legacy, 0, 10, 10)
    odd = tmp_path / "odd.csv"
    odd.write_text("v1,u1,v2,u2,verdict,tri1,tri2\n")
    with pytest.raises(ValueError, match="unexpected census schema"):
        es.parse_census_csv(odd, 0, 10, 10)
    bad = write_csv(tmp_path / "bad.csv",
                    [(0, 0, 3, 3, "sideways", 1.0, 2.0, 0, 0)])
    with pytest.raises(ValueError, match="unknown verdicts"):
        es.parse_census_csv(bad, 0, 10, 10)


def all_retained(nv: int, nu: int) -> np.ndarray:
    return np.ones((nv - 1, nu - 1), dtype=bool)


def test_coverage_from_rows_eight_corner_quad_level(tmp_path):
    """Round-24 BLOCKER 1: coverage is corners(quad1) U corners(quad2) --
    eight vertices including each quad's non-participant fourth corner."""
    csv = write_csv(tmp_path / "c_d0.csv", [
        (2, 3, 7, 9, "transverse", 1.5, 30.0, 0, 1),
        (4, 4, 8, 8, "grazing", 0.0, 2.0, 1, 0),
        (4, 4, 9, 9, "coplanar", 0.0, 0.0, 0, 0),
    ])
    rows = es.parse_census_csv(csv, 0, 20, 20)["rows"]
    cons = es.coverage_from_rows(rows, all_retained(20, 20))
    # only TRANSVERSE rows become constraints (grazing/coplanar are
    # contacts, not crossings)
    assert len(cons) == 1
    c = cons[0]
    assert c["key"] == (0, 2, 3, 0, 7, 9, 1, "transverse")
    assert (c["q1"], c["t1"], c["q2"], c["t2"]) == ((2, 3), 0, (7, 9), 1)
    assert c["coverage"] == sorted(
        {(2, 3), (3, 3), (2, 4), (3, 4), (7, 9), (8, 9), (7, 10), (8, 10)})
    assert len(c["coverage"]) == 8
    # d0 triangle 0 of a quad is (00, 01, 11); triangle 1 is (00, 11, 10)
    assert c["corners1"] == ((2, 3), (2, 4), (3, 4))
    assert c["corners2"] == ((7, 9), (8, 10), (8, 9))
    # the fourth corners are in coverage but are NOT triangle participants
    assert (3, 3) in c["coverage"] and (3, 3) not in c["participants"]
    assert (7, 10) in c["coverage"] and (7, 10) not in c["participants"]
    assert set(c["participants"]) < set(c["coverage"])


def test_coverage_from_rows_diagonal_specific_triangles(tmp_path):
    """d1 tessellation: triangle 0 is (00, 01, 10), triangle 1 is
    (01, 11, 10) -- the mapping must use the CSV's own diagonal."""
    csv = write_csv(tmp_path / "c_d1.csv",
                    [(2, 3, 7, 9, "transverse", 1.5, 30.0, 1, 0)])
    rows = es.parse_census_csv(csv, 1, 20, 20)["rows"]
    c = es.coverage_from_rows(rows, all_retained(20, 20))[0]
    assert c["diag"] == 1
    assert c["corners1"] == ((2, 4), (3, 4), (3, 3))
    assert c["corners2"] == ((7, 9), (7, 10), (8, 9))
    # coverage is diagonal-independent: quad-level validity
    assert c["coverage"] == sorted(
        {(2, 3), (3, 3), (2, 4), (3, 4), (7, 9), (8, 9), (7, 10), (8, 10)})


def test_coverage_from_rows_hard_errors(tmp_path):
    csv = write_csv(tmp_path / "c_d0.csv",
                    [(2, 3, 7, 9, "transverse", 1.5, 30.0, 0, 1)])
    rows = es.parse_census_csv(csv, 0, 20, 20)["rows"]
    Q = all_retained(20, 20)
    Q[7, 9] = False                       # engine censused it, Python did not
    with pytest.raises(ValueError, match="retention semantics diverge"):
        es.coverage_from_rows(rows, Q)
    # adjacency exclusion: quads sharing a corner must never be censused
    adj = write_csv(tmp_path / "adj_d0.csv",
                    [(2, 3, 3, 4, "transverse", 1.5, 30.0, 0, 1)])
    rows2 = es.parse_census_csv(adj, 0, 20, 20)["rows"]
    with pytest.raises(ValueError, match="adjacency exclusion violated"):
        es.coverage_from_rows(rows2, all_retained(20, 20))


def test_transverse_pair_records_dedups_quad_pairs(tmp_path):
    csv = write_csv(tmp_path / "c_d0.csv", [
        (2, 3, 7, 9, "transverse", 1.5, 30.0, 0, 1),
        (2, 3, 7, 9, "transverse", 0.5, 31.0, 1, 1),
        (2, 3, 7, 10, "transverse", 0.5, 31.0, 1, 1),
        (4, 4, 8, 8, "grazing", 0.0, 2.0, 1, 0),
    ])
    rows = es.parse_census_csv(csv, 0, 20, 20)["rows"]
    rec = es.transverse_pair_records(rows)
    assert len(rec) == 2                  # quad-level, grazing excluded
    assert sorted((int(r["v1"]), int(r["u1"]), int(r["v2"]), int(r["u2"]))
                  for r in rec) == [(2, 3, 7, 9), (2, 3, 7, 10)]


def test_witness_scopes_triangle_participant_vs_quad_retention():
    key = (0, 2, 3, 0, 7, 9, 1, "transverse")
    # (2,4) is a corner of quad1's d0 triangle 0 -> triangle_participant
    w = es.witness_for_key(key, {(2, 4)})
    assert w["witness_scope"] == "triangle_participant"
    assert w["destroyed_triangle"] == [0, 2, 3, 0]
    # (3,3) is quad1's fourth corner, absent from that triangle
    w = es.witness_for_key(key, {(3, 3)})
    assert w["witness_scope"] == "quad_retention"
    assert w["destroyed_triangle"] == [0, 2, 3, 0]
    assert es.witness_for_key(key, {(0, 0)}) is None


# --------------------------------------------- 2. shared-support detection
def rec_of(pairs) -> np.ndarray:
    return np.array([(a[0], a[1], b[0], b[1], 1.0, 30.0) for a, b in pairs],
                    dtype=PAIR_DTYPE)


def touching_pairs():
    """Three chained pairs whose parity-normalised regions END UP
    grid-adjacent: A = {(0,0),(0,1),(0,2)}, B = {(0,3),(0,4),(0,5)}, and
    (0,2) touches (0,3). Every emitted pair is itself at Chebyshev
    distance 3, so the engine's exclude=1 would have emitted all of them."""
    return [((0, 0), (0, 3)), ((0, 1), (0, 4)), ((0, 2), (0, 5))]


def separated_pairs():
    return [((0, 0), (0, 5)), ((0, 1), (0, 6))]


def test_shared_support_events_detects_third_class():
    evs = oriented_events(rec_of(touching_pairs()))
    assert len(evs) == 1 and evs[0]["self_touching"]
    hits = es.shared_support_events({0: evs, 1: []})
    assert len(hits) == 1
    assert hits[0]["diagonal"] == 0
    assert hits[0]["region_a"] and hits[0]["region_b"]

    clean = oriented_events(rec_of(separated_pairs()))
    assert len(clean) == 1 and not clean[0]["self_touching"]
    assert es.shared_support_events({0: clean, 1: clean}) == []


def test_shared_support_label_is_the_spec_wording():
    assert "shared-support topology" in es.SHARED_SUPPORT_LABEL
    assert "remesh/junction-excision required" in es.SHARED_SUPPORT_LABEL


# ------------------------------------------------- 3. the pre-reg budget gate
def test_budget_verdict_gate():
    ok = es.budget_verdict(1.0, 1000.0, 0.01)
    assert ok["within_budget"] and ok["excised_fraction"] == 0.001
    bad = es.budget_verdict(50.0, 1000.0, 0.01)
    assert not bad["within_budget"] and bad["excised_fraction"] == 0.05
    assert "REFUSED" in bad["label"]
    # exactly at the ceiling is inside it (the pre-registered comparison
    # is <=, recorded in the certificate as such)
    assert es.budget_verdict(10.0, 1000.0, 0.01)["within_budget"]
    assert es.budget_verdict(1.0, 0.0, 0.01)["within_budget"] is False


def test_area_block_partition_identity():
    areas = np.full((4, 4), 2.0)
    Q = np.ones((4, 4), dtype=bool)
    b = es.area_block(areas, Q, [(0, 0), (1, 1)], [(2, 2)])
    assert b["A_input"] == 32.0
    assert b["A_excised"] == 4.0
    assert b["A_unresolved"] == 2.0
    assert b["A_clean"] == 26.0
    assert b["identity_residual"] == 0.0
    assert b["excised_fraction"] == 4.0 / 32.0


def test_cut_boundary_only_edges_shared_with_retained_quads():
    # one removed quad (1,1) surrounded by retained quads on two sides
    edges = es.cut_boundary({(1, 1)}, {(0, 1), (1, 2)})
    assert edges == sorted([tuple(sorted(((1, 1), (1, 2)))),
                            tuple(sorted(((1, 2), (2, 2))))])
    assert es.cut_boundary({(1, 1)}, set()) == []


# ------------------------------------------------ end-to-end, census mocked
def make_segment(tmp_path: Path, nv: int = 3, nu: int = 10) -> tuple:
    """A synthetic corpus with one tiny, real, readable tifxyz mesh."""
    root = tmp_path / "corpus"
    seg = root / "segTEST"
    mesh = seg / "mesh" / f"segTEST-on-{VOL}-9.0um.tifxyz"
    mesh.mkdir(parents=True)
    P = np.zeros((nv, nu, 3), np.float32)
    for v in range(nv):
        for u in range(nu):
            P[v, u] = (float(v), float(u), 0.0)
    for i, ax in enumerate(("x", "y", "z")):
        tifffile.imwrite(mesh / f"{ax}.tif", P[..., i])
    (mesh / "meta.json").write_text(json.dumps({"scale": [1.0, 1.0],
                                                "bbox": [[0, 0, 0],
                                                         [1, 1, 1]]}))
    return [("TestCorpus", str(root), VOL, "out/unused")], mesh


def fake_census(pairs, transverse_tris=(0, 0)):
    """A stand-in for `crossing_census.census_one`: writes the CSVs the
    driver parses and returns the engine's counts dict shape."""
    def census_one(path, name, exclude, cell, threads, maxedge, work):
        work = Path(work)
        work.mkdir(parents=True, exist_ok=True)
        rows = [(a[0], a[1], b[0], b[1], "transverse", 1.0, 30.0,
                 transverse_tris[0], transverse_tris[1]) for a, b in pairs]
        counts = {"segment": name, "grid": None, "valid_cells": 0}
        for d in (0, 1):
            write_csv(work / f"{name[:40]}_d{d}.csv", rows)
            counts[f"d{d}"] = {"triangles": 0, "quads_dropped": 0,
                               "pairs_tested": 0, "transverse": len(rows),
                               "coplanar": 0, "grazing": 0}
        return counts
    return census_one


def run_main(monkeypatch, tmp_path, corpora, census_one, extra=(), argv=()):
    monkeypatch.setattr(es, "CORPORA", corpora)
    monkeypatch.setattr(es, "OUT", tmp_path / "out")
    monkeypatch.setattr(es, "census_one", census_one)
    for name, val in extra:
        monkeypatch.setattr(es, name, val)
    monkeypatch.setattr(sys, "argv",
                        ["excise_segment.py", "--segment", "segTEST",
                         *argv])
    return es.main()


def test_main_refuses_shared_support_cleanly(monkeypatch, tmp_path):
    """Round-23 third class, still reachable behind --refuse-shared-support:
    labelled non-output, exit 0, no mesh, and the refusal happens BEFORE any
    solve (the solver is booby-trapped)."""
    corpora, _mesh = make_segment(tmp_path)

    def boom(*a, **k):
        raise AssertionError("solver must never run on a refused input")

    code = run_main(monkeypatch, tmp_path, corpora,
                    fake_census(touching_pairs()),
                    extra=[("solve_global", boom)],
                    argv=["--refuse-shared-support"])
    assert code == 0
    rec = json.loads((tmp_path / "out"
                      / "segTEST_excision_refusal.json").read_text())
    assert rec["status"] == "refused_shared_support"
    assert rec["clean_claim"] is False
    assert rec["label"] == es.SHARED_SUPPORT_LABEL
    assert rec["output_mesh"] is None
    assert rec["shared_support_events"]
    assert not (tmp_path / "out" / "segTEST.tifxyz").exists()
    assert not (tmp_path / "out"
                / "segTEST_excision_certificate.json").exists()


def test_main_refuses_over_budget_solution(monkeypatch, tmp_path):
    """The pre-registration gate: a fabricated over-budget cut is refused
    BEFORE any mesh is emitted, and the refusal record carries the gate."""
    corpora, _mesh = make_segment(tmp_path)

    def fake_solve(constraints, P64, Q_in, protected, time_limit):
        # invalidate a whole column of vertices: removes ~2/16 of the quads,
        # far above the 1% pre-registered ceiling
        chosen = [(v, 4) for v in range(P64.shape[0])]
        return {"status": "optimal", "chosen": chosen, "lexicographic": True,
                "scipy_version": "fabricated",
                "reduction": {"n_raw": 1, "n_after_dedup": 1,
                              "n_after_dominance": 1, "n_components": 1,
                              "component_sizes": [1], "largest_component": 1},
                "records": [{"stage": 1, "purpose": "min_excised_area",
                             "status": "optimal", "raw_status": 0,
                             "message": "fabricated", "objective": 99.0,
                             "dual_bound": 99.0, "mip_gap": 0.0,
                             "solve_time_s": 0.0, "n_vertex_vars": 0,
                             "n_quad_vars": 0, "n_constraints": 0}]}

    code = run_main(
        monkeypatch, tmp_path, corpora, fake_census(separated_pairs()),
        extra=[("solve_global", fake_solve),
               ("event_bounds", lambda g, X, rec, ev: (
                   3.0, 2.9, {"allowance_vx": 1e-3}, 4))])
    assert code == 1
    rec = json.loads((tmp_path / "out"
                      / "segTEST_excision_refusal.json").read_text())
    assert rec["status"] == "refused_over_budget"
    assert rec["clean_claim"] is False
    gate = rec["budget_gate"]
    assert gate["within_budget"] is False
    assert gate["max_excision_fraction"] == es.MAX_EXCISION_FRACTION
    assert gate["excised_fraction"] > es.MAX_EXCISION_FRACTION
    assert gate["basis"] == es.MAX_EXCISION_FRACTION_BASIS
    # the pre-registration was recorded before the solve, and the mesh was
    # never emitted
    pre = rec["pre_registration"]
    assert pre["recorded_before_any_solve"] is True
    assert pre["max_excision_fraction"] == es.MAX_EXCISION_FRACTION
    assert pre["target_events"][0]["min_exit_bounds"]["L_safe_vx"] == 2.9
    assert pre["target_events"][0]["rigid_verdict"] == "certified_infeasible"
    assert (tmp_path / "out" / "segTEST_prereg.json").exists()
    assert not (tmp_path / "out" / "segTEST.tifxyz").exists()
    assert not (tmp_path / "out"
                / "segTEST_excision_certificate.json").exists()


# ---------------- 4. round-25/26: hybrid emission, naive reader, junction
def fake_census_then_clean(pairs, transverse_tris=(0, 0)):
    """Census stand-in that reports `pairs` on the BEFORE census and NOTHING
    on every later one -- i.e. the emitted cut actually worked."""
    def census_one(path, name, exclude, cell, threads, maxedge, work):
        work = Path(work)
        work.mkdir(parents=True, exist_ok=True)
        first = name.endswith("_before")
        rows = ([(a[0], a[1], b[0], b[1], "transverse", 1.0, 30.0,
                  transverse_tris[0], transverse_tris[1])
                 for a, b in pairs] if first else [])
        counts = {"segment": name, "grid": None, "valid_cells": 0}
        for d in (0, 1):
            write_csv(work / f"{name[:40]}_d{d}.csv", rows)
            counts[f"d{d}"] = {"triangles": 0, "quads_dropped": 0,
                               "pairs_tested": 0, "transverse": len(rows),
                               "coplanar": 0, "grazing": 0}
        return counts
    return census_one


def naive_read(mesh: Path):
    """A consumer implementing ONLY the x=y=z=-1 convention: it never opens
    mask.tif. This is the round-25 A1 interoperability test."""
    P = np.stack([np.asarray(tifffile.imread(mesh / f"{ax}.tif"),
                             dtype=np.float32) for ax in ("x", "y", "z")],
                 axis=-1)
    return P, ~np.all(P == -1.0, axis=-1)


def test_main_hybrid_invalidation_and_naive_consumer(monkeypatch, tmp_path):
    """END TO END: the emitted mesh must be invalidated BOTH ways, and a
    reader that ignores mask.tif must see exactly the same retained surface
    (round-25 A1). The naive-consumer census is part of the certificate."""
    corpora, mesh_in = make_segment(tmp_path)
    code = run_main(monkeypatch, tmp_path, corpora,
                    fake_census_then_clean(separated_pairs()),
                    argv=["--max-excision-fraction", "0.5", "--skip-bounds"])
    assert code == 0
    out = tmp_path / "out" / "segTEST.tifxyz"
    cert = json.loads((tmp_path / "out"
                       / "segTEST_excision_certificate.json").read_text())
    assert cert["status"] == "clean" and cert["clean_claim"] is True

    # BOTH carriers are present on disk
    assert (out / "mask.tif").exists()
    mask = np.asarray(tifffile.imread(out / "mask.tif")).astype(bool)
    P_out, V_naive = naive_read(out)
    P_in, V_in = naive_read(mesh_in)
    assert not V_in.all() or V_in.all()          # input has no stamps to hide
    excised = V_in & ~mask
    assert excised.any(), "the test must actually cut something"

    # the two conventions agree cell for cell -- that is the guarantee
    assert np.array_equal(V_naive, mask)
    # excised cells carry the -1 marker; retained coordinates are untouched
    assert (P_out[excised] == -1.0).all()
    assert P_out[mask].tobytes() == P_in[mask].tobytes()
    changed = np.any(P_out != P_in, axis=-1)
    assert np.array_equal(changed, excised)

    assert cert["retained_coordinate_bit_identity"] is True
    assert cert["coordinates_changed_only_at_excised_cells"] is True
    assert cert["reload_checks"]["excised_cells_stamped_missing"] is True
    nc = cert["invalidation"]["naive_consumer_check"]
    assert nc["has_mask_sidecar"] is False        # the naive copy had none
    assert nc["valid_equals_masked_reading"] is True
    assert nc["naive_reader_sees_no_transverse_crossing"] is True
    assert cert["acceptance"]["naive_consumer_transverse_clean"] is True
    assert not (tmp_path / "out" / "segTEST_naive.tifxyz").exists()


def test_main_segment_wide_one_milp_and_reduction_recorded(monkeypatch,
                                                           tmp_path):
    """Round-26: ONE MILP over every transverse row of both diagonals, one
    emission, one recensus -- and the reduction statistics are recorded."""
    corpora, _mesh = make_segment(tmp_path)
    code = run_main(monkeypatch, tmp_path, corpora,
                    fake_census_then_clean(separated_pairs()),
                    argv=["--max-excision-fraction", "0.5", "--skip-bounds",
                          "--skip-naive-check"])
    assert code == 0
    cert = json.loads((tmp_path / "out"
                       / "segTEST_excision_certificate.json").read_text())
    solver = cert["solver"]
    assert solver["iterations"] == 1              # emit once, recensus once
    assert "SEGMENT-WIDE" in solver["constraint_scope"]
    red = solver["constraint_reduction"]
    # 2 pairs x 2 diagonals = 4 rows in the census
    assert red["n_raw"] == 4
    assert red["n_after_dedup"] <= red["n_raw"]
    assert red["n_after_dominance"] <= red["n_after_dedup"]
    assert sum(red["component_sizes"]) == red["n_after_dominance"]
    assert solver["final_stage1_objective"] is not None
    assert solver["optimality_gap"] is not None
    for s in solver["solves"]:
        assert "objective" in s and "dual_bound" in s and "mip_gap" in s
    assert "attributable" not in json.dumps(cert["acceptance"])
    assert cert["acceptance"]["transverse_clean_both_diagonals"] is True


def test_main_shared_support_is_junction_excision_not_a_refusal(monkeypatch,
                                                                tmp_path):
    """Round-26 Q1: shared support is cut and LABELLED, not refused."""
    corpora, _mesh = make_segment(tmp_path)
    code = run_main(monkeypatch, tmp_path, corpora,
                    fake_census_then_clean(touching_pairs()),
                    argv=["--max-excision-fraction", "0.5", "--skip-bounds",
                          "--skip-naive-check"])
    assert code == 0
    assert not (tmp_path / "out"
                / "segTEST_excision_refusal.json").exists()
    cert = json.loads((tmp_path / "out"
                       / "segTEST_excision_certificate.json").read_text())
    assert cert["operation_label"].startswith("junction_excision")
    assert "NEVER a branch separation" in cert["operation_label"]
    assert cert["acceptance"]["shared_support_events"]
    assert cert["status"] == "clean"


def test_solve_totals_sums_components():
    recs = [{"stage": 1, "objective": 2.0, "dual_bound": 2.0,
             "mip_gap": 0.0, "solve_time_s": 1.0, "status": "optimal"},
            {"stage": 1, "objective": 3.0, "dual_bound": 1.0,
             "mip_gap": 0.66, "solve_time_s": 2.0, "status": "best_found"},
            {"stage": 2, "objective": 9.0, "dual_bound": 9.0,
             "mip_gap": 0.0, "solve_time_s": 0.5, "status": "optimal"}]
    t = es.solve_totals(recs, 1)
    assert t["objective"] == 5.0 and t["dual_bound"] == 3.0
    assert abs(t["mip_gap"] - 0.4) < 1e-12
    assert t["seconds"] == 3.0 and t["n_solves"] == 2
    assert t["statuses"] == ["best_found", "optimal"]
    assert es.solve_totals([], 1)["objective"] is None
