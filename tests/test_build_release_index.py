"""Tests for the release index generator (bench/build_release_index.py).

Everything here is SYNTHETIC and lives under tmp_path: a fabricated base
manifest, three fabricated certificates covering the three dispositions, a
fabricated driver summary and a fabricated recensus record.  No real corpus
is touched, so these pass on a checkout with no out/ directory at all and
they run in well under a second.

What is pinned:

  * the three terminal dispositions map to the three public names, and the
    non-censusable input survives into the index with its reason attached
    rather than being dropped;
  * both area denominators are carried through per segment, and the
    area-weighted aggregates use the denominator they claim to use;
  * a field the certificate does not record comes out null with a note,
    never filled in;
  * the archive block measures real bytes on disk.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
for _p in (REPO_ROOT / "src", REPO_ROOT / "bench"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import build_release_index as bri  # noqa: E402


# ------------------------------------------------------------------ helpers

def write_mesh(root: Path, segment: str, payload: bytes = b"plane") -> str:
    mesh = root / "out" / "excised" / "corpus" / f"{segment}_excised.tifxyz"
    mesh.mkdir(parents=True, exist_ok=True)
    for ax in ("x", "y", "z"):
        (mesh / f"{ax}.tif").write_bytes(payload + ax.encode())
    return f"out/excised/corpus/{segment}_excised.tifxyz"


def base_cert(segment: str, scroll: str, disposition: str) -> dict:
    return {
        "record_kind": "excision certificate",
        "segment": segment,
        "corpus": scroll,
        "volume": "vol",
        "voxel_um": 7.91,
        "terminal_disposition": disposition,
        "base_kind": "original",
        "base_mesh": f"data/{segment}.tifxyz",
        "base_hashes": {"x": "aa" * 32, "y": "bb" * 32, "z": "cc" * 32,
                        "mask": None},
        "original_mesh": f"data/{segment}.tifxyz",
        "original_hashes": {"x": "aa" * 32, "y": "bb" * 32, "z": "cc" * 32,
                            "mask": None},
        "geometry_key": "gk-" + segment,
        "is_canonical": True,
        "duplicate_of": None,
        "base_manifest_sha256": "dd" * 32,
        "policy_version": "test-policy",
        "policy_hash": "0" * 16,
        "grid_shape": [10, 10],
        "n_valid": 100,
    }


def build_fixture(tmp_path: Path) -> Path:
    """A three-segment corpus: transformed, already-clean, not censusable."""
    corpus = tmp_path / "out" / "excised" / "corpus"
    corpus.mkdir(parents=True)
    (corpus / "logs").mkdir()
    rows = []

    # --- transformed: 90% of a 1000-unit base, base IS the original -------
    seg = "seg_transformed"
    cert = base_cert(seg, "Scroll 1", "transformed")
    cert.update({
        "status": "transverse_clean",
        "output_mesh": write_mesh(tmp_path, seg),
        "output_mesh_hashes": {"x": "11" * 32, "y": "22" * 32,
                               "z": "33" * 32},
        "wall_seconds": 10.0,
        "input_area_canonical": 1000.0,
        "input_transverse_total": 4,
        "output_transverse_total": 0,
        "geometry_status": "transverse_clean_certified",
        "area": {"canonical": {"A_input": 1000.0, "A_retained": 900.0,
                               "A_excised": 100.0,
                               "retained_fraction": 0.9}},
        "headline_area": {"measurable": True,
                          "A_original_canonical": 1000.0,
                          "A_removed_priced_on_original": 100.0,
                          "headline_retained_fraction": 0.9},
        "excision": {"n_removed_quads": 7},
        "component_recovery": {
            "components_before": 2, "components_after": 3,
            "core_gate": {"n_input_components": 2, "n_core_components": 2,
                          "min_R_main_core": 0.5,
                          "n_core_components_below_gate": 1,
                          "core_gate_pass": False},
            "unthresholded": {"min_R_main_all_components": 0.5,
                              "area_weighted_R_main": 0.8,
                              "n_input_components_fully_destroyed": 0},
        },
        "core_gate_pass": False,
        "operational_retained_fraction": 0.9,
        "headline_retained_fraction": 0.9,
        "claimed_clean": True,
        "selection_status": "area_optimal",
    })
    rows.append(cert)

    # --- already clean, on a repaired base: no priced-area block ----------
    seg = "seg_clean"
    cert = base_cert(seg, "Scroll 1", "already_clean")
    cert.update({
        "status": "already_transverse_clean",
        "base_kind": "displacement_repaired",
        # the base that was censused is the REPAIRED mesh, not the published
        # original, so no original-coordinate area is recorded anywhere
        "base_mesh": f"out/repaired/{seg}_repaired.tifxyz",
        "input_mesh": f"out/repaired/{seg}_repaired.tifxyz",
        "base_hashes": {"x": "ee" * 32, "y": "ff" * 32, "z": "0e" * 32,
                        "mask": None},
        "output_mesh": None,
        "output_mesh_hashes": None,
        "wall_seconds": 0.5,
        "input_area_canonical": 500.0,
        "input_transverse_total": 0,
        "geometry_status": "transverse_clean_certified",
        "note": "no excision required; nothing to measure",
        "operational_retained_fraction": 1.0,
        "headline_retained_fraction": 1.0,
        "n_removed_quads": 0,
        "core_gate_pass": True,
        "claimed_clean": True,
    })
    rows.append(cert)

    # --- triangle-empty input --------------------------------------------
    seg = "seg_empty"
    cert = base_cert(seg, "PHerc0814", "triangle_empty_invalid")
    cert.update({
        "status": "triangle_empty_invalid",
        "output_mesh": None,
        "output_mesh_hashes": None,
        "wall_seconds": 0.07,
        "input_area_canonical": 0.0,
        "note": ("the input carries no triangles, so no census, no cut and "
                 "no cleanliness claim is defined on it"),
        "evidence": {"n_valid_vertices": 863, "n_retained_quads": 0,
                     "grid_shape": [30, 45], "maxedge": 60.0},
        "operational_retained_fraction": None,
        "headline_retained_fraction": None,
        "core_gate_pass": None,
    })
    rows.append(cert)

    summary_rows = []
    for cert in rows:
        seg = cert["segment"]
        cpath = corpus / f"{seg}_excision_certificate.json"
        cpath.write_text(json.dumps(cert, indent=1))
        (corpus / "logs" / f"{seg}.log").write_text("log\n")
        summary_rows.append({
            "segment": seg,
            "corpus": cert["corpus"],
            "terminal_disposition": cert["terminal_disposition"],
            "status": cert["status"],
            "operational_retained_fraction":
                cert["operational_retained_fraction"],
            "headline_retained_fraction": cert["headline_retained_fraction"],
            "core_gate_pass": cert["core_gate_pass"],
            "wall_seconds": cert["wall_seconds"],
            "certificate": f"out/excised/corpus/{seg}_excision_certificate.json",
            "certificate_sha256": bri.sha256_file(cpath),
            "output_mesh": cert["output_mesh"],
            "base_kind": cert["base_kind"],
            "geometry_key": cert["geometry_key"],
            "is_canonical": True,
            "duplicate_of": None,
            "log": f"out/excised/corpus/logs/{seg}.log",
        })

    (corpus / "corpus_summary.jsonl").write_text(
        "".join(json.dumps(r) + "\n" for r in summary_rows))
    (corpus / "verification.json").write_text(json.dumps({
        "schema": "verify_corpus/v1",
        "census_params": {"cell": 40.0, "diagonals": [0, 1]},
        "independence_note": "fresh workdir per segment",
        "summary": {"n_certificates_checked": 3, "n_meshes_rehashed": 1,
                    "n_meshes_rehash_ok": 1, "n_meshes_rehash_mismatch": 0,
                    "n_recensused": 2,
                    "n_recensus_clean_both_diagonals": 2,
                    "n_census_disagreements": 0,
                    "n_not_censusable_confirmed": 1, "n_failed": 0},
        "segments": [
            {"segment": "seg_transformed",
             "hash_check": {"checked": True, "status": "ok"},
             "recensus": {"kind": "output", "ran": True, "d0_transverse": 0,
                          "d1_transverse": 0, "clean_both_diagonals": True,
                          "comparison": {"transverse_disagrees": False}}},
            {"segment": "seg_clean",
             "hash_check": {"checked": False, "status": "not_applicable"},
             "recensus": {"kind": "input_base", "ran": True,
                          "d0_transverse": 0, "d1_transverse": 0,
                          "clean_both_diagonals": True,
                          "comparison": {"transverse_disagrees": False}}},
            {"segment": "seg_empty",
             "hash_check": {"checked": False, "status": "not_applicable"},
             "recensus": {"ran": False, "not_censusable": True}},
        ],
    }))

    (tmp_path / "out").mkdir(exist_ok=True)
    (tmp_path / "out" / "corpus_bases.json").write_text(json.dumps({
        "schema": "corpus_bases/v1", "n_segments": 3, "n_canonical": 3,
        "n_duplicate_aliases": 0, "n_displacement_repaired": 1,
        "n_original": 2,
        "entries": [{"segment": c["segment"]} for c in rows],
    }))
    return tmp_path


# -------------------------------------------------------------------- tests

def test_dispositions_and_roster(tmp_path):
    index = bri.build(build_fixture(tmp_path))
    assert index["schema"] == "release_index/v1"
    assert len(index["order"]) == 3
    assert index["summary"]["dispositions"] == {
        "transformed": 1, "already_clean": 1, "not_censusable": 1}


def test_non_censusable_is_listed_not_dropped(tmp_path):
    index = bri.build(build_fixture(tmp_path))
    listed = index["summary"]["not_censusable"]
    assert [x["segment"] for x in listed] == ["seg_empty"]
    assert "no triangles" in listed[0]["reason"]
    rec = index["segments"]["seg_empty"]
    assert rec["disposition"] == "not_censusable"
    assert rec["retention"]["operational_retained_fraction"] is None
    assert rec["retention"]["measured"] is False
    assert rec["fragmentation"]["min_R_main_all_components"] is None
    assert rec["fragmentation"]["note"]


def test_both_denominators_carried_and_aggregated(tmp_path):
    index = bri.build(build_fixture(tmp_path))
    rec = index["segments"]["seg_transformed"]["retention"]
    assert rec["operational_A_input_canonical"] == 1000.0
    assert rec["operational_A_excised_canonical"] == 100.0
    assert rec["headline_A_original_canonical"] == 1000.0
    assert rec["headline_A_removed_priced_on_original"] == 100.0

    ret = index["summary"]["retention"]
    # operational: (1000 + 500) cut down by 100 -> 1 - 100/1500
    assert ret["operational_area_weighted"] == 1.0 - 100.0 / 1500.0
    assert ret["operational_segments_counted"] == 2
    # headline: the already-clean segment sits on a repaired base, its
    # original mesh does not exist under tmp_path and no cached area covers
    # it, so it cannot be priced -- and an unpriceable segment is SAID to be
    # outside the denominator rather than silently dropped.
    assert ret["headline_segments_without_area"] == ["seg_clean"]
    assert ret["headline_area_weighted"] == 1.0 - 100.0 / 1000.0
    assert ret["headline_area_weighted_transformed_only"] == 0.9
    # the not-censusable input stays inside the population, at zero area
    assert ret["headline_segments_zero_area"] == ["seg_empty"]
    assert ret["headline_segments_priced"] == 1
    assert ret["headline_segments_counted"] == 2      # 3 pinned, 1 unpriced


def write_area_cache(root: Path, segment: str, area: float) -> Path:
    """The memoised original area, in the shape headline_decision writes.

    A row is trusted only while the segment still names the same original
    mesh with the same hashes, so the fixture's values are repeated here.
    """
    path = root / "out" / "headline_original_areas.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({
        "schema": "headline_original_areas/v1",
        "areas": {segment: {
            "A_original_canonical": area,
            "original_mesh": f"data/{segment}.tifxyz",
            "grid_shape": [10, 10],
            "n_original_retained_quads": 81,
            "maxedge": 60.0,
            "original_hashes": {"x": "aa" * 32, "y": "bb" * 32,
                                "z": "cc" * 32, "mask": None},
        }},
    }, indent=1))
    return path


def test_never_cut_segment_is_priced_from_the_shared_area_cache(tmp_path):
    """An already-clean segment on a repaired base joins the denominator.

    Its certificate records no headline area, because nothing was cut. That
    is a bookkeeping gap, not a missing measurement: the segment enters at
    exactly 1.0 against the full canonical area of its original published
    mesh, read from the same cache bench/headline_decision.py reads.
    """
    root = build_fixture(tmp_path)
    write_area_cache(root, "seg_clean", 600.0)
    index = bri.build(root)

    rec = index["segments"]["seg_clean"]["retention"]
    assert rec["headline_A_original_canonical"] == 600.0
    assert rec["headline_A_removed_priced_on_original"] == 0.0
    assert rec["headline_retained_fraction"] == 1.0
    assert "cached" in rec["headline_area_source"]

    ret = index["summary"]["retention"]
    assert ret["headline_segments_without_area"] == []
    assert ret["headline_area_weighted"] == 1.0 - 100.0 / 1600.0
    assert ret["headline_segments_priced"] == 2
    assert ret["headline_segments_counted"] == 3     # the whole roster
    # the operational denominator is untouched by any of this
    assert ret["operational_area_weighted"] == 1.0 - 100.0 / 1500.0
    assert ret["operational_segments_counted"] == 2

    ap = index["area_pricing"]
    assert ap["n_priced"] == 1 and ap["n_recomputed_from_mesh"] == 0
    assert ap["failures"] == []

    md = bri.render_markdown(index)
    assert "out/headline_original_areas.json" in md
    assert "over 3 pinned artifacts" in md


def test_area_cache_is_not_used_when_the_mesh_no_longer_matches(tmp_path):
    """A stale cache row is ignored, not silently believed."""
    root = build_fixture(tmp_path)
    path = write_area_cache(root, "seg_clean", 600.0)
    doc = json.loads(path.read_text())
    doc["areas"]["seg_clean"]["original_hashes"]["x"] = "99" * 32
    path.write_text(json.dumps(doc, indent=1))

    index = bri.build(root)
    assert index["summary"]["retention"][
        "headline_segments_without_area"] == ["seg_clean"]
    assert index["area_pricing"]["failures"]


def test_percentages_quoted_as_the_decision_rule_prints_them():
    """Truncated toward zero, never rounded up."""
    assert bri.pct_truncated(0.98999999) == "98.999%"
    assert bri.pct_truncated(0.9950546604963569) == "99.505%"
    assert bri.pct_truncated(None) == "not recorded"


def test_already_clean_fragmentation_is_not_invented(tmp_path):
    index = bri.build(build_fixture(tmp_path))
    frag = index["segments"]["seg_clean"]["fragmentation"]
    assert frag["measured"] is False
    assert frag["min_R_main_all_components"] is None
    assert frag["area_weighted_R_main"] is None
    assert "by construction" in frag["note"]
    assert frag["core_gate_pass"] is True


def test_fragmentation_failures_surface(tmp_path):
    index = bri.build(build_fixture(tmp_path))
    frag = index["summary"]["fragmentation"]
    assert frag["n_core_gate_fail"] == 1
    fail = frag["core_gate_failures"][0]
    assert fail["segment"] == "seg_transformed"
    assert fail["min_R_main_core"] == 0.5
    assert fail["area_weighted_R_main"] == 0.8


def test_hash_mapping_and_runtime(tmp_path):
    index = bri.build(build_fixture(tmp_path))
    rec = index["segments"]["seg_transformed"]
    assert rec["input_hashes"]["x"] == "aa" * 32
    assert rec["output_hashes"]["x"] == "11" * 32
    assert index["segments"]["seg_clean"]["output_hashes"] is None
    rt = index["summary"]["runtime_seconds"]
    assert rt["min"] == 0.07 and rt["max"] == 10.0 and rt["n"] == 3
    assert rt["n_over_gate"] == 0


def test_archive_block_measures_real_bytes(tmp_path):
    root = build_fixture(tmp_path)
    index = bri.build(root)
    arc = index["archive"]
    assert arc["total_bytes"] > 0
    assert set(arc["per_scroll"]) == {"Scroll 1", "PHerc0814"}
    assert arc["per_scroll"]["Scroll 1"]["n_meshes"] == 1
    assert arc["per_scroll"]["PHerc0814"]["n_meshes"] == 0
    assert {p["scroll"] for p in arc["proposed_archives"]} == {
        "Scroll 1", "PHerc0814"}
    assert sum(b["bytes"] for b in arc["per_scroll"].values()) \
        == arc["total_bytes"]


def test_markdown_renders(tmp_path):
    index = bri.build(build_fixture(tmp_path))
    md = bri.render_markdown(index)
    assert "represented surface retained" in md
    assert "seg_empty" in md
    assert "usable papyrus" not in md
    for banned in ("text recovered",):
        assert banned not in md
