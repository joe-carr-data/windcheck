"""Release provenance, and the index's refusal to publish a contradiction.

Two obligations:

  * provenance a public consumer can recompute from the published files --
    no commit sha anywhere in a published artifact;
  * index generation FAILS LOUDLY when a certificate's declared base_kind
    disagrees with the hash-verified answer, rather than quietly publishing
    the declaration.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "bench"))

from windcheck.manifest import BaseKindDisagreement            # noqa: E402
from windcheck.provenance import (release_provenance,          # noqa: E402
                                  source_files, source_tree_digest)

import build_release_index as bri                              # noqa: E402


# ------------------------------------------------------------- provenance

def test_source_tree_digest_is_stable_and_covers_the_code():
    assert source_tree_digest(REPO_ROOT) == source_tree_digest(REPO_ROOT)
    files = source_files(REPO_ROOT)
    assert "src/windcheck/manifest.py" in files
    assert "src/windcheck/provenance.py" in files
    assert "uv.lock" in files
    assert not any("__pycache__" in f for f in files)


def test_source_tree_digest_moves_when_a_source_file_moves(tmp_path):
    root = tmp_path / "tree"
    (root / "src" / "windcheck").mkdir(parents=True)
    (root / "pyproject.toml").write_text("[project]\n")
    f = root / "src" / "windcheck" / "excise.py"
    f.write_text("A = 1\n")
    before = source_tree_digest(root)
    f.write_text("A = 2\n")
    assert source_tree_digest(root) != before


def test_provenance_block_cites_no_repository_history():
    block = release_provenance(REPO_ROOT)
    assert block["code_version"]
    assert len(block["source_tree_digest"]) == 64
    blob = json.dumps(block).lower()
    for banned in ("commit", "branch", "git"):
        assert banned not in blob, f"provenance still cites {banned}"


def test_published_index_carries_no_commit_sha():
    index_path = REPO_ROOT / "out" / "release" / "index.json"
    if not index_path.exists():                       # not built in this tree
        pytest.skip("release index not present")
    doc = json.loads(index_path.read_text())
    assert "code_commit" not in doc
    prov = doc["provenance"]
    assert len(prov["source_tree_digest"]) == 64
    assert prov["code_version"]


def test_published_certificates_carry_no_commit_sha():
    corpus = REPO_ROOT / "out" / "excised" / "corpus"
    certs = sorted(corpus.glob("*_excision_certificate.json"))
    if not certs:
        pytest.skip("corpus certificates not present")
    for path in certs:
        doc = json.loads(path.read_text())
        assert "code_commit" not in doc, path.name
        prov = doc.get("code_provenance") or {}
        for banned in ("commit", "branch", "git_status_porcelain"):
            assert banned not in prov, f"{path.name} still cites {banned}"
        # what a reader CAN check is still there
        assert prov.get("engine_source_sha256")
        assert doc.get("policy_hash")


# ------------------------------------------------- base-kind disagreement

def _cert(base_kind: str, base_hashes: dict, original_hashes: dict) -> dict:
    return {"segment": "seg1", "base_kind": base_kind,
            "base_hashes": base_hashes,
            "original_hashes": original_hashes,
            "terminal_disposition": "transformed",
            "certificate": "out/excised/corpus/seg1_excision_certificate.json"}


ORIG = {"x": "a" * 64, "y": "b" * 64, "z": "c" * 64, "mask": None}
REPAIRED = {"x": "d" * 64, "y": "b" * 64, "z": "c" * 64, "mask": None}


def test_base_kind_check_agrees_on_a_truthful_certificate():
    chk = bri.base_kind_check(_cert("original", ORIG, ORIG), "seg1")
    assert chk["agrees"] and chk["verified"] == "original"
    chk = bri.base_kind_check(_cert("displacement_repaired", REPAIRED, ORIG),
                              "seg1")
    assert chk["agrees"] and chk["verified"] == "displacement_repaired"


def test_base_kind_check_catches_a_lying_certificate():
    chk = bri.base_kind_check(_cert("original", REPAIRED, ORIG), "seg1")
    assert chk["agrees"] is False
    assert chk["verified"] == "displacement_repaired"
    assert "seg1" in chk["reason"]


def test_legacy_original_published_spelling_still_agrees():
    chk = bri.base_kind_check(_cert("original_published", ORIG, ORIG), "seg1")
    assert chk["agrees"] and chk["verified"] == "original"


def build_fixture(tmp_path: Path, base_kind: str, base_hashes: dict) -> Path:
    """A minimal repo layout `build_release_index.build` can read."""
    corpus = tmp_path / "out" / "excised" / "corpus"
    corpus.mkdir(parents=True)
    (tmp_path / "out" / "release").mkdir(parents=True)
    cert_rel = "out/excised/corpus/seg1_excision_certificate.json"
    cert = {
        "segment": "seg1", "terminal_disposition": "transformed",
        "record_kind": "excision certificate",
        "base_kind": base_kind, "base_hashes": base_hashes,
        "original_hashes": ORIG, "base_mesh": "data/seg1",
        "original_mesh": "data/seg1", "policy_hash": "p" * 64,
        "output_mesh": None, "output_mesh_hashes": None,
    }
    (tmp_path / cert_rel).write_text(json.dumps(cert))
    (corpus / "corpus_summary.jsonl").write_text(json.dumps(
        {"segment": "seg1", "corpus": "Scroll 1", "certificate": cert_rel,
         "terminal_disposition": "transformed"}) + "\n")
    (tmp_path / "out" / "corpus_bases.json").write_text(json.dumps({
        "schema": "corpus_bases/v1",
        "entries": [{"segment": "seg1", "base_kind": base_kind,
                     "base_hashes": base_hashes, "original_hashes": ORIG,
                     "is_canonical": True, "duplicate_of": None}]}))
    return tmp_path


def test_index_build_fails_loudly_on_a_base_kind_disagreement(tmp_path):
    root = build_fixture(tmp_path, "original", REPAIRED)
    with pytest.raises(BaseKindDisagreement) as exc:
        bri.build(root)
    msg = str(exc.value)
    assert "seg1" in msg
    assert "displacement_repaired" in msg


def test_index_build_accepts_an_honest_declaration(tmp_path):
    root = build_fixture(tmp_path, "displacement_repaired", REPAIRED)
    index = bri.build(root)
    rec = index["segments"]["seg1"]
    assert rec["base_kind"] == "displacement_repaired"
    assert rec["base_kind_agrees"] is True
    assert "code_commit" not in index
    assert index["provenance"]["source_tree_digest"]
