"""Tests for the round-28 Q5 independent verifier (bench/verify_corpus.py).

Everything here is SYNTHETIC and lives under tmp_path: a fabricated base
manifest, fabricated certificates, and fabricated .tif planes whose only
property that matters is their bytes.  The C++ engine is never invoked --
the census is injected -- so these are fast and pass on a checkout with no
out/ directory at all.

The three things being pinned are the three ways the verifier is allowed to
catch a lie: an emitted mesh whose planes no longer hash to what the
certificate recorded, a certificate with no terminal_disposition, and a
recensus that disagrees with the certificate's recorded census.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
for _p in (REPO_ROOT / "src", REPO_ROOT / "bench"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import verify_corpus as vc                                       # noqa: E402


# ------------------------------------------------------------------ helpers

AXES = ("x", "y", "z")


def write_mesh(root: Path, segment: str, payload: bytes = b"plane") -> Path:
    """A stand-in tifxyz directory: three planes with known bytes."""
    mesh = root / f"{segment}_excised.tifxyz"
    mesh.mkdir(parents=True, exist_ok=True)
    for ax in AXES:
        (mesh / f"{ax}.tif").write_bytes(payload + ax.encode())
    return mesh


def plane_hashes(mesh: Path) -> dict[str, str]:
    return {ax: hashlib.sha256((mesh / f"{ax}.tif").read_bytes()).hexdigest()
            for ax in AXES}


def base_hashes(tag: str = "base") -> dict[str, str]:
    return {ax: hashlib.sha256(f"{tag}-{ax}".encode()).hexdigest()
            for ax in AXES}


def census_block(d0: int, d1: int, triangles: int = 100) -> dict:
    return {f"d{d}": {"triangles": triangles, "quads_dropped": 0,
                      "pairs_tested": 1000, "transverse": t, "coplanar": 0,
                      "grazing": 0}
            for d, t in ((0, d0), (1, d1))}


def make_cert(tmp_path: Path, segment: str, *, disposition="transformed",
              d0: int = 0, d1: int = 0, mesh: Path | None = None,
              output_hashes: dict | None = None,
              bases: dict | None = None,
              drop_terminal_disposition: bool = False) -> Path:
    cert = {
        "record_kind": "excision certificate",
        "segment": segment,
        "terminal_disposition": disposition,
        "terminal_dispositions_defined": list(vc.FALLBACK_DISPOSITIONS),
        "policy_hash": "26296e23cb4d08e4",
        "code_commit": "0" * 40,
        "base_kind": "original",
        "base_mesh": f"data/fake/{segment}.tifxyz",
        "base_hashes": dict(bases or base_hashes()),
        "census_before": census_block(0, 0),
    }
    if disposition == "transformed":
        cert["output_mesh"] = str(mesh)
        cert["output_mesh_hashes"] = output_hashes or plane_hashes(mesh)
        cert["census_after"] = census_block(d0, d1)
    if drop_terminal_disposition:
        cert.pop("terminal_disposition")
    p = tmp_path / f"{segment}{vc.CERT_SUFFIX}"
    p.write_text(json.dumps(cert, indent=1))
    return p


def entry_for(segment: str, bases: dict | None = None) -> dict:
    return {"segment": segment, "base_hashes": dict(bases or base_hashes()),
            "is_canonical": True}


def fake_census(counts: dict):
    """An injected census function returning a fixed engine row."""
    def _fn(mesh, tag):
        return dict(counts)
    return _fn


# -------------------------------------------------------------------- tests

def test_clean_transformed_segment_verifies(tmp_path):
    """The control: honest certificate, matching planes, clean recensus."""
    seg = "20240101000000"
    mesh = write_mesh(tmp_path, seg)
    cert = make_cert(tmp_path, seg, mesh=mesh)
    rec = vc.verify_segment(cert, entry_for(seg),
                            census_fn=fake_census(census_block(0, 0)),
                            repo_root=tmp_path)
    assert rec["problems"] == []
    assert rec["ok"] is True
    assert rec["hash_check"]["status"] == "ok"
    assert rec["recensus"]["clean_both_diagonals"] is True
    assert rec["recensus"]["comparison"]["transverse_disagrees"] is False


def test_output_hash_mismatch_is_detected(tmp_path):
    """A plane rewritten after the certificate was signed must be caught."""
    seg = "20240102000000"
    mesh = write_mesh(tmp_path, seg)
    cert = make_cert(tmp_path, seg, mesh=mesh)
    # tamper with the emitted mesh AFTER the certificate recorded its hashes
    (mesh / "y.tif").write_bytes(b"tampered")

    rec = vc.verify_segment(cert, entry_for(seg),
                            census_fn=fake_census(census_block(0, 0)),
                            repo_root=tmp_path)
    assert rec["ok"] is False
    assert rec["hash_check"]["status"] == "mismatch"
    assert rec["hash_check"]["mismatched_axes"] == ["y"]
    assert rec["hash_check"]["per_axis"]["x"]["match"] is True
    assert rec["hash_check"]["per_axis"]["z"]["match"] is True
    assert any("hash mismatch" in p for p in rec["problems"])

    summary = vc.summarise([rec], {"n_manifest_entries": 1})
    assert summary["n_meshes_rehash_mismatch"] == 1
    assert summary["hash_mismatch_segments"] == [seg]


def test_missing_terminal_disposition_is_flagged(tmp_path):
    """A certificate with no terminal_disposition is not a result."""
    seg = "20240103000000"
    mesh = write_mesh(tmp_path, seg)
    cert = make_cert(tmp_path, seg, mesh=mesh,
                     drop_terminal_disposition=True)
    rec = vc.verify_segment(cert, entry_for(seg),
                            census_fn=fake_census(census_block(0, 0)),
                            repo_root=tmp_path)
    assert rec["ok"] is False
    assert rec["disposition"] is None
    assert "no terminal_disposition" in rec["certificate_problems"]
    assert "no terminal_disposition" in rec["problems"]


def test_recorded_vs_recomputed_census_disagreement_is_flagged(tmp_path):
    """The certificate says 0/0; the independent recensus says otherwise."""
    seg = "20240104000000"
    mesh = write_mesh(tmp_path, seg)
    cert = make_cert(tmp_path, seg, mesh=mesh, d0=0, d1=0)
    rec = vc.verify_segment(cert, entry_for(seg),
                            census_fn=fake_census(census_block(2, 5)),
                            repo_root=tmp_path)
    assert rec["ok"] is False
    assert rec["recensus"]["d0_transverse"] == 2
    assert rec["recensus"]["d1_transverse"] == 5
    assert rec["recensus"]["clean_both_diagonals"] is False
    cmp_ = rec["recensus"]["comparison"]
    assert cmp_["transverse_disagrees"] is True
    assert {(f["diagonal"], f["recomputed"], f["recorded"])
            for f in cmp_["differing_fields"] if f["field"] == "transverse"} \
        == {(0, 2, 0), (1, 5, 0)}
    assert any("recensus NOT clean" in p for p in rec["problems"])
    assert any("disagrees with the certificate" in p for p in rec["problems"])

    summary = vc.summarise([rec], {"n_manifest_entries": 1})
    assert summary["n_recensus_clean_both_diagonals"] == 0
    assert summary["n_census_disagreements"] == 1
    assert summary["census_disagreement_segments"] == [seg]


def test_stale_base_hashes_are_flagged(tmp_path):
    """The manifest moved under the certificate: the cut is not replayable."""
    seg = "20240105000000"
    mesh = write_mesh(tmp_path, seg)
    cert = make_cert(tmp_path, seg, mesh=mesh, bases=base_hashes("old"))
    rec = vc.verify_segment(cert, entry_for(seg, base_hashes("new")),
                            census_fn=fake_census(census_block(0, 0)),
                            repo_root=tmp_path)
    assert rec["ok"] is False
    assert any(p.startswith("base hash ") for p in rec["certificate_problems"])


def test_already_clean_input_must_recensus_clean(tmp_path):
    """already_clean claims the INPUT is 0/0; a dirty input breaks it."""
    seg = "20240106000000"
    base = write_mesh(tmp_path, seg)
    cert_path = tmp_path / f"{seg}{vc.CERT_SUFFIX}"
    cert_path.write_text(json.dumps({
        "record_kind": "excision certificate",
        "segment": seg,
        "terminal_disposition": "already_clean",
        "terminal_dispositions_defined": list(vc.FALLBACK_DISPOSITIONS),
        "policy_hash": "26296e23cb4d08e4",
        "code_commit": "0" * 40,
        "base_mesh": str(base),
        "base_hashes": base_hashes(),
        "census_before": census_block(0, 0),
    }, indent=1))

    good = vc.verify_segment(cert_path, entry_for(seg),
                             census_fn=fake_census(census_block(0, 0)),
                             repo_root=tmp_path)
    assert good["ok"] is True
    assert good["recensus"]["kind"] == "input_base"

    bad = vc.verify_segment(cert_path, entry_for(seg),
                            census_fn=fake_census(census_block(0, 3)),
                            repo_root=tmp_path)
    assert bad["ok"] is False
    assert any("recensus NOT clean" in p for p in bad["problems"])


def test_triangle_empty_invalid_must_not_be_censusable(tmp_path):
    """The claim is that no census is defined; a census that succeeds
    contradicts it."""
    seg = "20240107000000"
    base = write_mesh(tmp_path, seg)
    body = {
        "record_kind": "excision certificate",
        "segment": seg,
        "terminal_disposition": "triangle_empty_invalid",
        "terminal_dispositions_defined": list(vc.FALLBACK_DISPOSITIONS),
        "policy_hash": "26296e23cb4d08e4",
        "code_commit": "0" * 40,
        "base_mesh": str(base),
        "base_hashes": base_hashes(),
    }
    p = tmp_path / f"{seg}{vc.CERT_SUFFIX}"
    p.write_text(json.dumps(body, indent=1))

    # engine declines the mesh (census_one returns None) -> claim upheld
    ok = vc.verify_segment(p, entry_for(seg), census_fn=lambda m, t: None,
                           repo_root=tmp_path)
    assert ok["ok"] is True
    assert ok["recensus"]["not_censusable"] is True

    # engine censuses it after all -> the claim was false
    bad = vc.verify_segment(p, entry_for(seg),
                            census_fn=fake_census(census_block(0, 0)),
                            repo_root=tmp_path)
    assert bad["ok"] is False
    assert any("triangle-empty/invalid" in x for x in bad["problems"])


def test_missing_output_mesh_is_a_failure(tmp_path):
    """A transformed certificate whose mesh is not on disk fails."""
    seg = "20240108000000"
    mesh = write_mesh(tmp_path, seg)
    cert = make_cert(tmp_path, seg, mesh=mesh)
    for ax in AXES:
        (mesh / f"{ax}.tif").unlink()
    mesh.rmdir()
    rec = vc.verify_segment(cert, entry_for(seg),
                            census_fn=fake_census(census_block(0, 0)),
                            repo_root=tmp_path)
    assert rec["ok"] is False
    assert rec["hash_check"]["status"] == "mesh_missing"
    assert any("mesh_missing" in p for p in rec["problems"])
    assert any("mesh directory missing" in p for p in rec["problems"])


def test_roster_gaps_both_directions(tmp_path):
    entries = [{"segment": "a"}, {"segment": "b"}]
    certs = [tmp_path / f"b{vc.CERT_SUFFIX}", tmp_path / f"c{vc.CERT_SUFFIX}"]
    roster = vc.build_roster(entries, certs)
    assert roster["manifest_entries_without_certificate"] == ["a"]
    assert roster["certificates_not_on_manifest"] == ["c"]


def test_normalise_hash_keys_accepts_both_spellings():
    assert vc.normalise_hash_keys({"x.tif": "1", "y": "2"}) == {"x": "1",
                                                                "y": "2"}
