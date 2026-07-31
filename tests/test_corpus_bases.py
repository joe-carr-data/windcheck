"""Tests for the round-28 base manifest builder (bench/corpus_bases.py).

Everything here runs on synthetic tifxyz directories under tmp_path -- the
real 185-segment corpus is never touched, so these pass on a checkout with
no data/ at all.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "bench"))

import corpus_bases as cb                                       # noqa: E402


# ------------------------------------------------------------------ helpers

def make_mesh(root: Path, name: str, payload: bytes) -> Path:
    """A minimal tifxyz directory whose three planes derive from payload."""
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    for axis in ("x", "y", "z"):
        (d / f"{axis}.tif").write_bytes(payload + axis.encode())
    (d / "meta.json").write_text('{"synthetic": true}')
    return d


def sha_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def write_cert(cert_root: Path, segment: str, out_hashes: dict) -> Path:
    cert_root.mkdir(parents=True, exist_ok=True)
    p = cert_root / f"{segment}_multi_certificate.json"
    p.write_text(json.dumps({
        "segment": segment,
        "hashes": {"input": {}, "output": out_hashes},
    }))
    return p


def row(segment: str, mesh: Path, corpus: str = "Synthetic",
        volume: str = "20240101000000") -> dict:
    return {"segment": segment, "corpus": corpus, "volume": volume,
            "mesh": mesh}


# --------------------------------------------------------------- voxel scale

def test_voxel_scale_parses_resolution_token():
    assert cb.voxel_scale(
        "20230702185753-on-20230205180739-7.91um.tifxyz") == 7.91
    assert cb.voxel_scale("seg-on-vol-15.4um.tifxyz") == 15.4
    assert cb.voxel_scale("seg-on-vol-8um.tifxyz") == 8.0


def test_voxel_scale_defaults_when_token_absent():
    # Same fallback bench/excise_shadow.py uses: RES_UM miss -> 7.91.
    assert cb.voxel_scale("no_resolution_here.tifxyz") == cb.DEFAULT_VOXEL_UM
    assert cb.voxel_scale("no_resolution_here.tifxyz") == 7.91
    # The token must be anchored to the tifxyz suffix, not matched anywhere.
    assert cb.voxel_scale("7.91um-in-the-middle.tifxyz") == 7.91
    assert cb.voxel_scale("mesh-3.24um.tifxyz.bak") == 7.91


# ------------------------------------------------------------- mesh hashing

def test_mesh_hashes_and_geometry_key(tmp_path):
    a = make_mesh(tmp_path, "a.tifxyz", b"same")
    b = make_mesh(tmp_path, "b.tifxyz", b"same")
    c = make_mesh(tmp_path, "c.tifxyz", b"different")

    ha, hb, hc = (cb.mesh_hashes(p) for p in (a, b, c))
    assert ha == hb and ha != hc
    assert ha["mask"] is None                     # no mask.tif written
    assert ha["x"] == sha_bytes(b"samex")
    assert cb.geometry_key(ha) == cb.geometry_key(hb)
    assert cb.geometry_key(ha) != cb.geometry_key(hc)

    # mask presence must not be allowed to change the geometry key silently
    (a / "mask.tif").write_bytes(b"mask")
    ha2 = cb.mesh_hashes(a)
    assert ha2["mask"] == sha_bytes(b"mask")
    assert cb.geometry_key(ha2) == cb.geometry_key(hb)


# ------------------------------------------------- repaired-mesh verification

def test_repaired_mesh_matching_certificate_is_used(tmp_path):
    seg = "seg_ok"
    orig = make_mesh(tmp_path / "orig", f"{seg}-on-v-7.91um.tifxyz", b"orig")
    rep_root = tmp_path / "meshes"
    rep = make_mesh(rep_root, f"{seg}_repaired.tifxyz", b"repaired")
    disk = cb.mesh_hashes(rep)
    write_cert(tmp_path / "cert", seg,
               {f"{a}.tif": disk[a] for a in ("x", "y", "z")})

    rec = cb.classify_base(seg, orig, cb.mesh_hashes(orig),
                           repaired_root=rep_root,
                           cert_root=tmp_path / "cert")

    assert rec["base_kind"] == "displacement_repaired"
    assert Path(rec["base_mesh"]) == rep
    assert rec["base_hashes"] == disk
    assert rec["repair_cert_output_hashes_verified"] is True
    assert rec["repair_certificate_sha256"] is not None
    assert "match" in rec["repair_verification_note"]


def test_repaired_mesh_hash_mismatch_falls_back_to_original(tmp_path):
    seg = "seg_bad"
    orig = make_mesh(tmp_path / "orig", f"{seg}-on-v-7.91um.tifxyz", b"orig")
    rep_root = tmp_path / "meshes"
    make_mesh(rep_root, f"{seg}_repaired.tifxyz", b"repaired")
    # Certificate records hashes of some OTHER content: a stale/tampered mesh.
    write_cert(tmp_path / "cert", seg,
               {f"{a}.tif": sha_bytes(b"somethingelse" + a.encode())
                for a in ("x", "y", "z")})

    orig_h = cb.mesh_hashes(orig)
    rec = cb.classify_base(seg, orig, orig_h, repaired_root=rep_root,
                           cert_root=tmp_path / "cert")

    assert rec["base_kind"] == "original"
    assert Path(rec["base_mesh"]) == orig
    assert rec["base_hashes"] == orig_h
    assert rec["repair_cert_output_hashes_verified"] is False
    assert rec["repair_verification_note"]
    assert "disagrees" in rec["repair_verification_note"]
    # partial mismatch (x/y right, z wrong) is still a rejection
    disk = cb.mesh_hashes(rep_root / f"{seg}_repaired.tifxyz")
    write_cert(tmp_path / "cert", seg,
               {"x.tif": disk["x"], "y.tif": disk["y"],
                "z.tif": sha_bytes(b"wrong")})
    rec2 = cb.classify_base(seg, orig, orig_h, repaired_root=rep_root,
                            cert_root=tmp_path / "cert")
    assert rec2["base_kind"] == "original"
    assert rec2["repair_cert_output_hashes_verified"] is False
    assert "z.tif" in rec2["repair_verification_note"]


def test_repaired_mesh_without_certificate_is_refused(tmp_path):
    seg = "seg_nocert"
    orig = make_mesh(tmp_path / "orig", f"{seg}-on-v-7.91um.tifxyz", b"orig")
    rep_root = tmp_path / "meshes"
    make_mesh(rep_root, f"{seg}_repaired.tifxyz", b"repaired")

    rec = cb.classify_base(seg, orig, cb.mesh_hashes(orig),
                           repaired_root=rep_root,
                           cert_root=tmp_path / "cert")

    assert rec["base_kind"] == "original"
    assert rec["repair_cert_output_hashes_verified"] is False
    assert rec["repair_certificate"] is None
    assert "no repair certificate" in rec["repair_verification_note"]


def test_certificate_without_output_hashes_is_refused(tmp_path):
    seg = "seg_nohash"
    orig = make_mesh(tmp_path / "orig", f"{seg}-on-v-7.91um.tifxyz", b"orig")
    rep_root = tmp_path / "meshes"
    make_mesh(rep_root, f"{seg}_repaired.tifxyz", b"repaired")
    cert_root = tmp_path / "cert"
    cert_root.mkdir()
    (cert_root / f"{seg}_multi_certificate.json").write_text(
        json.dumps({"segment": seg, "hashes": {"input": {"x.tif": "aa"}}}))

    rec = cb.classify_base(seg, orig, cb.mesh_hashes(orig),
                           repaired_root=rep_root, cert_root=cert_root)

    assert rec["base_kind"] == "original"
    assert rec["repair_cert_output_hashes_verified"] is False
    assert "hashes.output" in rec["repair_verification_note"]


def test_no_repaired_mesh_is_original_not_a_failure(tmp_path):
    seg = "seg_plain"
    orig = make_mesh(tmp_path / "orig", f"{seg}-on-v-7.91um.tifxyz", b"orig")

    rec = cb.classify_base(seg, orig, cb.mesh_hashes(orig),
                           repaired_root=tmp_path / "meshes",
                           cert_root=tmp_path / "cert")

    assert rec["base_kind"] == "original"
    assert rec["repair_cert_output_hashes_verified"] is None
    assert rec["repair_certificate"] is None
    assert "no displacement-repaired mesh" in rec["repair_verification_note"]


# ----------------------------------------------------------------- duplicates

def test_duplicate_aliasing_picks_alphabetically_first_canonical(tmp_path):
    # bbb and aaa share geometry; ccc is unique.
    orig = tmp_path / "orig"
    rows = [
        row("bbb", make_mesh(orig, "bbb-on-v-7.91um.tifxyz", b"shared")),
        row("aaa", make_mesh(orig, "aaa-on-v-7.91um.tifxyz", b"shared")),
        row("ccc", make_mesh(orig, "ccc-on-v-7.91um.tifxyz", b"unique")),
    ]
    doc = cb.build(rows, manifests={}, repaired_root=tmp_path / "meshes",
                   cert_root=tmp_path / "cert", code_commit="TEST")
    by = {e["segment"]: e for e in doc["entries"]}

    assert [e["segment"] for e in doc["entries"]] == ["aaa", "bbb", "ccc"]
    assert by["aaa"]["is_canonical"] is True
    assert by["aaa"]["duplicate_of"] is None
    assert by["bbb"]["is_canonical"] is False
    assert by["bbb"]["duplicate_of"] == "aaa"
    assert by["ccc"]["is_canonical"] is True
    assert by["ccc"]["duplicate_of"] is None

    assert doc["n_segments"] == 3
    assert doc["n_canonical"] == 2
    assert doc["n_duplicate_aliases"] == 1
    assert cb.duplicate_groups(doc["entries"]) == [
        (by["aaa"]["geometry_key"], ["aaa", "bbb"])]


def test_duplicate_grouping_is_on_the_base_not_the_original(tmp_path):
    """Two segments with identical ORIGINALS diverge once one of them gets a
    verified repaired base -- the aliasing must follow the base geometry."""
    orig = tmp_path / "orig"
    rep_root, cert_root = tmp_path / "meshes", tmp_path / "cert"
    rows = [
        row("aaa", make_mesh(orig, "aaa-on-v-7.91um.tifxyz", b"shared")),
        row("bbb", make_mesh(orig, "bbb-on-v-7.91um.tifxyz", b"shared")),
    ]
    rep = make_mesh(rep_root, "bbb_repaired.tifxyz", b"moved")
    d = cb.mesh_hashes(rep)
    write_cert(cert_root, "bbb", {f"{a}.tif": d[a] for a in ("x", "y", "z")})

    doc = cb.build(rows, manifests={}, repaired_root=rep_root,
                   cert_root=cert_root, code_commit="TEST")
    by = {e["segment"]: e for e in doc["entries"]}

    assert by["bbb"]["base_kind"] == "displacement_repaired"
    assert by["aaa"]["base_kind"] == "original"
    assert doc["n_canonical"] == 2 and doc["n_duplicate_aliases"] == 0
    assert all(e["is_canonical"] for e in doc["entries"])
    # ... but the ORIGINAL duplication is still visible to a reader
    assert by["aaa"]["original_geometry_key"] == \
        by["bbb"]["original_geometry_key"]
    assert cb.duplicate_groups(doc["entries"],
                               key="original_geometry_key") == [
        (by["aaa"]["original_geometry_key"], ["aaa", "bbb"])]


def test_assign_duplicates_three_way_group():
    entries = [{"segment": s, "geometry_key": k} for s, k in
               [("zed", "K"), ("mid", "K"), ("abe", "K"), ("solo", "J")]]
    cb.assign_duplicates(entries)
    by = {e["segment"]: e for e in entries}
    assert by["abe"]["is_canonical"] and by["abe"]["duplicate_of"] is None
    for s in ("mid", "zed"):
        assert by[s]["is_canonical"] is False
        assert by[s]["duplicate_of"] == "abe"
    assert by["solo"]["is_canonical"] and by["solo"]["duplicate_of"] is None


# ---------------------------------------------------------- published manifest

def test_published_manifest_lookup_and_agreement(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    orig = make_mesh(tmp_path / "orig", "seg1-on-v-7.91um.tifxyz", b"orig")
    disk = cb.mesh_hashes(orig)
    prefix = "seg1/mesh/seg1-on-v-7.91um.tifxyz"
    (data / "MANIFEST-Synth.json").write_text(json.dumps({
        "bucket": "b", "sample": "Synth", "volume": "V1",
        "files": [{"path": f"{prefix}/{a}.tif", "sha256": disk[a]}
                  for a in ("x", "y", "z")]}))
    index = cb.load_manifests(data)

    pub, src = cb.published_lookup(index, "V1", "seg1",
                                   "seg1-on-v-7.91um.tifxyz")
    assert src.endswith("MANIFEST-Synth.json")
    assert pub["x"] == disk["x"] and pub["mask"] is None
    assert cb.manifest_agrees(pub, disk) is True

    # a segment the manifest does not list
    assert cb.published_lookup(index, "V1", "other", "m.tifxyz") == (None, None)
    assert cb.manifest_agrees(None, disk) is None
    # disagreement is reported, not swallowed
    bad = dict(disk, x="0" * 64)
    assert cb.manifest_agrees(pub, bad) is False


# --------------------------------------------------------------- determinism

def test_build_is_deterministic(tmp_path):
    orig = tmp_path / "orig"
    rep_root, cert_root = tmp_path / "meshes", tmp_path / "cert"
    rows = [
        row("aaa", make_mesh(orig, "aaa-on-v-7.91um.tifxyz", b"shared")),
        row("bbb", make_mesh(orig, "bbb-on-v-7.91um.tifxyz", b"shared")),
        row("ccc", make_mesh(orig, "ccc-on-v-15.4um.tifxyz", b"unique")),
    ]
    rep = make_mesh(rep_root, "ccc_repaired.tifxyz", b"fixed")
    d = cb.mesh_hashes(rep)
    write_cert(cert_root, "ccc", {f"{a}.tif": d[a] for a in ("x", "y", "z")})

    d1 = cb.build(rows, manifests={}, repaired_root=rep_root,
                  cert_root=cert_root, code_commit="TEST")
    d2 = cb.build(list(reversed(rows)), manifests={}, repaired_root=rep_root,
                  cert_root=cert_root, code_commit="TEST")

    assert json.dumps(d1["entries"]) == json.dumps(d2["entries"])
    assert cb.stable(d1) == cb.stable(d2)
    assert cb.diff_docs(d1, d2) == []
    # generated_utc is the only wall-clock field and it is excluded from
    # the drift comparison
    d3 = dict(d1, generated_utc="1970-01-01T00:00:00Z")
    assert cb.diff_docs(d1, d3) == []
    assert d1["entries"][2]["voxel_um"] == 15.4


def test_diff_docs_reports_real_drift(tmp_path):
    orig = tmp_path / "orig"
    rows = [row("aaa", make_mesh(orig, "aaa-on-v-7.91um.tifxyz", b"one"))]
    d1 = cb.build(rows, manifests={}, repaired_root=tmp_path / "m",
                  cert_root=tmp_path / "c", code_commit="TEST")
    (orig / "aaa-on-v-7.91um.tifxyz" / "x.tif").write_bytes(b"two")
    d2 = cb.build(rows, manifests={}, repaired_root=tmp_path / "m",
                  cert_root=tmp_path / "c", code_commit="TEST")
    msgs = cb.diff_docs(d1, d2)
    assert msgs
    assert any("aaa.geometry_key" in m for m in msgs)


# --------------------------------------------------------------- entry schema

REQUIRED_ENTRY_FIELDS = [
    "segment", "corpus", "volume", "voxel_um", "original_mesh",
    "original_hashes", "published_manifest_hashes", "published_manifest_source",
    "published_manifest_agrees_with_disk", "base_kind", "base_mesh",
    "base_hashes", "repair_certificate", "repair_certificate_sha256",
    "repair_cert_output_hashes_verified", "repair_verification_note",
    "geometry_key", "original_geometry_key", "duplicate_of", "is_canonical",
]


def test_entry_schema_and_top_level_keys(tmp_path):
    orig = tmp_path / "orig"
    rows = [row("aaa", make_mesh(orig, "aaa-on-v-7.91um.tifxyz", b"one"))]
    doc = cb.build(rows, manifests={}, repaired_root=tmp_path / "m",
                   cert_root=tmp_path / "c", code_commit="TEST")
    assert doc["schema"] == "corpus_bases/v1"
    # provenance a public reader can recompute; never a commit sha
    assert len(doc["provenance"]["source_tree_digest"]) == 64
    assert "code_commit" not in doc
    for k in ("schema", "generated_utc", "provenance", "n_segments",
              "n_canonical", "n_duplicate_aliases", "n_displacement_repaired",
              "n_original", "counts_note", "entries"):
        assert k in doc
    e = doc["entries"][0]
    assert list(e) == REQUIRED_ENTRY_FIELDS
    assert e["base_kind"] in ("original", "displacement_repaired")
    assert isinstance(e["voxel_um"], float)
    assert isinstance(e["base_mesh"], str)
    # the document must round-trip through json unchanged
    assert json.loads(cb.dumps(doc)) == doc


def test_repo_path_is_stable_for_absolute_paths(tmp_path):
    # paths outside the repo stay absolute rather than raising
    p = tmp_path / "x.tif"
    p.write_text("x")
    assert cb.repo_path(p) == p.as_posix()
    assert cb.repo_path(REPO_ROOT / "bench" / "corpus_bases.py") == \
        "bench/corpus_bases.py"
