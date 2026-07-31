"""The mesh content manifest, its digest, and the one shared base-kind rule.

Three things the release depends on and that a regression would hide:

  * the manifest and its serialisation are stable — same bytes, same digest,
    independent of filesystem enumeration order;
  * changing `mask.tif` changes the directory digest, because the mask is
    what decides which triangles exist and it used to sit outside the hash;
  * a declared-but-absent file is represented explicitly, not skipped, so
    "no mask" and "a mask that was never looked at" cannot collide.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from windcheck.manifest import (ABSENT, SEMANTIC_FILES,          # noqa: E402
                                BaseKind, BaseKindDisagreement,
                                base_kind_from_manifests, digest, file_row,
                                mesh_manifest, normalise_base_kind,
                                require_base_kind, same_mesh, serialise,
                                verify_base_kind)


def make_mesh(root: Path, name: str, *, mask: bytes | None = None,
              meta: dict | None = None, planes: bytes = b"XYZ") -> Path:
    d = root / name
    d.mkdir(parents=True)
    for axis in ("x", "y", "z"):
        (d / f"{axis}.tif").write_bytes(planes + axis.encode())
    if mask is not None:
        (d / "mask.tif").write_bytes(mask)
    if meta is not None:
        (d / "meta.json").write_text(json.dumps(meta))
    return d


# ---------------------------------------------------------------- stability

def test_manifest_is_stable_across_calls(tmp_path):
    mesh = make_mesh(tmp_path, "m", mask=b"\x01\x00", meta={"shape": [2, 2]})
    a, b = mesh_manifest(mesh), mesh_manifest(mesh)
    assert a == b
    assert a["digest"] == b["digest"]


def test_manifest_covers_every_semantic_file(tmp_path):
    mesh = make_mesh(tmp_path, "m", mask=b"\x01", meta={"shape": [1, 1]})
    paths = [r["path"] for r in mesh_manifest(mesh)["files"]]
    assert set(SEMANTIC_FILES) <= set(paths)
    assert paths == sorted(paths)


def test_digest_is_independent_of_row_order(tmp_path):
    mesh = make_mesh(tmp_path, "m", mask=b"\x01", meta={"shape": [1, 1]})
    rows = mesh_manifest(mesh)["files"]
    assert digest(rows) == digest(list(reversed(rows)))


def test_serialisation_is_the_documented_one(tmp_path):
    mesh = make_mesh(tmp_path, "m", meta={"shape": [1, 1]})
    man = mesh_manifest(mesh)
    blob = serialise(man["files"]).decode("utf-8")
    lines = blob.splitlines()
    assert blob.endswith("\n")
    assert len(lines) == len(man["files"])
    for line, row in zip(lines, sorted(man["files"],
                                       key=lambda r: r["path"])):
        path, size, sha = line.split("\0")
        assert path == row["path"]
        if row["present"]:
            assert size == str(row["size"]) and sha == row["sha256"]
        else:
            assert size == ABSENT and sha == ABSENT
    assert man["digest"] == hashlib.sha256(
        blob.encode("utf-8")).hexdigest()


def test_manifest_carries_the_legacy_plane_view(tmp_path):
    mesh = make_mesh(tmp_path, "m", mask=b"\x01")
    man = mesh_manifest(mesh)
    for axis in ("x", "y", "z"):
        assert man[axis] == hashlib.sha256(
            (mesh / f"{axis}.tif").read_bytes()).hexdigest()
    assert man["mask"] == hashlib.sha256(b"\x01").hexdigest()


# -------------------------------------------------------- the mask matters

def test_changing_mask_changes_the_directory_digest(tmp_path):
    a = make_mesh(tmp_path, "a", mask=b"\x01\x01", meta={"shape": [2, 1]})
    before = mesh_manifest(a)
    (a / "mask.tif").write_bytes(b"\x01\x00")          # excise one cell
    after = mesh_manifest(a)
    assert before["digest"] != after["digest"], (
        "a mask edit must invalidate the recorded mesh digest: the mask is "
        "what decides which triangles exist")
    assert before["x"] == after["x"], "coordinates did not move"
    assert same_mesh(before, after) is False


def test_changing_meta_changes_the_directory_digest(tmp_path):
    a = make_mesh(tmp_path, "a", meta={"shape": [2, 1]})
    before = mesh_manifest(a)
    (a / "meta.json").write_text(json.dumps({"shape": [1, 2]}))
    assert mesh_manifest(a)["digest"] != before["digest"]


def test_adding_a_mask_changes_the_digest(tmp_path):
    a = make_mesh(tmp_path, "a")
    before = mesh_manifest(a)
    (a / "mask.tif").write_bytes(b"\x01")
    assert mesh_manifest(a)["digest"] != before["digest"]


# ------------------------------------------------ absent files are explicit

def test_absent_optional_file_is_represented_explicitly(tmp_path):
    mesh = make_mesh(tmp_path, "m")                    # no mask, no meta
    man = mesh_manifest(mesh)
    rows = {r["path"]: r for r in man["files"]}
    for name in ("mask.tif", "mask.png", "meta.json"):
        assert name in rows, f"{name} was silently skipped"
        assert rows[name] == {"path": name, "present": False,
                              "size": None, "sha256": None}
    blob = serialise(man["files"]).decode("utf-8")
    assert f"mask.tif\0{ABSENT}\0{ABSENT}" in blob


def test_absent_row_is_not_the_same_as_no_row(tmp_path):
    mesh = make_mesh(tmp_path, "m")
    rows = mesh_manifest(mesh)["files"]
    without = [r for r in rows if r["present"]]
    assert digest(rows) != digest(without), (
        "dropping the absent rows must change the digest, otherwise an "
        "undeclared file is indistinguishable from an absent one")


def test_file_row_of_a_missing_file(tmp_path):
    assert file_row(tmp_path, "nope.tif") == {
        "path": "nope.tif", "present": False, "size": None, "sha256": None}


# --------------------------------------------------------- semantic equality

def test_same_mesh_on_identical_content(tmp_path):
    a = make_mesh(tmp_path, "a", mask=b"\x01", meta={"s": 1})
    b = make_mesh(tmp_path, "b", mask=b"\x01", meta={"s": 1})
    assert same_mesh(mesh_manifest(a), mesh_manifest(b)) is True
    assert mesh_manifest(a)["digest"] == mesh_manifest(b)["digest"]


def test_same_mesh_is_none_when_undecidable():
    assert same_mesh({}, {"x": "0" * 64}) is None
    assert same_mesh(None, None) is None


def test_same_mesh_accepts_legacy_plane_dicts():
    legacy = {"x": "a" * 64, "y": "b" * 64, "z": "c" * 64, "mask": None}
    other = {"x.tif": "a" * 64, "y.tif": "b" * 64, "z.tif": "c" * 64}
    assert same_mesh(legacy, other) is True


# ------------------------------------------------------------- the base kind

def test_base_kind_from_manifests(tmp_path):
    orig = mesh_manifest(make_mesh(tmp_path, "orig", meta={"s": 1}))
    same = dict(orig)
    other = mesh_manifest(make_mesh(tmp_path, "rep", planes=b"QQQ",
                                    meta={"s": 1}))
    assert base_kind_from_manifests(same, orig) is BaseKind.ORIGINAL
    assert base_kind_from_manifests(other, orig) is \
        BaseKind.DISPLACEMENT_REPAIRED
    assert base_kind_from_manifests({}, orig) is None


def test_base_kind_ignores_paths(tmp_path):
    """Same bytes at two different paths is the SAME mesh."""
    a = make_mesh(tmp_path / "download", "seg", meta={"s": 1})
    b = make_mesh(tmp_path / "fresh_workdir", "seg", meta={"s": 1})
    assert a != b
    assert base_kind_from_manifests(mesh_manifest(a),
                                    mesh_manifest(b)) is BaseKind.ORIGINAL


def test_legacy_spellings_normalise():
    assert normalise_base_kind("original_published") is BaseKind.ORIGINAL
    assert normalise_base_kind("original") is BaseKind.ORIGINAL
    assert normalise_base_kind("displacement_repaired") is \
        BaseKind.DISPLACEMENT_REPAIRED
    assert normalise_base_kind("nonsense") is None
    assert normalise_base_kind(None) is None


def test_verify_base_kind_agreement(tmp_path):
    orig = mesh_manifest(make_mesh(tmp_path, "orig"))
    v = verify_base_kind("original_published", orig, orig)
    assert v["agrees"] and v["verified"] == "original"


def test_verify_base_kind_disagreement_is_loud(tmp_path):
    orig = mesh_manifest(make_mesh(tmp_path, "orig"))
    rep = mesh_manifest(make_mesh(tmp_path, "rep", planes=b"QQQ"))
    v = verify_base_kind("original", rep, orig, label="seg1")
    assert v["agrees"] is False
    assert v["verified"] == "displacement_repaired"
    assert "seg1" in v["reason"] and "disagrees" in v["reason"]
    with pytest.raises(BaseKindDisagreement):
        require_base_kind("original", rep, orig, label="seg1")


def test_verify_base_kind_undecidable_never_reads_as_agreement(tmp_path):
    orig = mesh_manifest(make_mesh(tmp_path, "orig"))
    v = verify_base_kind("original", {}, orig)
    assert v["agrees"] is False and v["decidable"] is False
    with pytest.raises(BaseKindDisagreement):
        require_base_kind("original", {}, orig)


def test_a_mask_only_difference_makes_it_a_repaired_base(tmp_path):
    """The regression the manifest exists for: identical coordinates, a
    different mask, and the base is NOT the original mesh."""
    orig = make_mesh(tmp_path, "orig", mask=b"\x01\x01")
    cut = make_mesh(tmp_path, "cut", mask=b"\x01\x00")
    assert mesh_manifest(orig)["x"] == mesh_manifest(cut)["x"]
    assert base_kind_from_manifests(mesh_manifest(cut),
                                    mesh_manifest(orig)) is \
        BaseKind.DISPLACEMENT_REPAIRED
