"""What the user-facing command line promises, pinned against the real engine.

These run the actual C++ census on tiny synthetic surfaces rather than mocking
it, because the promises being tested are about what lands on disk after a real
run: that `check` is report-only, that its two artifacts exist and parse, and
that a mesh appears only when `transform` was typed.

The fixtures are a flat plane (clean) and a strip swept along a self-crossing
plane curve (one real transverse crossing, far apart in the grid), both small
enough that the whole file runs in seconds.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
import tifffile

from windcheck import cli, pipeline

ENGINE = pipeline.ENGINE
needs_engine = pytest.mark.skipif(
    not ENGINE.exists(), reason="engines/selfcross is not built")


# ----------------------------------------------------------------- fixtures

def write_tifxyz(path: Path, P: np.ndarray) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    for i, ax in enumerate("xyz"):
        tifffile.imwrite(path / f"{ax}.tif", P[..., i].astype(np.float32))
    (path / "meta.json").write_text(json.dumps(
        {"scale": [1.0, 1.0], "bbox": [[0, 0, 0], [1, 1, 1]]}))
    return path


def flat_mesh(root: Path, nv: int = 60, nu: int = 100) -> Path:
    """A plane. Cannot cross itself, so the verdict must be clean."""
    v, u = np.meshgrid(np.arange(nv), np.arange(nu), indexing="ij")
    P = np.stack([u * 4.0, v * 4.0, np.zeros_like(u, float)], -1)
    return write_tifxyz(root / "segFLAT" / "mesh" / "segFLAT.tifxyz", P)


def crossing_mesh(root: Path, nv: int = 40, nu: int = 160) -> Path:
    """A strip swept along a lemniscate: the sheet passes through itself once.

    The curve `(cos t, sin t cos t)` returns to the origin at two well
    separated parameters, so the two branches meet at a large angle and are far
    apart in the grid -- a genuine non-adjacent transverse crossing, not a
    near-touch that adjacency exclusion would have to arbitrate.
    """
    t = np.linspace(0.0, 2.0 * np.pi, nu)
    cx, cz = 120.0 * np.cos(t), 120.0 * np.sin(t) * np.cos(t)
    y = np.arange(nv) * 5.0
    P = np.empty((nv, nu, 3), np.float64)
    P[..., 0] = cx[None, :]
    P[..., 1] = y[:, None]
    P[..., 2] = cz[None, :]
    return write_tifxyz(root / "segCROSS" / "mesh" / "segCROSS.tifxyz", P)


def tree_hash(root: Path) -> dict:
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


# ------------------------------------------------------------ check: shape

@needs_engine
def test_check_reports_a_real_crossing(tmp_path, capsys):
    mesh = crossing_mesh(tmp_path)
    out = tmp_path / "report"
    assert cli.main(["check", str(mesh), "--out", str(out)]) == 0
    text = capsys.readouterr().out
    assert "NOT clean" in text
    assert "transverse contacts" in text
    assert "Report only" in text


@needs_engine
def test_check_calls_a_plane_clean(tmp_path, capsys):
    mesh = flat_mesh(tmp_path)
    out = tmp_path / "report"
    assert cli.main(["check", str(mesh), "--out", str(out)]) == 0
    assert "clean" in capsys.readouterr().out
    cert = json.loads(next(out.glob("*_check_certificate.json")).read_text())
    assert cert["clean"] is True
    assert cert["measurements"]["transverse_both_diagonals"] == 0


# ------------------------------------------------------- check: report-only

@needs_engine
def test_check_does_not_modify_the_input(tmp_path):
    mesh = crossing_mesh(tmp_path)
    seg_root = mesh.parent.parent
    before = tree_hash(seg_root)
    cli.main(["check", str(mesh), "--out", str(tmp_path / "report")])
    assert tree_hash(seg_root) == before


@needs_engine
def test_check_writes_no_mesh_anywhere(tmp_path):
    mesh = crossing_mesh(tmp_path)
    out = tmp_path / "report"
    cli.main(["check", str(mesh), "--out", str(out)])
    assert [p for p in out.rglob("*.tifxyz")] == []
    assert [p for p in out.rglob("*.tif")] == []


def test_check_has_no_flag_that_transforms():
    """Transformation is unreachable from `check`: it has no option for it."""
    actions = {a.dest for a in cli.build_parser()._subparsers._group_actions[0]
               .choices["check"]._actions}
    for forbidden in ("displaced", "no_displacement", "transform", "repair",
                      "excise", "fix", "write_mesh", "in_place"):
        assert forbidden not in actions


def test_transform_must_be_named_and_needs_an_out_dir(capsys):
    """`transform` is never the default, and never guesses where to write."""
    with pytest.raises(SystemExit):
        cli.main([])                       # no subcommand at all
    with pytest.raises(SystemExit):
        cli.main(["transform", "some/path"])   # --out is required
    assert "--out" in capsys.readouterr().err


# ---------------------------------------------------------- check: artifacts

@needs_engine
def test_certificate_parses_and_carries_the_census(tmp_path):
    mesh = crossing_mesh(tmp_path)
    out = tmp_path / "report"
    cli.main(["check", str(mesh), "--out", str(out)])
    certs = sorted(out.glob("*_check_certificate.json"))
    assert len(certs) == 1
    cert = json.loads(certs[0].read_text())

    assert cert["tool"] == "windcheck check"
    assert cert["report_only"] is True
    assert cert["clean"] is False
    m = cert["measurements"]
    # Both canonical triangulations are censused and reported separately.
    assert m["transverse_d0"] > 0 and m["transverse_d1"] > 0
    assert m["transverse_both_diagonals"] == m["transverse_d0"] + m["transverse_d1"]
    assert m["crossing_events"] >= 1
    assert cert["census"]["parameters"]["diagonals"] == [0, 1]
    assert cert["census"]["parameters"]["exclude"] == 1
    assert set(cert["census"]["csv"]) == {"d0", "d1"}
    assert cert["mesh"]["grid_shape"] == [40, 160]
    assert cert["wall_seconds"] >= 0.0
    assert cert["caveats"]


@needs_engine
def test_point_collection_is_loadable_and_well_formed(tmp_path):
    mesh = crossing_mesh(tmp_path)
    out = tmp_path / "report"
    cli.main(["check", str(mesh), "--out", str(out)])
    doc = json.loads(next(out.glob("*_points.json")).read_text())

    assert list(doc) == ["1"]
    coll = doc["1"]
    assert coll["id"] == 1
    assert coll["name"].startswith("windcheck-crossings-")
    assert len(coll["color"]) == 3
    assert coll["points"], "a crossing was found, so there must be a site"
    for key, pt in coll["points"].items():
        assert str(pt["id"]) == key
        assert len(pt["p"]) == 3
        assert all(isinstance(c, float) for c in pt["p"])
        assert pt["collection_id"] == 1
    # The sites are volume coordinates, so they must sit inside the mesh bounds.
    P = np.stack([tifffile.imread(mesh / f"{ax}.tif") for ax in "xyz"], -1)
    lo, hi = P.reshape(-1, 3).min(0), P.reshape(-1, 3).max(0)
    for pt in coll["points"].values():
        assert np.all(np.asarray(pt["p"]) >= lo - 1e-6)
        assert np.all(np.asarray(pt["p"]) <= hi + 1e-6)


# ------------------------------------------------------------- transform

@needs_engine
def test_transform_writes_a_mesh_and_a_certificate(tmp_path):
    mesh = crossing_mesh(tmp_path)
    out = tmp_path / "fixed"
    rc = cli.main(["transform", str(mesh), "--out", str(out)])

    meshes = sorted(out.glob("*_transformed.tifxyz"))
    assert len(meshes) == 1
    for ax in "xyz":
        assert (meshes[0] / f"{ax}.tif").exists()
    assert (meshes[0] / "mask.tif").exists(), "hybrid invalidation needs both"

    cert = json.loads(next(out.glob("*_transform_certificate.json")).read_text())
    assert cert["tool"] == "windcheck transform"
    assert cert["policy_version"]
    assert cert["policy_hash"]
    assert cert["census_before"]["d0"]["transverse"] > 0
    assert cert["output_mesh"] == str(meshes[0])
    # The clean claim, if made, is made by recensusing the emitted arrays.
    assert cert["status"] in ("clean", "residual_contacts")
    if cert["status"] == "clean":
        assert rc == 0
        assert cert["claimed_clean"] is True
        assert cert["census_after"]["d0"]["transverse"] == 0
        assert cert["census_after"]["d1"]["transverse"] == 0
        assert 0.0 < cert["excised_area_fraction"] < 0.5


@needs_engine
def test_transform_keeps_retained_coordinates_bit_identical(tmp_path):
    mesh = crossing_mesh(tmp_path)
    out = tmp_path / "fixed"
    cli.main(["transform", str(mesh), "--out", str(out)])
    dst = next(out.glob("*_transformed.tifxyz"))
    keep = np.asarray(tifffile.imread(dst / "mask.tif")).astype(bool)
    for ax in "xyz":
        a = np.asarray(tifffile.imread(mesh / f"{ax}.tif"))
        b = np.asarray(tifffile.imread(dst / f"{ax}.tif"))
        assert np.array_equal(a[keep], b[keep])
        assert np.all(b[~keep] == -1.0)


@needs_engine
def test_transform_of_a_clean_mesh_still_emits_one_artifact(tmp_path):
    mesh = flat_mesh(tmp_path)
    out = tmp_path / "fixed"
    assert cli.main(["transform", str(mesh), "--out", str(out)]) == 0
    cert = json.loads(next(out.glob("*_transform_certificate.json")).read_text())
    assert cert["status"] == "already_clean"
    assert cert["n_removed_quads"] == 0
    assert Path(cert["output_mesh"]).is_dir()


@needs_engine
def test_transform_does_not_modify_the_input(tmp_path):
    mesh = crossing_mesh(tmp_path)
    seg_root = mesh.parent.parent
    before = tree_hash(seg_root)
    cli.main(["transform", str(mesh), "--out", str(tmp_path / "fixed")])
    assert tree_hash(seg_root) == before


@needs_engine
def test_transform_prefers_a_named_displacement_base(tmp_path):
    """`--displaced` is what makes the two-stage play explicit and auditable."""
    mesh = crossing_mesh(tmp_path)
    moved = crossing_mesh(tmp_path / "alt")
    out = tmp_path / "fixed"
    cli.main(["transform", str(mesh), "--out", str(out),
              "--displaced", str(moved)])
    cert = json.loads(next(out.glob("*_transform_certificate.json")).read_text())
    assert cert["displacement"]["applied"] is True
    assert cert["displacement"]["base_mesh"] == str(moved)
    assert cert["input_mesh"] == str(moved)
    assert cert["original_mesh"] == str(mesh)


# --------------------------------------------------------------- discovery

def test_find_mesh_accepts_a_segment_directory(tmp_path):
    mesh = flat_mesh(tmp_path)
    seg = mesh.parent.parent
    assert pipeline.find_mesh(seg) == mesh
    assert pipeline.find_mesh(mesh) == mesh
    assert pipeline.find_mesh(tmp_path / "nothing") is None
    assert pipeline.segment_name(mesh) == "segFLAT"
