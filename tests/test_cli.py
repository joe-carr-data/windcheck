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

from windcheck import cli, pipeline, tifxyz

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

    # The envelope is not cosmetic. PointCollections::loadFromJSON refuses any
    # file lacking this exact key/value, then reads `collections`; without both
    # it logs "incorrect version or missing version info" and loads nothing.
    # An earlier version of this test asserted our own bare {"1": ...} shape
    # and passed while every emitted overlay was unopenable in VC3D -- it
    # checked that we wrote what we meant to write, not that the viewer reads
    # it. Assert the consumer's contract instead.
    assert doc["vc_pointcollections_json_version"] == "1"
    assert set(doc) == {"vc_pointcollections_json_version", "collections"}
    assert list(doc["collections"]) == ["1"]
    coll = doc["collections"]["1"]
    # Required by from_json(Collection): .at() on any of these throws.
    for required in ("name", "points", "metadata", "color"):
        assert required in coll, f"loadFromJSON requires {required!r}"
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


def test_certificate_discloses_the_upstream_validity_gap(tmp_path):
    """We keep cells the upstream loader drops; the certificate must say so.

    QuadSurface::load rewrites any point with z <= 0 to the invalid sentinel
    before masking. Our reader does not, so our valid set is a superset and a
    crossing confined to those cells is one the pipeline never sees. The rule
    is disclosed rather than adopted, because adopting it would silently
    restate every published count.
    """
    import numpy as np
    import tifffile
    from windcheck import tifxyz

    d = tmp_path / "zfloor.tifxyz"
    d.mkdir()
    x = np.ones((4, 4), np.float32)
    y = np.ones((4, 4), np.float32)
    z = np.ones((4, 4), np.float32)
    z[0, 0] = -3.0     # below the floor: we keep it, the loader does not
    z[1, 1] = 0.0      # exactly zero is discarded upstream too (z <= 0)
    x[3, 3] = y[3, 3] = z[3, 3] = -1.0          # sentinel: invalid for both
    for ax, arr in (("x", x), ("y", y), ("z", z)):
        tifffile.imwrite(d / f"{ax}.tif", arr)

    s = tifxyz.read(d)
    assert s.n_valid == 15                       # only the sentinel is dropped
    assert s.n_valid_pipeline == 13              # ...plus the two z <= 0 cells
    assert s.z_floor_cells == 2
    assert s.valid[0, 0] and not s.valid_pipeline[0, 0]


def test_written_meta_carries_a_recomputed_bbox(tmp_path):
    """A bbox inherited from the source describes the wrong surface.

    Consumers filter on it, so a stale bbox drops a mesh from their inputs
    without an error anywhere. This is the failure tifxyz-repair exists to
    fix upstream; our writers must not reintroduce it.
    """
    import json
    import numpy as np
    from windcheck import tifxyz

    src, dst = tmp_path / "src", tmp_path / "dst"
    src.mkdir(); dst.mkdir()
    (src / "meta.json").write_text(json.dumps(
        {"format": "tifxyz", "scale": [2.0, 3.0],
         "bbox": [[0, 0, 0], [9999, 9999, 9999]], "uuid": "keep-me"}))

    pts = np.zeros((3, 3, 3), np.float32)
    pts[..., 0] = 5.0; pts[..., 1] = 6.0; pts[..., 2] = 7.0
    pts[0, 0] = (1.0, 2.0, 3.0)
    valid = np.ones((3, 3), bool)
    valid[2, 2] = False
    pts[2, 2] = (1e6, 1e6, 1e6)          # invalid: must not reach the bbox

    tifxyz.write_meta(src, dst, pts, valid)
    m = json.loads((dst / "meta.json").read_text())
    assert m["bbox"] == [[1.0, 2.0, 3.0], [5.0, 6.0, 7.0]]
    assert m["scale"] == [2.0, 3.0] and m["uuid"] == "keep-me"


# ------------------------------------------------------------- check-pairs

@needs_engine
def test_check_pairs_separates_conflicting_from_clean(tmp_path):
    """Two surfaces that cross must be flagged; two that do not must not.

    The pair is tested by writing both into one grid separated by invalid
    rows, so the risk is that adjacency exclusion hides a real contact, or
    that a surface's own self-intersection is blamed on the pair. Both are
    checked here: the crossing fixture self-intersects on its own, and it
    still comes back clean when paired with a distant plane.
    """
    import numpy as np
    from windcheck import pairs

    def plane(name, x0, z0):
        v, u = np.meshgrid(np.arange(30), np.arange(40), indexing="ij")
        P = np.stack([u * 4.0 + x0, v * 4.0, np.full_like(u, z0, float)], -1)
        return write_tifxyz(tmp_path / name, P)

    def vertical(name, x0):
        v, u = np.meshgrid(np.arange(30), np.arange(40), indexing="ij")
        # A sheet driven through where `plane` sits. It is deliberately
        # TILTED and offset. An axis-aligned sheet at a whole multiple of
        # the grid pitch meets the plane exactly along its mesh edges, and
        # the engine calls that grazing rather than transverse -- correctly,
        # since an exact edge-on meeting is a touch and not a penetration.
        # Getting a real crossing out of a synthetic fixture means avoiding
        # every such coincidence.
        P = np.stack([x0 + 62.0 + (v - 15) * 0.7, v * 4.0,
                      (u - 20) * 4.0 + 2.0], -1)
        return write_tifxyz(tmp_path / name, P)

    flat = plane("flat.tifxyz", 0.0, 0.0)
    far = plane("far.tifxyz", 0.0, 5000.0)
    thru = vertical("thru.tifxyz", 0.0)
    work = tmp_path / "w"
    work.mkdir()

    apart = pairs.classify(flat, far, work, "apart")
    assert apart.verdict in (pairs.NO_CONFLICT, pairs.NOT_TESTABLE)
    assert apart.transverse_both == 0

    crossed = pairs.classify(flat, thru, work, "crossed")
    assert crossed.verdict == pairs.CONFLICT, crossed.reason
    assert crossed.transverse_both > 0
    # the contacts are BETWEEN the two, not either surface's own
    assert crossed.self_a == 0

    # a clean verdict must never read as a merge-safety guarantee
    assert "not a statement that they are compatible" in \
        pairs.VERDICT_NOTE[pairs.NO_CONFLICT]


def _pair_fixtures(tmp_path):
    """Four surfaces: a plane, the same plane far away, one driven through
    it, and the lemniscate strip that passes through ITSELF. The last is
    what makes the fixtures able to tell a surface's own contacts apart
    from the pair's -- without it every self count is trivially zero."""
    import numpy as np

    def plane(name, z0):
        v, u = np.meshgrid(np.arange(30), np.arange(40), indexing="ij")
        P = np.stack([u * 4.0, v * 4.0, np.full_like(u, z0, float)], -1)
        return write_tifxyz(tmp_path / name, P)

    v, u = np.meshgrid(np.arange(30), np.arange(40), indexing="ij")
    # Tilted and offset: see the note in the fixture above. An axis-aligned
    # sheet at a whole multiple of the grid pitch only grazes.
    thru = write_tifxyz(tmp_path / "thru.tifxyz", np.stack(
        [62.0 + (v - 15) * 0.7, v * 4.0, (u - 20) * 4.0 + 2.0], -1))
    return (plane("flat.tifxyz", 0.0), plane("far.tifxyz", 5000.0), thru,
            crossing_mesh(tmp_path / "lem"))


@needs_engine
def test_check_pairs_batch_matches_stitch(tmp_path):
    """The batched census and the stitched one must agree on every field.

    These are two genuinely different ways of putting a pair in front of
    the engine. `classify` lays both surfaces out in one grid separated by
    invalid rows and never tells the engine there is more than one surface;
    `classify_many` loads them as distinct surfaces from one atlas and
    reads the pair off the contact's surface ids. If the surface-id
    bookkeeping, the adjacency scoping or the batch's shared broad-phase
    grid were wrong, the two would part company here.
    """
    from windcheck import pairs

    flat, far, thru, lem = _pair_fixtures(tmp_path)
    work = tmp_path / "w"
    work.mkdir()
    edges = [(flat, far), (flat, thru), (far, thru), (flat, lem)]

    batched = pairs.classify_many(edges, work, tag="t")
    single = [pairs.classify(a, b, work, f"s{i}")
              for i, (a, b) in enumerate(edges)]

    for got, want, (a, b) in zip(batched, single, edges):
        assert got.as_dict() == want.as_dict(), f"{a.name} x {b.name}"

    # the two planes are 5,000 voxels apart, so there is nothing to decide
    assert batched[0].verdict == pairs.NOT_TESTABLE
    assert batched[2].verdict == pairs.NOT_TESTABLE
    # the tilted sheet really is driven through the plane
    assert batched[1].verdict == pairs.CONFLICT
    assert batched[1].transverse_both > 0 and batched[1].self_a == 0
    # and the lemniscate's own crossing is charged to the lemniscate
    assert batched[3].self_a == 0 and batched[3].self_b > 0


@needs_engine
def test_batch_census_recovers_each_surfaces_own_verdict(tmp_path):
    """A surface's self-census must not depend on what shares its atlas.

    This is what makes batching legitimate. The engine's broad phase owns
    each triangle pair at the cell holding the minimum corner of their
    overlapping boxes, so a pair is tested exactly once however the grid's
    origin falls -- and the grid's origin does move, because a batch's
    bounding box spans every surface in it.
    """
    from windcheck import pairs, pipeline

    meshes = _pair_fixtures(tmp_path)
    work = tmp_path / "w"
    work.mkdir()

    stats = pairs.census_batch(list(meshes), work, "own")
    for i, mesh in enumerate(meshes):
        alone = sum(pipeline.run_engine(mesh, f"a{i}", work, d)[1]["transverse"]
                    for d in (0, 1))
        assert stats.get((i, i), {}).get("total", 0) == alone, mesh.name
    assert stats[(3, 3)]["total"] > 0     # the lemniscate is a real positive


@needs_engine
def test_check_pairs_verdicts_do_not_depend_on_batching(tmp_path):
    """How the edge list is cut into batches must not change any answer.

    Batching is a performance decision -- it bounds how much geometry the
    engine holds at once -- and a performance decision that moved a verdict
    would be a bug of the worst kind, since it would only show up on inputs
    large enough to split. Forcing one edge per batch takes the shared
    broad-phase grid away entirely.
    """
    from windcheck import pairs

    flat, far, thru, lem = _pair_fixtures(tmp_path)
    work = tmp_path / "w"
    work.mkdir()
    edges = [(flat, far), (flat, thru), (far, thru), (flat, lem), (thru, lem)]

    whole = pairs.classify_many(edges, work, tag="whole")
    split = pairs.classify_many(edges, work, tag="split", max_batch_surfaces=2)
    tiny = pairs.classify_many(edges, work, tag="tiny",
                               max_batch_triangles=1)
    for a, b, c in zip(whole, split, tiny):
        assert a.as_dict() == b.as_dict() == c.as_dict()


@needs_engine
def test_check_pairs_sees_through_an_aliased_surface(tmp_path):
    """One surface named twice is one surface, not a pair.

    Loaded as two atlas entries it would be its own perfect duplicate:
    every contact between the copies is coplanar, never transverse, so the
    edge would come back `no_transverse_conflict` -- a clean verdict for a
    pair that was never tested. It has to be refused instead.
    """
    from windcheck import pairs

    flat, _far, _thru, lem = _pair_fixtures(tmp_path)
    alias = tmp_path / "lem_alias.tifxyz"
    alias.symlink_to(lem, target_is_directory=True)
    work = tmp_path / "w"
    work.mkdir()

    same, other = pairs.classify_many([(lem, alias), (flat, lem)], work,
                                      tag="al")
    assert same.verdict == pairs.NOT_TESTABLE
    assert "same surface" in same.reason
    assert same.transverse_both == 0
    # and the aliased edge did not disturb the real one sharing its batch
    assert other.as_dict() == pairs.classify(flat, lem, work, "ref").as_dict()


def _grid4(tmp_path, name, P, valid=None):
    import numpy as np
    P = np.asarray(P, float).copy()
    if valid is not None:
        P[~valid] = -1.0
    return write_tifxyz(tmp_path / name, P)


@needs_engine
def test_a_found_contact_outranks_the_overlap_gate(tmp_path):
    """Two planes crossing through each other's interiors must be flagged.

    `min_overlap_points` counts vertices inside the two surfaces' shared
    bounding box. Orthogonal planes cross with that box degenerate in two
    axes and no vertex of either inside it, so the count is zero -- while
    the engine finds twelve transverse contacts. Applying the gate before
    the census threw those away and returned `not_testable`, which is the
    one failure this tool must not have: a real conflict reported as
    something nobody needs to look at.
    """
    import numpy as np
    from windcheck import pairs

    g = np.array([0.0, 10.0, 20.0, 30.0])
    X, Y = np.meshgrid(g, g, indexing="ij")
    flat = _grid4(tmp_path, "h.tifxyz", np.stack([X, Y, np.zeros_like(X)], -1))
    # crosses flat at x = 15, z = 0 -- between grid lines of both, so
    # neither contributes a vertex to the shared box
    Yv, Zv = np.meshgrid(g, np.array([-15.0, -5.0, 5.0, 15.0]), indexing="ij")
    vert = _grid4(tmp_path, "v.tifxyz",
                  np.stack([np.full_like(Yv, 15.0), Yv, Zv], -1))
    work = tmp_path / "w"
    work.mkdir()

    batched = pairs.classify_many([(flat, vert)], work, tag="g")[0]
    stitched = pairs.classify(flat, vert, work, "gs")
    assert batched.overlap_points == 0        # the gate would have refused
    assert batched.verdict == pairs.CONFLICT
    assert batched.transverse_both > 0
    assert stitched.as_dict() == batched.as_dict()


@needs_engine
def test_check_pairs_refuses_a_surface_with_no_triangles(tmp_path):
    """A surface can be all-valid and still contribute nothing to a census.

    A checkerboard of valid cells has no quad with four valid corners, so
    it builds no triangle. Reporting such an edge as free of contact would
    be a clean verdict on a surface that was never tested; and when every
    surface in a batch is like that, the engine used to exit without its
    summary and the caller raised IndexError on empty stdout.
    """
    import numpy as np
    from windcheck import pairs

    g = np.arange(4) * 10.0
    X, Y = np.meshgrid(g, g, indexing="ij")
    P = np.stack([X, Y, np.zeros_like(X)], -1)
    board = (np.arange(4)[:, None] + np.arange(4)[None, :]) % 2 == 0
    empty = _grid4(tmp_path, "e1.tifxyz", P, board)
    empty2 = _grid4(tmp_path, "e2.tifxyz", P + np.array([1.0, 1.0, 0.0]), board)
    solid = _grid4(tmp_path, "s.tifxyz", P)
    work = tmp_path / "w"
    work.mkdir()

    assert pairs.n_triangles(tifxyz.read(empty)) == 0
    assert pairs.n_triangles(tifxyz.read(solid)) > 0

    both, one = pairs.classify_many([(empty, empty2), (empty, solid)], work,
                                    tag="e")
    for r in (both, one, pairs.classify(empty, solid, work, "es")):
        assert r.verdict == pairs.NOT_TESTABLE
        assert "no triangle" in r.reason


@needs_engine
def test_engine_reports_an_empty_census_rather_than_saying_nothing(tmp_path):
    """Zero triangles is a result, not a failure to produce one."""
    import numpy as np
    from windcheck import pipeline

    g = np.arange(4) * 10.0
    X, Y = np.meshgrid(g, g, indexing="ij")
    board = (np.arange(4)[:, None] + np.arange(4)[None, :]) % 2 == 0
    empty = _grid4(tmp_path, "z.tifxyz",
                   np.stack([X, Y, np.zeros_like(X)], -1), board)
    csvp, counts = pipeline.run_engine(empty, "z", tmp_path / "w", 0)
    assert counts["triangles"] == 0 and counts["transverse"] == 0
    assert csvp.read_text().strip() == pipeline.SCHEMA_V2_HEADER


@needs_engine
def test_check_pairs_refuses_duplicate_geometry(tmp_path):
    """One surface under two names is not a pair, however it is spelled.

    Path resolution catches a symlink; two separately written but identical
    directories it cannot. Censused as two atlas entries the copies agree
    everywhere, so every contact between them is coplanar and never
    counted, and the edge would come back clean without having been tested.
    """
    import numpy as np
    from windcheck import pairs

    g = np.arange(4) * 10.0
    X, Y = np.meshgrid(g, g, indexing="ij")
    P = np.stack([X, Y, np.zeros_like(X)], -1)
    one = _grid4(tmp_path, "one.tifxyz", P)
    copy = _grid4(tmp_path, "copy.tifxyz", P)   # same bytes, different path
    work = tmp_path / "w"
    work.mkdir()

    for r in (pairs.classify_many([(one, copy)], work, tag="d")[0],
              pairs.classify(one, copy, work, "ds")):
        assert r.verdict == pairs.NOT_TESTABLE
        assert "same geometry" in r.reason


def test_check_pairs_is_report_only():
    """Like `check`, there is no flag on `check-pairs` that writes geometry."""
    actions = {a.dest for a in cli.build_parser()._subparsers._group_actions[0]
               .choices["check-pairs"]._actions}
    for forbidden in ("transform", "repair", "excise", "fix", "write_mesh",
                      "in_place", "merge"):
        assert forbidden not in actions


def test_transaction_refuses_unknown_sidecars(tmp_path):
    """The transaction never guesses sidecar semantics (exit 2)."""
    d = tmp_path / "seg"
    d.mkdir()
    for n in ("x.tif", "y.tif", "z.tif"):
        (d / n).write_bytes(b"")
    (d / "meta.json").write_text("{}")
    (d / "mystery.bin").write_bytes(b"x")
    rc = cli.main(["transaction", str(d), "--out", str(tmp_path / "out"),
                   "--report", str(tmp_path / "r.json")])
    assert rc == 2
    import json
    assert "mystery.bin" in json.loads((tmp_path / "r.json").read_text())["note"]


def test_transaction_refuses_existing_out(tmp_path):
    d = tmp_path / "seg"
    d.mkdir()
    for n in ("x.tif", "y.tif", "z.tif"):
        (d / n).write_bytes(b"")
    out = tmp_path / "out"
    out.mkdir()
    rc = cli.main(["transaction", str(d), "--out", str(out)])
    assert rc == 2


def test_transaction_r45_contract(tmp_path):
    """R45 contract: meta.json required (2); report inside input refused
    (2); nonexistent out parent works; certificate preserved in output."""
    import json as _json
    import numpy as np
    import tifffile
    d = tmp_path / "seg"
    d.mkdir()
    H, W = 80, 80
    x, yv = np.meshgrid(np.arange(W, dtype=np.float32),
                        np.arange(H, dtype=np.float32))
    z = np.full((H, W), 5.0, np.float32)
    tifffile.imwrite(d / "x.tif", x)
    tifffile.imwrite(d / "y.tif", yv)
    tifffile.imwrite(d / "z.tif", z)
    # missing meta.json -> 2
    assert cli.main(["transaction", str(d),
                     "--out", str(tmp_path / "o1")]) == 2
    (d / "meta.json").write_text(_json.dumps(
        {"scale": [1.0, 1.0], "uuid": "t", "type": "seg",
         "format": "tifxyz"}))
    # report inside input -> 2
    assert cli.main(["transaction", str(d), "--out", str(tmp_path / "o2"),
                     "--report", str(d / "r.json")]) == 2
    # nonexistent out parent + certificate preserved
    out = tmp_path / "deep" / "nested" / "final"
    rc = cli.main(["transaction", str(d), "--out", str(out)])
    assert rc == 0
    assert (out / "windcheck_transaction" / "certificate.json").exists()
    # R46: the authoritative report is committed inside the output
    rep = _json.loads((out / "windcheck_transaction" /
                       "transaction_report.json").read_text())
    assert rep["committed"] is True and rep["exit_code"] == 0


def test_transaction_ignores_sibling_facemap(tmp_path):
    """A facemap next to (not inside) the input is never consumed."""
    import json as _json
    import numpy as np
    import tifffile
    d = tmp_path / "seg"
    d.mkdir()
    H, W = 80, 80
    x, yv = np.meshgrid(np.arange(W, dtype=np.float32),
                        np.arange(H, dtype=np.float32))
    tifffile.imwrite(d / "x.tif", x)
    tifffile.imwrite(d / "y.tif", yv)
    tifffile.imwrite(d / "z.tif", np.full((H, W), 5.0, np.float32))
    (d / "meta.json").write_text(_json.dumps(
        {"scale": [1.0, 1.0], "uuid": "t", "type": "seg",
         "format": "tifxyz"}))
    (tmp_path / "unrelated_facemap.i32").write_bytes(b"\x00" * 16)
    out = tmp_path / "out"
    rc = cli.main(["transaction", str(d), "--out", str(out),
                   "--adapter", "scrollfiesta",
                   "--report", str(tmp_path / "r.json")])
    assert rc == 0
    rep = _json.loads((tmp_path / "r.json").read_text())
    assert "unrelated_facemap.i32" not in rep["output_files_sha256"]
    assert not (out / "unrelated_facemap.i32").exists()


def _clean_seg(tmp_path, name="seg"):
    """A minimal transverse-clean tifxyz input (exit-0 path)."""
    import json as _json
    import numpy as np
    import tifffile
    d = tmp_path / name
    d.mkdir()
    H, W = 80, 80
    x, yv = np.meshgrid(np.arange(W, dtype=np.float32),
                        np.arange(H, dtype=np.float32))
    tifffile.imwrite(d / "x.tif", x)
    tifffile.imwrite(d / "y.tif", yv)
    tifffile.imwrite(d / "z.tif", np.full((H, W), 5.0, np.float32))
    (d / "meta.json").write_text(_json.dumps(
        {"scale": [1.0, 1.0], "uuid": "t", "type": "seg",
         "format": "tifxyz"}))
    return d


def test_transaction_report_copy_failure_never_reverses_commit(tmp_path):
    """R46: --report is a post-commit copy; its failure warns but the
    committed exit code stands and the promoted output survives."""
    d = _clean_seg(tmp_path)
    blocker = tmp_path / "not_a_dir"
    blocker.write_bytes(b"")          # report parent is a FILE
    out = tmp_path / "out"
    rc = cli.main(["transaction", str(d), "--out", str(out),
                   "--report", str(blocker / "r.json")])
    assert rc == 0
    assert out.is_dir()
    assert (out / "windcheck_transaction" /
            "transaction_report.json").exists()


def test_transaction_resolves_validator_through_path(tmp_path, monkeypatch):
    """R46: a bare validator name resolves via PATH once; the resolved
    binary is what gets invoked and hashed."""
    import hashlib
    import json as _json
    d = _clean_seg(tmp_path)
    bindir = tmp_path / "bin"
    bindir.mkdir()
    fake = bindir / "fakecross"
    fake.write_text(
        "#!/bin/sh\n"
        "# args: <surface> -o <report>\n"
        'printf \'{"clean_of_transverse_self_intersection": true,'
        ' "census": []}\' > "$3"\n')
    fake.chmod(0o755)
    monkeypatch.setenv("PATH",
                       f"{bindir}:{__import__('os').environ['PATH']}")
    out = tmp_path / "out"
    rc = cli.main(["transaction", str(d), "--out", str(out),
                   "--official-validator", "fakecross",
                   "--report", str(tmp_path / "r.json")])
    assert rc == 0
    rep = _json.loads((tmp_path / "r.json").read_text())
    assert rep["official_validator_path"] == str(fake.resolve())
    want = hashlib.sha256(fake.read_bytes()).hexdigest()
    assert rep["official_validator"]["binary_sha256"] == want


def test_transaction_refuses_unresolvable_validator(tmp_path):
    d = _clean_seg(tmp_path)
    rc = cli.main(["transaction", str(d), "--out", str(tmp_path / "out"),
                   "--official-validator", "no_such_validator_xyz"])
    assert rc == 2
    assert not (tmp_path / "out").exists()


def test_transaction_adapter_exit_disagreement_is_internal(tmp_path,
                                                           monkeypatch):
    """R46: a transform report claiming handover cannot outrank a nonzero
    adapter exit; the disagreement is exit 1 and nothing is promoted."""
    import json as _json
    import types
    from windcheck import transaction as tx
    d = _clean_seg(tmp_path)

    def fake_run(cmd, **kw):
        rp = Path(cmd[cmd.index("--report") + 1])
        rp.parent.mkdir(parents=True, exist_ok=True)
        cert = rp.parent / "cert.json"
        cert.write_text("{}")
        rp.write_text(_json.dumps(
            {"handed_over": True, "terminal_disposition": "already_clean",
             "certificate": str(cert)}))
        return types.SimpleNamespace(returncode=1, stdout="", stderr="boom")

    monkeypatch.setattr(tx.subprocess, "run", fake_run)
    rc = cli.main(["transaction", str(d), "--out", str(tmp_path / "out"),
                   "--report", str(tmp_path / "r.json")])
    assert rc == 1
    assert not (tmp_path / "out").exists()
    rep = _json.loads((tmp_path / "r.json").read_text())
    assert "disagreement" in rep["note"]
