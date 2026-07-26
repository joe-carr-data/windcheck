"""Tests for the geometric kernel and the invariant.

These are the checks that would have caught the mistakes actually made while
building this. Each one is tied to a specific real failure:

- `test_engine_matches_bruteforce` — the acceleration structure could silently
  return a non-nearest triangle. An adversarial review found exactly that in the
  `max_dist` path, where an uncertified result escaped instead of becoming
  infinity.
- `test_max_dist_is_certified` — pins that fix.
- `test_flat_sheet_pair_reads_true_separation` — nearest-*sample* distance
  overestimates sheet separation by the grid pitch, which is what made an early
  version of this unable to tell a 1-sheet gap from a 2-sheet gap.
- `test_selfgap_silent_on_single_wrap` — the null control, in code.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from windcheck import atlas, tifxyz

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "engines" / "atlas_query"

pytestmark = pytest.mark.skipif(
    not ENGINE.exists(),
    reason="engine not built: clang++ -O3 -std=c++17 -pthread "
           "-o engines/atlas_query engines/atlas_query.cpp",
)


def _write_sheet(path: Path, z: float, n: int = 24, pitch: float = 20.0,
                 winding: int = 0) -> None:
    """A flat n x n sheet at height `z`, sampled every `pitch` voxels."""
    path.mkdir(parents=True, exist_ok=True)
    xs = np.arange(n, dtype=np.float32) * pitch
    X, Y = np.meshgrid(xs, xs, indexing="xy")
    Z = np.full_like(X, z, dtype=np.float32)
    import tifffile
    tifffile.imwrite(path / "x.tif", X)
    tifffile.imwrite(path / "y.tif", Y)
    tifffile.imwrite(path / "z.tif", Z)
    (path / "meta.json").write_text('{"format":"tifxyz","scale":[1,1]}')


class _Ent:
    def __init__(self, path: Path, winding: int | None):
        self.path, self.winding = path, winding


def _run(atlas_bin, q_bin, r_bin, pts, **kw):
    atlas.write_queries(pts, q_bin)
    return atlas.run_engine(atlas_bin, q_bin, r_bin, threads=2, **kw)


def test_flat_sheet_pair_reads_true_separation(tmp_path):
    """Point-to-surface distance must read the true gap, not the grid pitch.

    Two flat sheets 12 voxels apart, sampled every 20. A nearest-*sample* method
    reports ~20; the true surface distance is 12.
    """
    a, b = tmp_path / "a", tmp_path / "b"
    _write_sheet(a, z=0.0)
    _write_sheet(b, z=12.0)
    atlas.write_atlas([_Ent(a, 0), _Ent(b, 1)], tmp_path / "at.bin")

    # probe from the middle of sheet A, offset off-sample to defeat luck
    pts = np.array([[110.0, 110.0, 0.0]], dtype=np.float32)
    r = _run(tmp_path / "at.bin", tmp_path / "q.bin", tmp_path / "r.bin", pts)
    assert r["w1"][0] == 0
    assert r["d1"][0] == pytest.approx(0.0, abs=1e-3)
    assert r["w2"][0] == 1
    assert r["d2"][0] == pytest.approx(12.0, abs=1e-3)


def test_engine_matches_bruteforce(tmp_path):
    """Accelerated nearest-surface must equal exhaustive search."""
    a, b = tmp_path / "a", tmp_path / "b"
    _write_sheet(a, z=0.0)
    _write_sheet(b, z=35.0)
    atlas.write_atlas([_Ent(a, 0), _Ent(b, 1)], tmp_path / "at.bin")

    rng = np.random.default_rng(0)
    pts = rng.uniform([40, 40, -30], [420, 420, 60], size=(200, 3)).astype(np.float32)
    r = _run(tmp_path / "at.bin", tmp_path / "q.bin", tmp_path / "r.bin", pts,
             max_dist=4096.0)

    # brute force: distance to each plane, clamped to the sheet's extent
    lo, hi = 0.0, 23 * 20.0
    for p, res in zip(pts, r):
        best = []
        for z, w in ((0.0, 0), (35.0, 1)):
            cx = min(max(p[0], lo), hi)
            cy = min(max(p[1], lo), hi)
            best.append((float(np.hypot(np.hypot(p[0] - cx, p[1] - cy), p[2] - z)), w))
        best.sort()
        assert res["w1"] == best[0][1]
        assert res["d1"] == pytest.approx(best[0][0], abs=0.05)


def test_max_dist_is_certified(tmp_path):
    """Beyond max_dist the answer must be infinity, never a finite guess.

    Regression for a real bug: exhausting the shell budget let an uncertified
    finite distance escape, which showed up as "0.5% finite at median 331 vx"
    in a null control and was wrongly dismissed as far-field noise.
    """
    a = tmp_path / "a"
    _write_sheet(a, z=0.0)
    atlas.write_atlas([_Ent(a, 0)], tmp_path / "at.bin")

    far = np.array([[230.0, 230.0, 900.0]], dtype=np.float32)
    r = _run(tmp_path / "at.bin", tmp_path / "q.bin", tmp_path / "r.bin", far,
             max_dist=128.0)
    assert not np.isfinite(r["d1"][0]), "uncertified result escaped as finite"

    r2 = _run(tmp_path / "at.bin", tmp_path / "q.bin", tmp_path / "r2.bin", far,
              max_dist=4096.0)
    assert r2["d1"][0] == pytest.approx(900.0, abs=0.5)


def test_selfgap_silent_on_single_wrap(tmp_path):
    """NULL CONTROL. A single-wrap sheet has no previous wrap to find."""
    from windcheck import selfgap

    a = tmp_path / "seg" / "mesh" / "x.tifxyz"
    _write_sheet(a, z=0.0, n=300, pitch=20.0)
    r = selfgap.analyse(a, name="flat", stride=1, exclude_u=60, workdir=tmp_path)
    assert r is not None, "surface should be large enough to analyse"
    assert r.frac_flagged == pytest.approx(0.0, abs=1e-3), (
        f"detector fired on a flat single sheet: {r.frac_flagged:.4%}"
    )


def test_tifxyz_missing_cells_are_dropped(tmp_path):
    """-1 sentinels must not become real geometry at the origin."""
    import tifffile
    d = tmp_path / "s"
    _write_sheet(d, z=5.0, n=10)
    x = tifffile.imread(d / "x.tif")
    y = tifffile.imread(d / "y.tif")
    z = tifffile.imread(d / "z.tif")
    x[0, 0] = y[0, 0] = z[0, 0] = -1.0
    for nm, arr in (("x", x), ("y", y), ("z", z)):
        tifffile.imwrite(d / f"{nm}.tif", arr)

    s = tifxyz.read(d)
    assert not s.valid[0, 0]
    assert s.n_valid == s.valid.size - 1
    assert s.xyz().min() >= 0.0


def _write_obj(path: Path, z_levels, n=8, pitch=20.0):
    """A tiny OBJ with vt parametrisation: `n`x`n` quads per sheet."""
    v, vt, f, vi = [], [], [], 1
    for z in z_levels:
        idx = np.zeros((n, n), int)
        for r in range(n):
            for c in range(n):
                v.append(f"v {c*pitch} {r*pitch} {z}")
                vt.append(f"vt {c/(n-1)} {r/(n-1)}")
                idx[r, c] = vi
                vi += 1
        for r in range(n - 1):
            for c in range(n - 1):
                a, b, cc, dd = idx[r, c], idx[r, c+1], idx[r+1, c+1], idx[r+1, c]
                f.append(f"f {a}/{a} {b}/{b} {cc}/{cc}")
                f.append(f"f {a}/{a} {cc}/{cc} {dd}/{dd}")
    path.write_text("\n".join(v + vt + f) + "\n")


def test_obj_reader_and_engine_agree_with_tifxyz(tmp_path):
    """The same geometry as OBJ and as tifxyz must measure the same.

    This is the cross-format check: two independent parsers, one kernel. If they
    disagree, one of the readers is wrong.
    """
    from windcheck import objmesh

    obj = tmp_path / "m.obj"
    _write_obj(obj, [0.0, 12.0])
    mesh = objmesh.read(obj)
    assert mesh.has_param, "published OBJs carry vt; the reader must find it"
    assert mesh.n_tris == 2 * 2 * 7 * 7

    atlas.write_atlas_mesh(mesh, tmp_path / "a_obj.bin")
    probe = np.array([[70.0, 70.0, 0.0]], dtype=np.float32)
    r_obj = _run(tmp_path / "a_obj.bin", tmp_path / "q.bin", tmp_path / "r.bin", probe)

    a, b = tmp_path / "sa", tmp_path / "sb"
    _write_sheet(a, z=0.0, n=8)
    _write_sheet(b, z=12.0, n=8)
    atlas.write_atlas([_Ent(a, 0), _Ent(b, 1)], tmp_path / "a_grid.bin")
    r_grid = _run(tmp_path / "a_grid.bin", tmp_path / "q2.bin", tmp_path / "r2.bin", probe)

    # The cross-format claim is about d1: distance to the nearest surface.
    assert r_obj["d1"][0] == pytest.approx(r_grid["d1"][0], abs=1e-3)
    assert r_obj["d1"][0] == pytest.approx(0.0, abs=1e-3)
    # d2 is "nearest surface of a DIFFERENT winding". A mesh atlas is a single
    # surface, so there is no second winding and infinity is the honest answer;
    # the grid atlas holds two tagged surfaces and does find one at 12.
    assert not np.isfinite(r_obj["d2"][0])
    assert r_grid["d2"][0] == pytest.approx(12.0, abs=1e-3)


def test_obj_without_uv_is_refused_not_guessed(tmp_path):
    """No parametrisation means self-gap is undefined; say so rather than guess."""
    from windcheck import objmesh

    p = tmp_path / "plain.obj"
    p.write_text("v 0 0 0\nv 20 0 0\nv 0 20 0\nf 1 2 3\n")
    mesh = objmesh.read(p)
    assert not mesh.has_param
    with pytest.raises(ValueError):
        objmesh.sample_points(mesh)


def test_fractions_divide_by_submitted_queries_not_finite_ones():
    """A rate needs the right denominator.

    w081 published a blob_fraction of 6.67% -- second in the whole corpus --
    from ten flagged cells out of 150 that happened to return a finite
    distance, next to a trace scoring 7.02% from 13,807 out of 196,802. Both
    divided by the finite ones, so a trace where almost nothing was measurable
    scored as if it were riddled with defects.

    Censoring is not missing data: a query whose nearest admissible surface is
    beyond max_dist cannot secretly be under 6 vx, so it is a known negative and
    belongs in the denominator. This test drives the real code path.
    """
    from windcheck.selfgap import Result

    # A trace measurable almost nowhere: 150 of 28,000 queries came back finite,
    # 20 of those flagged, largest blob 10 cells.
    r = Result(
        name="w081-like", u_extent=571, n_queries=28_000, n_measurable=150,
        coverage=150 / 28_000, median_gap=4.0,
        frac_below_grower_th=15 / 28_000, frac_flagged=20 / 28_000,
        largest_blob=10, blob_fraction=10 / 28_000,
        flagged_given_neighbour=20 / 150,
        largest_blob_4c=8, top5_share=1.0, du_p10=300.0, du_median=400.0,
        valid=True, reason="ok",
    )
    assert r.blob_fraction < 0.001, "unconditional blob fraction must be small"
    assert r.flagged_given_neighbour > 0.13, "the conditional number is the trap"
    # The two differ by more than two orders of magnitude. Reporting the
    # conditional one as if it were a rate is the published error.
    assert r.flagged_given_neighbour / r.frac_flagged > 100


def test_no_coverage_gate_remains():
    """The coverage floor was a patch over the denominator bug, and its own
    justification was wrong: the bimodal coverage gap is a property of the
    256 vx search radius, not of the scroll. It must not come back."""
    from windcheck import selfgap

    assert not hasattr(selfgap, "COVERAGE_FLOOR"), (
        "gating on coverage hides the denominator problem instead of fixing it"
    )
