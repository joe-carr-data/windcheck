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
