"""Pins for the W1 shadow-measurement harness (bench/excise_shadow.py).

The shadow runs themselves need real data and the C++ engine; what these
tests pin is the MEASUREMENT logic, because a measurement instrument that
is wrong is worse than no measurement:

1. `component_recovery` -- round-26 Q3's fragmentation statistic. The
   headline failure mode is a hairline cut that scores ~100% on area while
   splitting a sheet in half, so the 50/50 split is pinned explicitly.
2. `cut_boundary_arrays` -- the vectorised boundary must agree EXACTLY
   with the scalar `cut_boundary` the certificate path uses.
3. `pick_input_mesh` -- the displacement-repaired mesh is preferred, the
   original is the labelled fallback.
"""
from __future__ import annotations

import sys
from math import fsum

import numpy as np

sys.path.insert(0, "bench")
import excise_shadow as sh                                      # noqa: E402
from excise_segment import cut_boundary, quad_area_grids        # noqa: E402


def test_component_recovery_uncut_is_identity():
    Q = np.ones((4, 6), bool)
    A = np.ones((4, 6))
    r = sh.component_recovery(Q, Q, A)
    assert r["components_before"] == r["components_after"] == 1
    assert r["min_R_main_any_component"] == 1.0
    assert r["area_weighted_R_main"] == 1.0
    assert r["per_component"][0]["retained_fraction"] == 1.0
    assert r["per_component"][0]["n_descendants"] == 1
    assert r["n_input_components_fully_destroyed"] == 0


def test_component_recovery_hairline_split_is_caught_by_R_main():
    """A cut that removes ONE column out of 21 keeps 95.2% of the area but
    halves the sheet: area says 0.952, R_main says 0.476. That gap is the
    entire reason Q3 made R_main the primary statistic."""
    Q_in = np.ones((4, 21), bool)
    Q_out = Q_in.copy()
    Q_out[:, 10] = False
    A = np.ones((4, 21))
    r = sh.component_recovery(Q_in, Q_out, A)
    assert r["components_before"] == 1
    assert r["components_after"] == 2
    retained = r["per_component"][0]["retained_fraction"]
    assert abs(retained - 80 / 84) < 1e-12          # ~0.952 on area
    assert abs(r["min_R_main_any_component"] - 40 / 84) < 1e-12
    assert r["n_components_R_main_below_0.9"] == 1
    assert r["per_component"][0]["n_descendants"] == 2
    assert r["per_component"][0]["fragments_over_1pct"] == 2


def test_component_recovery_fully_destroyed_component():
    Q_in = np.zeros((4, 9), bool)
    Q_in[:, :3] = True                              # big component
    Q_in[:, 6:] = True                              # second component
    Q_out = Q_in.copy()
    Q_out[:, 6:] = False                            # the second one is gone
    A = np.ones((4, 9))
    r = sh.component_recovery(Q_in, Q_out, A)
    assert r["components_before"] == 2
    assert r["components_after"] == 1
    assert r["n_input_components_fully_destroyed"] == 1
    assert r["min_R_main_any_component"] == 0.0
    assert r["area_weighted_R_main"] == 0.5


def test_cut_boundary_arrays_matches_the_scalar_certificate_path():
    rng = np.random.default_rng(7)
    P = rng.normal(size=(9, 11, 3))
    Q = np.ones((8, 10), bool)
    Q[3, 3] = False                                 # a pre-existing hole
    removed = np.zeros_like(Q)
    removed[2:5, 4:6] = True
    removed &= Q
    kept = Q & ~removed
    n, length = sh.cut_boundary_arrays(P, removed, kept)

    rset = {(int(v), int(u)) for v, u in zip(*np.nonzero(removed))}
    kset = {(int(v), int(u)) for v, u in zip(*np.nonzero(kept))}
    edges = cut_boundary(rset, kset)
    assert n == len(edges)
    scalar = fsum(float(np.linalg.norm(P[a] - P[b])) for a, b in edges)
    assert abs(length - scalar) < 1e-9


def test_cut_boundary_arrays_zero_when_nothing_is_cut():
    P = np.zeros((5, 5, 3))
    Q = np.ones((4, 4), bool)
    assert sh.cut_boundary_arrays(P, np.zeros_like(Q), Q) == (0, 0.0)


def test_pick_input_mesh_prefers_the_repaired_mesh(tmp_path, monkeypatch):
    monkeypatch.setattr(sh, "REPAIRED", tmp_path / "meshes")
    original = tmp_path / "orig.tifxyz"
    original.mkdir()
    got = sh.pick_input_mesh("segX", original)
    assert got["source"] == "original_published"
    assert got["path"] == original
    rep = tmp_path / "meshes" / "segX_repaired.tifxyz"
    rep.mkdir(parents=True)
    (rep / "x.tif").write_bytes(b"")
    got = sh.pick_input_mesh("segX", original)
    assert got["source"] == "displacement_repaired"
    assert got["path"] == rep


def test_quad_area_grids_matches_a_hand_computed_unit_square():
    P = np.zeros((2, 2, 3))
    P[0, 1] = (0, 1, 0)
    P[1, 0] = (1, 0, 0)
    P[1, 1] = (1, 1, 0)
    a = quad_area_grids(P)
    assert abs(float(a[0][0, 0]) - 1.0) < 1e-12
    assert abs(float(a[1][0, 0]) - 1.0) < 1e-12
    assert abs(float(a["canonical"][0, 0]) - 1.0) < 1e-12
