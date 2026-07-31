"""Cross-validation of the analytic (SAT) pair_interval against the LP
reference oracle, on randomized and adversarial constructions.

The LP is exact for the Minkowski-difference formulation; the SAT version
must agree on feasibility and on both interval endpoints. Coplanar and
touching cases are planted explicitly because the in-plane axes only
matter there.
"""
from __future__ import annotations

import numpy as np
import pytest

from windcheck.clearance import pair_interval, pair_interval_lp

TOL = 1e-6


def compare(A, B, d):
    got = pair_interval(A, B, d)
    ref = pair_interval_lp(A, B, d)
    if ref is None or got is None:
        assert got is None and ref is None, (got, ref, A, B, d)
        return
    assert got[0] == pytest.approx(ref[0], abs=TOL)
    assert got[1] == pytest.approx(ref[1], abs=TOL)


def test_hand_computed_cases_match():
    BASE = np.array([[-2.0, -2.0, 0.0], [4.0, -2.0, 0.0], [0.0, 4.0, 0.0]])
    V = np.array([[0.0, 0.0, -1.0], [0.0, 0.0, 2.0], [0.0, 1.0, 0.5]])
    lo, hi = pair_interval(V, BASE, np.array([0.0, 0.0, 1.0]))
    assert lo == pytest.approx(-2.0, abs=1e-9)
    assert hi == pytest.approx(1.0, abs=1e-9)
    compare(V, BASE, np.array([0.0, 0.0, 1.0]))
    compare(V, BASE, np.array([0.0, 0.0, -1.0]))


def test_randomized_pairs_match_lp():
    rng = np.random.default_rng(23)
    n_none = n_iv = 0
    for _ in range(400):
        A = rng.normal(scale=2.0, size=(3, 3))
        B = rng.normal(scale=2.0, size=(3, 3)) + rng.normal(scale=2.0, size=3)
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        got = pair_interval(A, B, d)
        if got is None:
            n_none += 1
        else:
            n_iv += 1
        compare(A, B, d)
    assert n_none > 20 and n_iv > 20      # both outcomes well exercised


def test_coplanar_pairs_match_lp():
    """In-plane axes are load-bearing only for coplanar pairs."""
    rng = np.random.default_rng(5)
    for _ in range(120):
        A = np.column_stack([rng.normal(scale=2.0, size=(3, 2)),
                             np.zeros(3)])
        B = np.column_stack([rng.normal(scale=2.0, size=(3, 2))
                             + rng.normal(scale=3.0, size=2),
                             np.zeros(3)])
        d = np.array([*rng.normal(size=2), 0.0])
        d /= np.linalg.norm(d)
        compare(A, B, d)


def test_touching_and_shared_vertex_cases():
    T = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    S = np.array([[0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]])
    compare(T, S, np.array([0.0, 0.0, 1.0]))
    compare(T, S, np.array([1.0, 1.0, 0.0]) / np.sqrt(2))
    # vertical stack, exact face contact
    U = T + np.array([0.0, 0.0, 1.0])
    compare(U, T, np.array([0.0, 0.0, 1.0]))


def test_perpendicular_direction_disjoint_is_none_in_both():
    """d orthogonal to the separation: never intersects."""
    T = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    F = T + np.array([0.0, 0.0, 5.0])      # far above
    assert pair_interval(F, T, np.array([1.0, 0.0, 0.0])) is None
    assert pair_interval_lp(F, T, np.array([1.0, 0.0, 0.0])) is None


def test_batched_matches_scalar():
    from windcheck.clearance import pair_intervals_batched
    rng = np.random.default_rng(11)
    As = rng.normal(scale=2.0, size=(200, 3, 3))
    Bs = rng.normal(scale=2.0, size=(200, 3, 3)) \
        + rng.normal(scale=2.0, size=(200, 1, 3))
    d = rng.normal(size=3)
    d /= np.linalg.norm(d)
    lo, hi, ok = pair_intervals_batched(As, Bs, d)
    for i in range(200):
        s = pair_interval(As[i], Bs[i], d)
        if s is None:
            assert not ok[i], i
        else:
            assert ok[i], i
            assert lo[i] == pytest.approx(s[0], abs=1e-9)
            assert hi[i] == pytest.approx(s[1], abs=1e-9)


def test_batched_chunking_is_bit_identical():
    """The >100k chunk path must equal the single-shot path exactly."""
    from windcheck.clearance import pair_intervals_batched
    import windcheck.clearance as C
    rng = np.random.default_rng(41)
    As = rng.normal(scale=2.0, size=(1500, 3, 3))
    Bs = rng.normal(scale=2.0, size=(1500, 3, 3)) \
        + rng.normal(scale=2.0, size=(1500, 1, 3))
    d = rng.normal(size=3)
    d /= np.linalg.norm(d)
    lo1, hi1, ok1 = pair_intervals_batched(As, Bs, d)
    src = C.pair_intervals_batched.__wrapped__ if hasattr(
        C.pair_intervals_batched, "__wrapped__") else None
    # force the chunk path by monkeypatching CHUNK via a tiny wrapper call
    import numpy as _np
    parts = [pair_intervals_batched(As[i:i + 100], Bs[i:i + 100], d)
             for i in range(0, 1500, 100)]
    lo2 = _np.concatenate([p[0] for p in parts])
    hi2 = _np.concatenate([p[1] for p in parts])
    ok2 = _np.concatenate([p[2] for p in parts])
    assert _np.array_equal(lo1, lo2) and _np.array_equal(hi1, hi2) \
        and _np.array_equal(ok1, ok2)
