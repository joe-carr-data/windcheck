"""Validate the exact self-intersection predicate on geometry we construct.

The whole reason this tool exists is that proximity is ambiguous and crossing is
not. So the tests that matter most are the NEGATIVE ones: a surface folded back
on itself 3 voxels away must produce zero crossings. If it does not, the tool
inherits exactly the ambiguity it was built to escape and the census is worthless.

Each case is a hairpin: the sheet runs out along x, turns, and comes back. The
two arms are far apart in the u parameter, so the adjacency exclusion does not
suppress them, which is precisely the nonlocal situation we care about.
"""
from __future__ import annotations

import json
import struct
import subprocess
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "engines" / "selfcross"

pytestmark = pytest.mark.skipif(
    not ENGINE.exists(),
    reason="build first: clang++ -O3 -std=c++17 -pthread -o engines/selfcross "
           "engines/selfcross.cpp",
)


def write_atlas(path: Path, points: np.ndarray, valid: np.ndarray) -> None:
    """Same WCAT container the distance engine reads, single surface."""
    rows, cols = valid.shape
    with path.open("wb") as fh:
        fh.write(b"WCAT")
        fh.write(struct.pack("<II", 1, 1))
        fh.write(struct.pack("<iII", -1, rows, cols))
        fh.write(np.ascontiguousarray(points, dtype="<f4").tobytes())
        fh.write(np.ascontiguousarray(valid, dtype=np.uint8).tobytes())


def run(path: Path, out: Path, exclude: int = 2, diagonal: int = 0) -> dict:
    r = subprocess.run(
        [str(ENGINE), str(path), str(out), "4", "40", str(exclude), str(diagonal)],
        capture_output=True, text=True, check=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def hairpin(arm: int = 60, width: int = 8, gap: float = 3.0,
            cross: bool = False, coplanar: bool = False,
            pitch: float = 6.0) -> tuple[np.ndarray, np.ndarray]:
    """A sheet that runs out, turns, and comes back `gap` voxels above itself.

    cross=True tilts the return arm through the outbound one, so a transverse
    intersection is guaranteed by construction.
    coplanar=True brings the return arm exactly onto the outbound plane.
    """
    cap = 10
    cols = arm + cap + arm
    pts = np.zeros((width, cols, 3), dtype=np.float32)
    for j in range(width):
        y = j * pitch
        for i in range(arm):                       # outbound, z = 0
            pts[j, i] = (i * pitch, y, 0.0)
        for k in range(cap):                       # semicircular turn
            th = np.pi * k / (cap - 1)
            cx = arm * pitch
            pts[j, arm + k] = (cx + np.sin(th) * gap * 0.5, y,
                               gap * 0.5 * (1 - np.cos(th)))
        for i in range(arm):                       # return
            x = (arm - 1 - i) * pitch
            if coplanar:
                z = 0.0
            elif cross:
                # sweep from +gap down through the outbound plane to -gap
                z = gap - 2.0 * gap * (i / max(1, arm - 1))
            else:
                z = gap
            pts[j, arm + cap + i] = (x, y, z)
    return pts, np.ones((width, cols), dtype=np.uint8)


def test_flat_sheet_has_no_crossings(tmp_path):
    """The null control. A plane cannot cross itself."""
    n = 40
    v, u = np.mgrid[0:n, 0:n]
    pts = np.stack([u * 6.0, v * 6.0, np.zeros_like(u, dtype=float)], -1)
    a, o = tmp_path / "a.bin", tmp_path / "o.csv"
    write_atlas(a, pts, np.ones((n, n), np.uint8))
    r = run(a, o)
    assert r["transverse"] == 0, f"flat sheet reported crossings: {r}"
    assert r["coplanar"] == 0, f"flat sheet reported coplanar overlap: {r}"


def test_tight_fold_three_voxels_apart_does_not_cross(tmp_path):
    """THE critical negative. Tight packing must not read as a crossing.

    Three voxels is closer than any real wrap spacing (calibration: 17.05) and
    well inside the 6 vx flag the proximity detector uses, so this is exactly
    the configuration Paul described. It must come back clean.
    """
    pts, valid = hairpin(gap=3.0, cross=False)
    a, o = tmp_path / "a.bin", tmp_path / "o.csv"
    write_atlas(a, pts, valid)
    r = run(a, o)
    assert r["transverse"] == 0, (
        f"tight but disjoint fold reported {r['transverse']} crossings; the "
        f"predicate is measuring proximity, not intersection: {r}")


def test_sub_voxel_fold_still_does_not_cross(tmp_path):
    """Half a voxel apart. Still no crossing, so still no report."""
    pts, valid = hairpin(gap=0.5, cross=False)
    a, o = tmp_path / "a.bin", tmp_path / "o.csv"
    write_atlas(a, pts, valid)
    r = run(a, o)
    assert r["transverse"] == 0, f"0.5 vx separation reported crossings: {r}"


def test_self_crossing_fold_is_detected(tmp_path):
    """The positive control: a return arm swept through the outbound plane."""
    pts, valid = hairpin(gap=6.0, cross=True)
    a, o = tmp_path / "a.bin", tmp_path / "o.csv"
    write_atlas(a, pts, valid)
    r = run(a, o)
    assert r["transverse"] > 0, f"planted crossing not found: {r}"
    rows = (tmp_path / "o.csv").read_text().strip().splitlines()
    assert len(rows) > 1
    assert any("transverse" in line for line in rows[1:])


def test_coincident_sheets_report_coplanar_not_transverse(tmp_path):
    """Two sheets flat against each other overlap without passing through.

    This is the case that must never be called a crossing, because a crushed
    scroll produces it legitimately.
    """
    pts, valid = hairpin(gap=6.0, coplanar=True)
    a, o = tmp_path / "a.bin", tmp_path / "o.csv"
    write_atlas(a, pts, valid)
    r = run(a, o)
    assert r["transverse"] == 0, (
        f"coincident coplanar sheets reported as transverse: {r}")


def test_adjacency_exclusion_is_load_bearing(tmp_path):
    """With exclusion off, neighbouring triangles share edges and report.

    This documents why the grid-index window exists: without it every adjacent
    pair grazes, and the census would be pure noise.
    """
    pts, valid = hairpin(gap=3.0, cross=False)
    a = tmp_path / "a.bin"
    write_atlas(a, pts, valid)
    off = run(a, tmp_path / "off.csv", exclude=-1)
    on = run(a, tmp_path / "on.csv", exclude=2)
    assert off["grazing"] + off["coplanar"] > on["grazing"] + on["coplanar"]


def test_both_diagonals_agree_on_a_clean_fold(tmp_path):
    """Tessellation sensitivity check.

    A real crossing should not depend on which way the quads are split. Running
    both diagonals is the diagnostic that separates genuine crossings from
    twisted-quad artifacts, so the two must agree on an unambiguous case.
    """
    pts, valid = hairpin(gap=3.0, cross=False)
    a = tmp_path / "a.bin"
    write_atlas(a, pts, valid)
    d0 = run(a, tmp_path / "d0.csv", diagonal=0)
    d1 = run(a, tmp_path / "d1.csv", diagonal=1)
    assert d0["transverse"] == d1["transverse"] == 0


def test_counts_are_deterministic_across_threads_and_cells(tmp_path):
    """The bug that invalidated the first census.

    A triangle pair can land in several spatial cells. With a per-thread `seen`
    set it was tested once per thread that received one of those cells, so the
    same trace returned 11,211 crossings on one thread and 21,213 on twelve and
    the number moved with scheduling. Each pair now has exactly one owning cell.

    Counts must therefore be identical across thread counts, and across cell
    sizes -- the grid is an acceleration structure, not a parameter of the
    answer.
    """
    pts, valid = hairpin(gap=6.0, cross=True)
    a = tmp_path / "a.bin"
    write_atlas(a, pts, valid)

    counts = set()
    for threads in (1, 2, 8):
        for cellsize in ("20", "40", "80"):
            r = subprocess.run(
                [str(ENGINE), str(a), str(tmp_path / f"o{threads}_{cellsize}.csv"),
                 str(threads), cellsize, "1", "0"],
                capture_output=True, text=True, check=True)
            j = json.loads(r.stdout.strip().splitlines()[-1])
            counts.add(j["transverse"])
    assert len(counts) == 1, (
        f"crossing count depends on threads or cell size: {sorted(counts)}")
    assert counts.pop() > 0, "positive control found nothing"


def test_grazing_count_is_also_deterministic(tmp_path):
    """Grazing must not track the broad-phase cell size either.

    It did: two triangles hundreds of voxels apart were classified GRAZING
    whenever one's vertex lay near the other's INFINITE plane, so a larger cell
    tested more distant pairs and manufactured more grazing (24,013 at cell 20
    against 266,337 at cell 80 on one real trace). An AABB separation test fixed
    it. This is the check that caught it.
    """
    pts, valid = hairpin(gap=3.0, cross=False)
    a = tmp_path / "a.bin"
    write_atlas(a, pts, valid)
    seen = set()
    for cellsize in ("20", "40", "80", "160"):
        r = subprocess.run(
            [str(ENGINE), str(a), str(tmp_path / f"g{cellsize}.csv"),
             "4", cellsize, "1", "0"],
            capture_output=True, text=True, check=True)
        seen.add(json.loads(r.stdout.strip().splitlines()[-1])["grazing"])
    assert len(seen) == 1, f"grazing count tracks cell size: {sorted(seen)}"


def test_endpoint_touch_is_not_a_crossing(tmp_path):
    """Intervals meeting at a single point are a touch, not a crossing.

    Two triangles can share a point of their intersection lines without either
    interior passing through the other. The first implementation reported that
    as TRANSVERSE, which overstates what was found.
    """
    # Two triangles meeting exactly along a shared plane line, no penetration.
    n, cols = 6, 40
    pts = np.zeros((n, cols, 3), dtype=np.float32)
    for j in range(n):
        for i in range(cols // 2):                      # flat arm at z=0
            pts[j, i] = (i * 6.0, j * 6.0, 0.0)
        for i in range(cols // 2, cols):                # arm that only touches
            k = i - cols // 2
            pts[j, i] = ((cols // 2 - 1) * 6.0, j * 6.0, k * 6.0)
    a, o = tmp_path / "a.bin", tmp_path / "o.csv"
    write_atlas(a, pts, np.ones((n, cols), np.uint8))
    r = run(a, o, exclude=1)
    assert r["transverse"] == 0, (
        f"a touching, non-penetrating join was reported as a crossing: {r}")


def test_penetration_is_reported_and_positive(tmp_path):
    """Every reported crossing must carry a margin, so marginal ones are visible."""
    pts, valid = hairpin(gap=6.0, cross=True)
    a, o = tmp_path / "a.bin", tmp_path / "o.csv"
    write_atlas(a, pts, valid)
    run(a, o, exclude=1)
    rows = [l.split(",") for l in o.read_text().strip().splitlines()[1:]]
    tv = [r for r in rows if r[4] == "transverse"]
    assert tv, "no transverse rows to check"
    pens = [float(r[5]) for r in tv]
    assert min(pens) > 0.0, "a transverse row reported non-positive penetration"
