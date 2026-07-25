"""The winding reference atlas, and the bridge to the C++ query engine.

The 44 single-winding segments (`w052`-`w095`, contiguous) say where each sheet
physically is. Together they form an atlas: any 3D point can be assigned the
winding of the nearest reference surface, with a confidence taken from how much
closer that surface is than the nearest surface of a *different* winding.

Heavy lifting is in `engines/atlas_query.cpp` -- ~17M query points against
~33M reference triangles is not a Python workload. Python owns IO, assembly and
reporting; C++ owns the kernel. Data crosses as flat binary.
"""

from __future__ import annotations

import re
import struct
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import tifxyz

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / "engines" / "atlas_query"
_WINDING_RE = re.compile(r"-w(\d{3})_")

MAGIC_ATLAS = b"WCAT"
MAGIC_QUERY = b"WCQP"
VERSION = 1

# The volume-space variant to read. Both published variants are the same
# geometry in one frame (catalog gives overlap_ratio 1.0 and identical bboxes),
# so this is a choice of encoding, not of data. 839 is stored uncompressed.
VOLUME = "20241024131839"


@dataclass
class Entry:
    winding: int | None
    long_id: str
    path: Path

    @property
    def is_auto_grown(self) -> bool:
        return self.winding is None


def discover(data_dir: Path) -> list[Entry]:
    """Find every downloaded segment and its volume-space tifxyz directory."""
    out: list[Entry] = []
    for seg in sorted(data_dir.iterdir()):
        if not seg.is_dir():
            continue
        matches = sorted(seg.glob(f"mesh/*-on-{VOLUME}-*.tifxyz"))
        if not matches:
            continue
        m = _WINDING_RE.search(seg.name)
        out.append(Entry(int(m.group(1)) if m else None, seg.name, matches[0]))
    out.sort(key=lambda e: (e.winding is None, e.winding if e.winding is not None else 0))
    return out


def write_atlas(entries: list[Entry], path: Path) -> int:
    """Serialise reference surfaces for the engine. Returns surface count."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        fh.write(MAGIC_ATLAS)
        fh.write(struct.pack("<II", VERSION, len(entries)))
        for e in entries:
            s = tifxyz.read(e.path)
            rows, cols = s.shape
            fh.write(struct.pack("<iII", e.winding if e.winding is not None else -1, rows, cols))
            fh.write(np.ascontiguousarray(s.points, dtype="<f4").tobytes())
            fh.write(np.ascontiguousarray(s.valid, dtype=np.uint8).tobytes())
    return len(entries)


def write_queries(points: np.ndarray, path: Path) -> int:
    """Serialise an (N, 3) float array of query points."""
    pts = np.ascontiguousarray(points, dtype="<f4").reshape(-1, 3)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        fh.write(MAGIC_QUERY)
        fh.write(struct.pack("<II", VERSION, len(pts)))
        fh.write(pts.tobytes())
    return len(pts)


def write_queries_grouped(points: np.ndarray, groups: np.ndarray, path: Path) -> int:
    """Serialise query points tagged with a group id (self-gap mode).

    The group is the point's u-index along the trace. The engine refuses any
    triangle whose own u-index is within `exclude_u` of it, so what it measures
    is the distance to the trace's *neighbouring wrap* rather than to the patch
    of surface the point is already sitting on.
    """
    pts = np.ascontiguousarray(points, dtype="<f4").reshape(-1, 3)
    grp = np.ascontiguousarray(groups, dtype="<i4").reshape(-1)
    if len(grp) != len(pts):
        raise ValueError(f"group count {len(grp)} != point count {len(pts)}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        fh.write(b"WCQ2")
        fh.write(struct.pack("<II", VERSION, len(pts)))
        fh.write(pts.tobytes())
        fh.write(grp.tobytes())
    return len(pts)


RESULT_DTYPE = np.dtype([("w1", "<i4"), ("d1", "<f4"), ("w2", "<i4"), ("d2", "<f4")])


def read_results(path: Path) -> np.ndarray:
    return np.fromfile(path, dtype=RESULT_DTYPE)


def run_engine(atlas: Path, queries: Path, out: Path, threads: int = 0,
               cell: float = 32.0, max_dist: float = 256.0,
               exclude_u: int = 0, quiet: bool = True) -> np.ndarray:
    """Invoke the C++ engine and return its results.

    `exclude_u > 0` switches on self-gap mode, which also requires the query
    file to have been written by `write_queries_grouped`.
    """
    if not ENGINE.exists():
        raise FileNotFoundError(
            f"engine not built: {ENGINE}\n"
            "  clang++ -O3 -std=c++17 -pthread -o engines/atlas_query engines/atlas_query.cpp"
        )
    cmd = [str(ENGINE), str(atlas), str(queries), str(out),
           str(threads or 0), str(cell), str(max_dist), str(exclude_u)]
    subprocess.run(cmd, check=True,
                   stderr=subprocess.DEVNULL if quiet else None)
    return read_results(out)


def sample_points(entry: Entry, stride: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Valid surface points of a segment, with their (v, u) grid indices.

    `stride` subsamples the grid. Winding index varies along `u` (around the
    scroll) and should be constant in `v` (up the scroll), so `v` can be
    subsampled hard without losing the signal.
    """
    s = tifxyz.read(entry.path)
    valid = s.valid[::stride, ::stride]
    pts = s.points[::stride, ::stride]
    vi, ui = np.nonzero(valid)
    return pts[vi, ui], np.stack([vi * stride, ui * stride], axis=1)
