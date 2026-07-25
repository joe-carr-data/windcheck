"""Reader for the `tifxyz` surface-mesh format used by VC3D.

A tifxyz surface is a 2D grid of 3D vertex coordinates: the (u, v) array shape
gives the mesh's grid topology, and the values give its shape in space. It is
stored as three float32 TIFFs -- `x.tif`, `y.tif`, `z.tif` -- plus a `meta.json`
carrying `bbox` and `scale`.

Missing cells are marked either by setting each coordinate to -1, or via an
optional sidecar mask image. Both are handled here; a cell is invalid if either
convention says so.

Deliberately implemented from the format description rather than by calling the
upstream `vesuvius` package. This tool exists to disagree with that pipeline
when the pipeline is wrong, so it must not inherit its reader's assumptions.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile

MISSING = -1.0


@dataclass
class Surface:
    """A tifxyz surface: (u, v) grid -> (x, y, z) volume coordinates."""

    points: np.ndarray  # float32, shape (v, u, 3)
    valid: np.ndarray  # bool, shape (v, u)
    scale: tuple[float, float]
    bbox: tuple[tuple[float, float, float], tuple[float, float, float]]
    path: Path

    @property
    def shape(self) -> tuple[int, int]:
        """(v, u) — the grid dimensions, in numpy row/col order."""
        return self.points.shape[:2]

    @property
    def n_valid(self) -> int:
        return int(self.valid.sum())

    @property
    def coverage(self) -> float:
        return self.n_valid / self.valid.size if self.valid.size else 0.0

    def xyz(self) -> np.ndarray:
        """All valid points as an (N, 3) array."""
        return self.points[self.valid]


def read(directory: str | Path) -> Surface:
    """Read a tifxyz surface from a directory containing x/y/z.tif + meta.json."""
    d = Path(directory)
    arrays = []
    for axis in ("x", "y", "z"):
        f = d / f"{axis}.tif"
        if not f.exists():
            raise FileNotFoundError(f"tifxyz missing {axis}.tif: {d}")
        arrays.append(np.asarray(tifffile.imread(f), dtype=np.float32))

    shapes = {a.shape for a in arrays}
    if len(shapes) != 1:
        raise ValueError(f"x/y/z grids disagree in shape at {d}: {shapes}")

    points = np.stack(arrays, axis=-1)  # (v, u, 3)

    # Convention 1: a cell is missing when every coordinate is -1.
    valid = ~np.all(points == MISSING, axis=-1)
    # NaNs are never usable, whatever the convention says.
    valid &= np.all(np.isfinite(points), axis=-1)

    # Convention 2: optional sidecar mask. Non-zero = keep.
    for name in ("mask.tif", "mask.png"):
        m = d / name
        if m.exists():
            mask = np.asarray(tifffile.imread(m)) if m.suffix == ".tif" else None
            if mask is not None and mask.shape == valid.shape:
                valid &= mask.astype(bool)
            break

    meta_path = d / "meta.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    scale = tuple(meta.get("scale", (1.0, 1.0)))
    bbox = meta.get("bbox", [[0.0] * 3, [0.0] * 3])

    return Surface(
        points=points,
        valid=valid,
        scale=(float(scale[0]), float(scale[1])),
        bbox=(tuple(bbox[0]), tuple(bbox[1])),  # type: ignore[arg-type]
        path=d,
    )
