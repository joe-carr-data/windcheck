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
    # The upstream loader's stricter rule, carried alongside rather than
    # substituted for `valid`. See `read` for why both exist.
    valid_pipeline: np.ndarray | None = None

    @property
    def shape(self) -> tuple[int, int]:
        """(v, u) — the grid dimensions, in numpy row/col order."""
        return self.points.shape[:2]

    @property
    def n_valid(self) -> int:
        return int(self.valid.sum())

    @property
    def n_valid_pipeline(self) -> int:
        """Cells the upstream loader would keep. Never larger than n_valid."""
        if self.valid_pipeline is None:
            return self.n_valid
        return int(self.valid_pipeline.sum())

    @property
    def z_floor_cells(self) -> int:
        """Cells we keep that `QuadSurface` discards for z <= 0.

        Zero on most surfaces. Where it is not, it is worth reporting: a
        crossing confined to these cells is one the pipeline never sees.
        """
        return self.n_valid - self.n_valid_pipeline

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

    # Convention 3, upstream only: `QuadSurface::load` rewrites every point
    # with z <= 0 to the -1 sentinel BEFORE it applies the mask, so its valid
    # set is a strict subset of the one above. We do not adopt that rule as
    # the default, because every published count was measured without it and
    # silently tightening a definition is how numbers stop meaning what they
    # said. We carry it alongside instead, so a certificate can state both and
    # a reader can see the size of the difference. Measured across the pinned
    # corpus the gap is 44 of 185 meshes and 1.49% of transverse rows, and it
    # changes no segment's verdict; see bench/zfloor_impact.py.
    valid_pipeline = valid & (points[..., 2] > 0)

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
        valid_pipeline=valid_pipeline,
    )


def write_meta(src: str | Path, dst: str | Path,
               points: np.ndarray, valid: np.ndarray) -> None:
    """Carry `meta.json` forward with a bbox recomputed from what was emitted.

    Every field but `bbox` is preserved verbatim: scale, uuid and anything
    else the pipeline put there is not ours to reinterpret. `bbox` is the
    one field that describes the points rather than the surface's identity,
    so it is the one field a writer must not inherit.
    """
    src, dst = Path(src), Path(dst)
    meta_path = src / "meta.json"
    meta = (json.loads(meta_path.read_text()) if meta_path.exists()
            else {"format": "tifxyz", "scale": [1, 1]})
    meta["bbox"] = bbox_of(points, valid)
    (dst / "meta.json").write_text(json.dumps(meta, indent=1) + "\n")


def bbox_of(points: np.ndarray, valid: np.ndarray) -> list[list[float]]:
    """[[xmin,ymin,zmin],[xmax,ymax,zmax]] over the valid points.

    `meta.json` carries a bbox that consumers filter on, and a stale one
    silently drops a surface from their inputs rather than failing loudly.
    An emitted mesh must therefore carry a bbox computed from the points
    it actually contains, never one inherited from the mesh it came from.
    """
    if not valid.any():
        return [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    P = points[valid]
    return [[float(v) for v in P.min(axis=0)],
            [float(v) for v in P.max(axis=0)]]
