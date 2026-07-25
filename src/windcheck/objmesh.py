"""Reader for Wavefront OBJ triangle meshes.

The check was built on `tifxyz`, VC3D's native grid format, but the requirement
is to accept standard community formats. Every published segment also ships
`.obj`, and those carry texture coordinates -- faces are indexed `v/vt`, so a
2D parametrisation is present.

That parametrisation is what makes self-gap well defined on an unstructured
mesh. On a `tifxyz` grid the exclusion is "at least N columns away in u"; here
it is the same statement against the `vt` u-coordinate, rescaled to comparable
units so a caller's `exclude_u` means the same thing in both formats.

A mesh with no `vt` block can still be used as a reference surface (how close
does some other surface come to it?), but not for self-gap, because there is no
coordinate along which to exclude a point's own neighbourhood. That case raises
rather than guessing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Mesh:
    """A triangle mesh with an optional 2D parametrisation."""

    verts: np.ndarray          # (V, 3) float32, volume coordinates
    tris: np.ndarray           # (T, 3) int32, indices into verts
    tri_u: np.ndarray | None   # (T,) int32, per-triangle u tag, or None

    @property
    def n_tris(self) -> int:
        return len(self.tris)

    @property
    def has_param(self) -> bool:
        return self.tri_u is not None

    def triangle_xyz(self) -> np.ndarray:
        """(T, 3, 3) explicit triangle corners."""
        return self.verts[self.tris]


def read(path: str | Path, u_scale: float | None = None) -> Mesh:
    """Parse an OBJ file.

    `u_scale` maps the `vt` u-coordinate onto integer tags. If omitted, it is
    chosen so the full u-range spans the same number of tags as a tifxyz grid
    of the same mesh would have columns, which keeps `exclude_u` meaningful
    across formats.
    """
    path = Path(path)
    verts: list[tuple[float, float, float]] = []
    uvs: list[float] = []
    faces: list[tuple[int, int, int]] = []
    face_uv: list[tuple[int, int, int]] = []

    with path.open() as fh:
        for line in fh:
            if not line or line[0] not in "vf":
                continue
            parts = line.split()
            if not parts:
                continue
            tag = parts[0]
            if tag == "v":
                verts.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif tag == "vt":
                uvs.append(float(parts[1]))
            elif tag == "f":
                idx, uvidx = [], []
                for token in parts[1:4]:
                    bits = token.split("/")
                    idx.append(int(bits[0]) - 1)
                    if len(bits) > 1 and bits[1]:
                        uvidx.append(int(bits[1]) - 1)
                if len(idx) == 3:
                    faces.append(tuple(idx))
                    if len(uvidx) == 3:
                        face_uv.append(tuple(uvidx))

    if not faces:
        raise ValueError(f"no triangular faces in {path}")

    V = np.asarray(verts, dtype=np.float32)
    T = np.asarray(faces, dtype=np.int32)

    tri_u = None
    if uvs and len(face_uv) == len(faces):
        U = np.asarray(uvs, dtype=np.float32)
        fu = U[np.asarray(face_uv, dtype=np.int32)].mean(axis=1)   # per-triangle u
        lo, hi = float(fu.min()), float(fu.max())
        span = max(hi - lo, 1e-9)
        if u_scale is None:
            # Match a tifxyz grid's column pitch: those are ~20 voxels apart, so
            # a mesh of this physical width would carry roughly this many columns.
            width = float(np.linalg.norm(V.max(axis=0) - V.min(axis=0)))
            u_scale = max(width / 20.0, 1.0) / span
        tri_u = np.round((fu - lo) * u_scale).astype(np.int32)

    return Mesh(verts=V, tris=T, tri_u=tri_u)


def sample_points(mesh: Mesh, stride: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Triangle centroids and their u tags, as query points.

    Centroids rather than vertices: a vertex is shared between triangles that
    may carry different u tags, which would make its exclusion ambiguous.
    """
    if not mesh.has_param:
        raise ValueError("mesh has no vt parametrisation; self-gap is undefined on it")
    cent = mesh.triangle_xyz()[::stride].mean(axis=1)
    return cent.astype(np.float32), mesh.tri_u[::stride]
