"""Per-voxel sheet labels for Scroll 1, and how to line them up with a mesh.

The Vesuvius Challenge publishes 80 hand-annotated 256^3 cubes for Scroll 1 at

    dl.ash2txt.org/full-scrolls/Scroll1/PHercParis4.volpkg/
        volumetric-instance-labels/instance-labels-harmonized/

Each voxel of papyrus carries an integer id naming the physical sheet it belongs
to, and the ids are global: 434 distinct ids in the range 1-489 appear across the
80 cubes, roughly 13 per cube, with a median 31% of voxels labelled. This is
ground truth of a kind nothing else in this project has had. Everything measured
so far is a property of the mesh; these labels are a statement about the scroll.

**They are usable directly.** Scroll 1 publishes every segment on
`20230205180739-7.91um`, which is the same volume the labels annotate, so no
registration or resampling is involved -- the mesh coordinates and the label
voxels are the same grid.

## The axis mapping, and how it was established rather than assumed

A `.nrrd` cube declares `space origin` equal to its directory name and an
identity `space directions`, but the file says `left-posterior-superior`, which
would ordinarily mean the axes are (x, y, z), while the coordinate ranges suggest
otherwise. Rather than pick one, all six permutations were tested against the
one thing that must be true if the mapping is right: **a mesh point lies on
papyrus**, so it should land on a labelled voxel far more often than a random
point in the same cube does.

    permutation      cubes overlapped   exact hit   within 2 vx   enrichment
    (z, y, x)              69             0.437        0.699        1.43x
    (z, x, y)              32             0.321        0.492        1.03x
    (y, x, z)              13             0.286        0.441        0.99x

The two alternatives sit exactly at chance. `(z, y, x)` is the mapping, and it
overlaps twice as many cubes, which is a second signature of being right.

**The enrichment is 1.43x, not 10x, and that is not a bug to hide.** Only 44% of
mesh points from these segments land exactly on labelled papyrus. Some of that is
label coverage and sub-voxel placement -- 70% land within two voxels. Some of it
may be tracing error, which is the thing under study. So the agreement rate is
reported per segment as a measurement in its own right, and no comparison in this
module relies on the global rate being high: the tests that matter are
*within* a segment, comparing flagged locations against unflagged ones, where a
constant registration quality cancels.
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

import numpy as np

# Established empirically above: mesh coordinate axes (x, y, z) index the cube
# array as (z, y, x). Written as the permutation applied to a mesh point.
MESH_TO_CUBE = (2, 1, 0)
CUBE = 256
DEFAULT_ROOT = Path("data/scroll1_labels")
_NAME_RE = re.compile(r"(\d+)_(\d+)_(\d+)_mask\.nrrd$")


class LabelCubes:
    """The 80 cubes, indexed by origin, loaded on demand.

    Cubes are ~800 KB compressed and 32 MB expanded, so they are read lazily and
    cached; a full corpus pass touches each one a handful of times.
    """

    def __init__(self, root: Path = DEFAULT_ROOT):
        self.origins: dict[tuple[int, int, int], Path] = {}
        for f in sorted(glob.glob(str(Path(root) / "*_mask.nrrd"))):
            m = _NAME_RE.search(Path(f).name)
            if m:
                self.origins[tuple(int(g) for g in m.groups())] = Path(f)
        self._cache: dict[tuple[int, int, int], np.ndarray] = {}

    def __len__(self) -> int:
        return len(self.origins)

    def data(self, origin) -> np.ndarray:
        if origin not in self._cache:
            import nrrd
            self._cache[origin] = nrrd.read(str(self.origins[origin]))[0]
        return self._cache[origin]

    def lookup(self, points: np.ndarray) -> np.ndarray:
        """Sheet id for each mesh point; 0 where unlabelled or outside a cube.

        `points` are (N, 3) mesh coordinates in the 7.91 um volume, exactly as
        `tifxyz.read` returns them. Points outside every cube get 0, which the
        caller must distinguish from "inside a cube but unlabelled" -- use
        `covered` for that, because conflating the two turns an unasked question
        into a negative answer, and that error has already cost this project a
        published claim.
        """
        q = np.asarray(points)[:, MESH_TO_CUBE]
        out = np.zeros(len(q), dtype=np.int32)
        for origin in self.origins:
            lo = np.array(origin)
            sel = np.all((q >= lo) & (q < lo + CUBE), axis=1)
            if not sel.any():
                continue
            idx = (q[sel] - lo).astype(int)
            out[sel] = self.data(origin)[idx[:, 0], idx[:, 1], idx[:, 2]]
        return out

    def covered(self, points: np.ndarray) -> np.ndarray:
        """True where the point falls inside some annotated cube at all."""
        q = np.asarray(points)[:, MESH_TO_CUBE]
        out = np.zeros(len(q), dtype=bool)
        for origin in self.origins:
            lo = np.array(origin)
            out |= np.all((q >= lo) & (q < lo + CUBE), axis=1)
        return out
