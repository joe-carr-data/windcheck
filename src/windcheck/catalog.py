"""Read the public Vesuvius open-data catalog.

The catalog is a single gzipped JSON blob at the bucket root. It is the only
authoritative index of samples -> volumes -> segments, and it is small (~1.1 MB
compressed), so we fetch it whole and cache it.
"""

from __future__ import annotations

import gzip
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from . import CATALOG_KEY, S3_BUCKET

_WINDING_RE = re.compile(r"-w(\d{3})_")


@dataclass(frozen=True)
class Segment:
    """One published segmentation of a scroll."""

    id: str
    long_id: str
    sample_id: str
    winding: int | None  # parsed from the `wNNN` name; None for auto_grown traces
    grid: tuple[int, int]  # tifxyz (u, v) dimensions
    bbox: tuple[tuple[float, float, float], tuple[float, float, float]]
    original_volume_id: str

    @property
    def is_auto_grown(self) -> bool:
        return "auto_grown" in self.long_id

    @property
    def s3_prefix(self) -> str:
        return f"{self.sample_id}/segments/{self.long_id}/"


def fetch_catalog(dest: Path) -> Path:
    """Download the gzipped catalog to `dest` if not already present."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        subprocess.run(
            ["aws", "s3", "cp", f"s3://{S3_BUCKET}/{CATALOG_KEY}", str(dest),
             "--no-sign-request", "--only-show-errors"],
            check=True,
        )
    return dest


def load_segments(catalog_path: Path, sample_id: str) -> list[Segment]:
    """Parse all segments for one sample, sorted by winding number.

    `auto_grown` traces carry no winding label and sort last.
    """
    with gzip.open(catalog_path, "rt") as fh:
        catalog = json.load(fh)

    out: list[Segment] = []
    for entry in catalog["samples"][sample_id]["segments"].values():
        long_id = entry["long_id"]
        meta = entry["creation"].get("metadata", {})
        match = _WINDING_RE.search(long_id)
        grid = meta.get("tiff_dimensions") or [0, 0]
        bbox = meta.get("bbox") or [[0.0] * 3, [0.0] * 3]
        out.append(
            Segment(
                id=entry["id"],
                long_id=long_id,
                sample_id=entry["sample_id"],
                winding=int(match.group(1)) if match else None,
                grid=(int(grid[0]), int(grid[1])),
                bbox=(tuple(bbox[0]), tuple(bbox[1])),  # type: ignore[arg-type]
                original_volume_id=entry.get("original_volume_id", ""),
            )
        )

    out.sort(key=lambda s: (s.winding is None, s.winding if s.winding is not None else 0))
    return out


def winding_run(segments: list[Segment]) -> tuple[int, int, list[int]]:
    """Return (lo, hi, gaps) for the labelled winding numbers.

    A gap means a missing sheet between lo and hi, which weakens the atlas: the
    monotone-index argument needs a contiguous reference run.
    """
    windings = sorted(s.winding for s in segments if s.winding is not None)
    if not windings:
        return (0, 0, [])
    lo, hi = windings[0], windings[-1]
    gaps = [w for w in range(lo, hi + 1) if w not in set(windings)]
    return lo, hi, gaps
