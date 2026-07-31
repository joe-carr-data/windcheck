"""Emit results in a form VC3D can open, not just a table a human can read.

An audit that ends in a report is something people cite. An audit that ends in a
file the existing tool loads is something people use, and use is the criterion
that matters most here.

Two outputs:

**A point collection** in volume-cartographer's own `PointCollections` JSON
schema, so it drops into VC3D's point-collection widget and each overlap becomes
a clickable location. Keys taken from `core/src/PointCollections.cpp`: a
collection carries `name`, `color`, `metadata` and `points`; a point carries `p`
as [x,y,z] in volume coordinates.

**A certificate** summarising the trace in physical units, with enough
provenance that a reader can reproduce it: input hash, engine parameters, and
the counts each threshold produced.

Coordinates are volume voxels of the volume named in the mesh path, which is what
VC3D navigates in, so no transform is needed at the far end.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Sequence

import numpy as np

from . import classify

SCHEMA_NOTE = ("volume-cartographer PointCollections JSON; keys from "
               "core/src/PointCollections.cpp")


def _sha256_head(path: Path, limit: int = 8 << 20) -> str:
    """Cheap PARTIAL fingerprint of the input -- not a content hash.

    A `.tifxyz` surface is a DIRECTORY of x/y/z tiffs plus meta.json, not a
    file, so hashing the path directly raises IsADirectoryError. For a directory
    the member names, sizes and a bounded head of each are folded in, which pins
    the content without reading hundreds of megabytes.

    The bound is the point and also the limitation: two surfaces differing only
    past the head collide. The certificate says so in its caveats rather than
    letting the field read as a full hash. Callers wanting a true content hash
    should read the whole file; that is a deliberate trade, not an oversight.
    """
    h = hashlib.sha256()
    if path.is_dir():
        for f in sorted(p for p in path.iterdir() if p.is_file()):
            h.update(f.name.encode())
            h.update(str(f.stat().st_size).encode())
            with f.open("rb") as fh:
                h.update(fh.read(min(limit, 1 << 20)))
    else:
        with path.open("rb") as fh:
            h.update(fh.read(limit))
    return h.hexdigest()[:16]


# `PointCollections::loadFromJSON` refuses any file whose top level does not
# carry this key with this exact value, and then reads the collections out of
# `collections`. Writing the bare `{"1": {...}}` map -- which is what the
# collections VALUE looks like -- parses as JSON and loads as nothing: the
# loader throws "incorrect version or missing version info" and returns false.
# An overlay nobody can open is worth less than no overlay at all, so the
# envelope is asserted by a test rather than trusted.
PC_VERSION_KEY = "vc_pointcollections_json_version"
PC_VERSION = "1"


def point_collection(name: str, points: Sequence[Sequence[float]],
                     color: Sequence[float] = (1.0, 0.2, 0.2),
                     metadata: dict | None = None) -> dict:
    """One VC3D-loadable collection. `points` are [x, y, z] volume voxels."""
    pts = {}
    for i, p in enumerate(points):
        pts[str(i + 1)] = {
            "id": i + 1,
            "p": [float(p[0]), float(p[1]), float(p[2])],
            "wind_a": None,
            "collection_id": 1,
        }
    return {
        PC_VERSION_KEY: PC_VERSION,
        "collections": {
            "1": {
                "id": 1,
                "name": name,
                "color": [float(c) for c in color],
                "points": pts,
                "metadata": metadata or {},
            }
        },
    }


def write_collection(path: Path, name: str, points, **kw) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(point_collection(name, points, **kw), indent=1))
    return len(points)


# Kept as a name other modules import; it is a LITERAL filter boundary on a
# measured separation, not a claim about what a segment is. See classify.py.
WRAP_SCALE_CUT_REV = classify.SEP_WIDE


def certificate(*, segment: str, mesh_path: Path, voxel_um: float,
                params: dict, total_area_mm2: float,
                participating_quad_area_mm2: float,
                parameter_span_mm_estimate: float,
                separation_revolutions: float | None,
                covering_span_revolutions: float | None,
                revolution_period_columns: dict,
                period_status: str = "unavailable",
                n_pairs: int, n_events: int,
                events_beyond_cut: int, median_penetration_vx: float,
                extra: dict | None = None) -> dict:
    """Everything needed to reproduce and to judge the result.

    Field names say what was measured rather than what one would like it to
    mean. A transverse intersection is a *curve*, with zero two-dimensional
    area, so there is no "duplicated area" to report -- what the number sums is
    the area of the quads that take part. Likewise the span is an estimate along
    the parameterisation, not an intrinsic geodesic distance on the surface.
    """
    return {
        "tool": "windcheck selfcross",
        "schema": 2,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "segment": segment,
        "mesh": {"path": str(mesh_path), "sha256_head": _sha256_head(mesh_path),
                 "voxel_um": voxel_um},
        "parameters": params,
        "measurements": {
            "total_area_mm2": round(total_area_mm2, 1),
            "participating_quad_area_mm2": round(participating_quad_area_mm2, 2),
            "participating_quad_fraction": (
                round(participating_quad_area_mm2 / total_area_mm2, 6)
                if total_area_mm2 else 0.0),
            "parameter_span_mm_estimate": round(parameter_span_mm_estimate, 1),
            # The unit that ports between scrolls. Millimetres do not: voxel
            # size, scroll diameter and crushing all differ.
            "separation_revolutions": (round(separation_revolutions, 3)
                                       if separation_revolutions is not None
                                       else None),
            "separation_band": classify.separation_band(separation_revolutions),
            "crossing_status": classify.crossing_status(n_pairs),
            "period_status": period_status,
            "covering_span_revolutions": (round(covering_span_revolutions, 3)
                                          if covering_span_revolutions is not None
                                          else None),
            "revolution_period_columns": revolution_period_columns,
            "median_penetration_vx": round(median_penetration_vx, 3),
            "intersecting_pairs": int(n_pairs),
            "crossing_events": int(n_events),
            "events_beyond_cut": int(events_beyond_cut),
        },
        "separation_filter_revolutions": {"nearby": classify.SEP_NEARBY,
                                          "wide": classify.SEP_WIDE},
        "verdict": classify.verdict(classify.crossing_status(n_pairs),
                                    period_status, separation_revolutions,
                                    covering_span_revolutions),
        "note": classify.describe(separation_revolutions,
                                  covering_span_revolutions),
        # Stated on the artifact itself so it cannot be separated from the number.
        "caveats": [
            "Deterministic floating-point validator, not exact predicates.",
            "Detects that the surface overlaps itself; does not identify which "
            "branch is wrong, nor establish a cause.",
            "Event counts depend on the clustering rule and on sampling "
            "density; the span in revolutions does not.",
            "Local overlaps are common across the published corpus and are "
            "associated with elevated quad twist.",
            "The revolution period is measured from the surface's own geometry "
            "and carries its own error; two independent estimators agree to "
            "within about 15% where both are measurable.",
            "Separations are reported as measured distances. No band in this "
            "output asserts a cause; three attempts to connect crossings to "
            "sheet misassignment did not succeed.",
            "A segment covering more than one revolution has two ends in the "
            "same angular sector and can meet itself there without error.",
            "The mesh fingerprint is a partial one: file names, sizes and a "
            "bounded head of each file, not full content.",
        ],
        **(extra or {}),
    }
