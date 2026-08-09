"""Tie each census report to the operand bytes it describes.

A census report on its own says "this surface is clean". It does not say
WHICH surface. Two traces with the same grid produce reports that are
interchangeable to any checker that only compares dimensions -- so a
clean report can be paired with a dirty candidate and the candidate is
reported clean.

This writes, per report, the sha256 of every file the loader consumes
from the operand, the sha256 of the report itself, and the identity of
the validator that produced it. `topology_eval.py --census-binding`
refuses a report whose binding does not match the candidate in front of
it.

What this CANNOT do: make a self-reported census trustworthy. Whoever
produced the report can produce the binding. It makes an accidental or
casual swap impossible, and a census this tool did not run is never
marked authoritative regardless.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

LOADED_FILES = ("x.tif", "y.tif", "z.tif", "mask.tif", "meta.json")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def operand_hashes(d: Path) -> dict:
    # Absence is recorded as null, not omitted: "no mask" and "a mask I
    # did not look at" must not be the same binding.
    return {n: (sha256(d / n) if (d / n).is_file() else None)
            for n in LOADED_FILES}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--reports", required=True)
    ap.add_argument("--image", default=None)
    ap.add_argument("--image-digest", default=None)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--villa", default=None)
    a = ap.parse_args()

    validator = {"image": a.image, "image_id": a.image_digest,
                 "binary_in_image": a.bin, "villa_mount": a.villa}
    if a.villa and a.bin:
        # The binary as it exists on the host mount, when reachable.
        host = Path(a.villa) / a.bin.lstrip("/").replace("villa/", "", 1)
        if host.is_file():
            validator["binary_sha256"] = sha256(host)
    if "binary_sha256" not in validator:
        validator["binary_sha256"] = None
        validator["note"] = ("the terminal binary was not hashable from "
                             "the host; the image id is the identity of "
                             "record")

    reports = Path(a.reports)
    written, missing = 0, []
    for line in Path(a.manifest).read_text().splitlines():
        if not line.strip():
            continue
        name, _, sdir = line.partition("\t")
        rep = reports / f"{name}.json"
        d = Path(sdir)
        if not rep.is_file() or not d.is_dir():
            missing.append(name)
            continue
        binding = {
            "schema": "windcheck_census_binding/v1",
            "name": name,
            "report": rep.name,
            "report_sha256": sha256(rep),
            "operand_sha256": operand_hashes(d),
            "validator": validator,
            "what_this_proves": (
                "this report was produced from an operand with exactly "
                "these bytes. It does not make a self-reported census "
                "authoritative."),
        }
        (reports / f"{name}.binding.json").write_text(
            json.dumps(binding, indent=1) + "\n")
        written += 1

    print(f"bound {written} reports"
          + (f"; {len(missing)} could not be bound: {missing[:5]}"
             if missing else ""))
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
