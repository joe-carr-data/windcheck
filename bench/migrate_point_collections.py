"""Wrap already-emitted point collections in the envelope VC3D requires.

Every overlay written before this fix is a bare `{"1": {...}}` map: the
value of what should have been the `collections` field, with no version
key. `PointCollections::loadFromJSON` rejects that outright, so those
files parse as JSON and load as nothing.

This rewrites them in place, idempotently: a file already carrying the
version key is left untouched, and anything that looks like neither form
is reported and skipped rather than guessed at.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from windcheck.certificate import PC_VERSION, PC_VERSION_KEY


def looks_like_bare_collections(doc: object) -> bool:
    """A bare map of numeric-string ids to collection objects."""
    if not isinstance(doc, dict) or not doc:
        return False
    for key, val in doc.items():
        if not key.isdigit() or not isinstance(val, dict):
            return False
        if "points" not in val or "name" not in val:
            return False
    return True


def migrate(path: Path, dry_run: bool) -> str:
    try:
        doc = json.loads(path.read_text())
    except json.JSONDecodeError as e:
        return f"UNPARSEABLE ({e})"

    if isinstance(doc, dict) and PC_VERSION_KEY in doc:
        return "already-current"
    if not looks_like_bare_collections(doc):
        return "UNRECOGNISED shape, skipped"

    if not dry_run:
        path.write_text(json.dumps(
            {PC_VERSION_KEY: PC_VERSION, "collections": doc}, indent=1))
    return "migrated"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("roots", nargs="+", help="directories to walk")
    ap.add_argument("--glob", default="*_points.json")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    counts: dict[str, int] = {}
    problems = []
    for root in a.roots:
        for p in sorted(Path(root).rglob(a.glob)):
            verdict = migrate(p, a.dry_run)
            counts[verdict] = counts.get(verdict, 0) + 1
            if verdict not in ("migrated", "already-current"):
                problems.append(f"{p}: {verdict}")

    for k in sorted(counts):
        print(f"{counts[k]:6d}  {k}")
    for line in problems:
        print(f"  {line}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
