"""Release provenance a public consumer can verify from the published files.

A commit sha in a published certificate is a promise the reader cannot
check: it points into a repository they do not have. So nothing public
cites one. What a certificate records instead is

  * the code version and the frozen policy version + policy hash, and
  * a SOURCE TREE DIGEST -- a content digest over every published file that
    determines a result, computed with the same frozen manifest
    serialisation used for mesh identity (see `windcheck.manifest`).

Both are recomputable from the release tree alone:

    uv run python -m windcheck.provenance

prints the digest and the file count; `--verify <digest>` exits non-zero
when the tree does not match, and `--manifest` prints every row that went
into it.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .manifest import (MANIFEST_SCHEMA, SERIALISATION_RULE, digest, file_row,
                       serialise)

PROVENANCE_SCHEMA = "windcheck_provenance/v1"

REPO_ROOT = Path(__file__).resolve().parents[2]

# Every published path whose bytes can change a result. Deliberately
# excludes tests, notes, docs and generated outputs: a consumer verifying a
# certificate is verifying the code that produced it.
SOURCE_TREE_GLOBS: tuple[str, ...] = (
    "pyproject.toml",
    "uv.lock",
    "src/windcheck/*.py",
    "engines/*.cpp",
    "engines/*.h",
    "bench/*.py",
)

SOURCE_TREE_RULE = (
    "sha256 over the canonical file manifest of the published source tree: "
    f"paths matching {list(SOURCE_TREE_GLOBS)} relative to the repository "
    "root, excluding __pycache__, serialised by the frozen rule below. "
    "Recompute with `uv run python -m windcheck.provenance`; no repository "
    "history is required.")


def source_files(root: str | Path = REPO_ROOT) -> list[str]:
    """Sorted repo-relative posix paths of every source file in the digest."""
    root = Path(root)
    names: set[str] = set()
    for pattern in SOURCE_TREE_GLOBS:
        for p in root.glob(pattern):
            if not p.is_file() or "__pycache__" in p.parts:
                continue
            names.add(p.relative_to(root).as_posix())
    return sorted(names)


def source_manifest(root: str | Path = REPO_ROOT) -> list[dict]:
    root = Path(root)
    return [file_row(root, rel) for rel in source_files(root)]


def source_tree_digest(root: str | Path = REPO_ROOT) -> str:
    return digest(source_manifest(root))


def release_provenance(root: str | Path = REPO_ROOT,
                       policy_version: str | None = None,
                       policy_hash: str | None = None) -> dict:
    """The provenance block every public artifact carries.

    Self-contained by construction: no commit sha, no branch name, nothing
    that points outside the published files.
    """
    rows = source_manifest(root)
    block = {
        "schema": PROVENANCE_SCHEMA,
        "code_version": __version__,
        "source_tree_digest": digest(rows),
        "source_tree_n_files": len(rows),
        "source_tree_globs": list(SOURCE_TREE_GLOBS),
        "manifest_schema": MANIFEST_SCHEMA,
        "serialisation": SERIALISATION_RULE,
        "rule": SOURCE_TREE_RULE,
        "verify_command": "uv run python -m windcheck.provenance",
    }
    if policy_version is None or policy_hash is None:
        try:
            from .excise import FROZEN_POLICY_VERSION, frozen_policy_hash
            policy_version = policy_version or FROZEN_POLICY_VERSION
            policy_hash = policy_hash or frozen_policy_hash()
        except Exception:                                    # pragma: no cover
            pass
    if policy_version is not None:
        block["policy_version"] = policy_version
    if policy_hash is not None:
        block["policy_hash"] = policy_hash
    return block


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m windcheck.provenance",
        description="Recompute the release provenance of this source tree.")
    ap.add_argument("--root", default=str(REPO_ROOT))
    ap.add_argument("--manifest", action="store_true",
                    help="print every file row that feeds the digest")
    ap.add_argument("--serialisation", action="store_true",
                    help="print the exact bytes the digest is taken over")
    ap.add_argument("--json", action="store_true",
                    help="print the whole provenance block as JSON")
    ap.add_argument("--verify", metavar="DIGEST",
                    help="exit non-zero unless the tree hashes to DIGEST")
    a = ap.parse_args(argv)
    root = Path(a.root)

    rows = source_manifest(root)
    d = digest(rows)
    if a.serialisation:
        sys.stdout.write(serialise(rows).decode("utf-8"))
        return 0
    if a.manifest:
        for row in sorted(rows, key=lambda r: r["path"]):
            mark = row["sha256"] if row["present"] else "absent"
            print(f"{mark}  {row['size'] if row['present'] else '-':>9}  "
                  f"{row['path']}")
    if a.json:
        print(json.dumps(release_provenance(root), indent=2, sort_keys=True))
    elif not a.manifest:
        print(f"source_tree_digest {d}")
        print(f"files              {len(rows)}")
        print(f"code_version       {__version__}")
    if a.verify:
        if a.verify.strip().lower() != d:
            print(f"MISMATCH: tree hashes to {d}, expected {a.verify.strip()}",
                  file=sys.stderr)
            return 1
        print("OK: source tree matches the recorded digest")
    return 0


if __name__ == "__main__":                                   # pragma: no cover
    raise SystemExit(main())
