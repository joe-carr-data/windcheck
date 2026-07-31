#!/usr/bin/env python3
"""Make already-published certificates self-contained about provenance.

The corpus excision certificates under `out/excised/corpus/` were written
before release provenance was made verifiable-from-the-release. Each one
cites a commit sha and a branch name from a working repository. A public
consumer holds neither, so those fields are claims they cannot check, and
they point at history that is not part of the release.

This pass removes exactly those unverifiable citations and nothing else:

    code_commit                          (top level)
    code_provenance.commit
    code_provenance.branch
    code_provenance.git_status_porcelain

Every measured field is untouched, and so is every provenance field a
reader CAN verify from the published files -- the driver, excision-module,
engine-binary, engine-source and lockfile sha256s already on the
certificate, and the frozen policy version and hash. A short note records
that the removal happened and why, so the transformation is visible on the
artifact rather than silent.

Nothing is invented: no digest is back-dated onto a certificate produced by
a different tree.

    uv run python bench/publish_provenance.py --check
    uv run python bench/publish_provenance.py --apply
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DIR = REPO_ROOT / "out" / "excised" / "corpus"

TOP_LEVEL_DROP = ("code_commit",)
NESTED_DROP = ("commit", "branch", "git_status_porcelain")

NOTE_KEY = "code_provenance_note"
NOTE = (
    "This certificate was produced before release provenance was made "
    "self-contained. The repository commit sha, branch name and working-tree "
    "status it carried have been removed: a reader of the release cannot "
    "resolve them, so they were claims rather than evidence. Everything a "
    "reader CAN verify is still here -- the sha256 of the driver, the "
    "excision module, the census engine binary and source and the lockfile, "
    "the frozen policy version and hash, and the input and output mesh "
    "hashes. The release-wide source-tree digest is in "
    "out/release/index.json under `provenance`; recompute it with "
    "`uv run python -m windcheck.provenance`. No measured value on this "
    "certificate was changed.")


def strip_doc(doc: dict) -> tuple[dict, list[str]]:
    """Return (document, removed field names). Pure; does not touch disk."""
    removed: list[str] = []
    out = dict(doc)
    for key in TOP_LEVEL_DROP:
        if key in out:
            removed.append(key)
            out.pop(key)
    prov = out.get("code_provenance")
    if isinstance(prov, dict):
        prov = dict(prov)
        for key in NESTED_DROP:
            if key in prov:
                removed.append(f"code_provenance.{key}")
                prov.pop(key)
        out["code_provenance"] = prov
    if removed:
        out[NOTE_KEY] = NOTE
    return out, removed


def process(path: Path, apply: bool) -> list[str]:
    doc = json.loads(path.read_text())
    out, removed = strip_doc(doc)
    if removed and apply:
        path.write_text(json.dumps(out, indent=2) + "\n")
    return removed


def refresh_certificate_hashes(root: Path) -> int:
    """Re-point the records that INDEX the certificates at their new bytes.

    `corpus_summary.jsonl` and `verification.json` each carry the sha256 of
    the certificate file they describe. Removing an unverifiable provenance
    citation changes those bytes, so those two indexes are re-pointed at the
    files as they now stand. Nothing measured is touched, and the base
    manifest's own `base_manifest_sha256` -- which pins WHAT WAS CUT -- is
    left exactly as recorded.
    """
    corpus = root / "out" / "excised" / "corpus"
    changed = 0

    summary = corpus / "corpus_summary.jsonl"
    if summary.is_file():
        lines = []
        for line in summary.read_text().splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            cert = row.get("certificate")
            path = root / cert if cert else None
            if path and path.is_file() and row.get("certificate_sha256"):
                new = sha256_file(path)
                if new != row["certificate_sha256"]:
                    row["certificate_sha256"] = new
                    changed += 1
            lines.append(json.dumps(row))
        summary.write_text("\n".join(lines) + "\n")

    ver = corpus / "verification.json"
    if ver.is_file():
        doc = json.loads(ver.read_text())
        for rec in doc.get("segments", []):
            cert = rec.get("certificate")
            path = root / cert if cert else None
            if path and path.is_file() and rec.get("certificate_sha256"):
                new = sha256_file(path)
                if new != rec["certificate_sha256"]:
                    rec["certificate_sha256"] = new
                    changed += 1
        ver.write_text(json.dumps(doc, indent=2) + "\n")
    return changed


def sha256_file(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(DEFAULT_DIR))
    ap.add_argument("--also", action="append", default=[],
                    help="an extra published JSON document to strip "
                         "(e.g. out/corpus_bases.json)")
    ap.add_argument("--glob", default="*.json")
    ap.add_argument("--apply", action="store_true",
                    help="rewrite the files (default is a dry run)")
    a = ap.parse_args(argv)

    files = sorted(Path(a.dir).glob(a.glob)) + [Path(x) for x in a.also]
    touched = fields = 0
    for f in files:
        removed = process(f, a.apply)
        if removed:
            touched += 1
            fields += len(removed)
    if a.apply:
        n = refresh_certificate_hashes(REPO_ROOT)
        print(f"re-pointed {n} certificate index hash(es)")
    verb = "removed" if a.apply else "would remove"
    print(f"{len(files)} file(s) scanned; {verb} {fields} unverifiable "
          f"provenance field(s) across {touched} file(s)")
    return 0


if __name__ == "__main__":                                   # pragma: no cover
    raise SystemExit(main())
