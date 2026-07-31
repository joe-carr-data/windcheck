"""Canonical content manifests for tifxyz meshes, and the one shared answer
to "is this base the original mesh?".

Two things live here because they are the same idea twice.

1.  A mesh is not three coordinate planes. A reader consumes `x.tif`,
    `y.tif`, `z.tif`, an optional `mask.tif` (or `mask.png`) that decides
    which triangles exist at all, and `meta.json`. Hashing only the
    coordinate planes leaves the mask free to move underneath a recorded
    hash — and the mask is exactly the carrier of an excision. So mesh
    identity is a manifest over every file the reader consumes
    semantically, with a single directory digest computed over an
    explicitly defined serialisation of that manifest.

2.  "Is the base the original published mesh?" must be decided by content,
    not by path. Paths change in downloaded archives and fresh working
    directories; the bytes do not. One enum, one predicate, used by the
    pipeline, the corpus base manifest, the release index and the headline
    decision rule alike.

SERIALISATION (frozen, part of the schema):

    a manifest is a list of file rows, each

        {"path": str, "present": bool, "size": int|None, "sha256": str|None}

    sorted by `path` (byte order of the UTF-8 encoding). Each row
    serialises to one line

        present:  f"{path}\\0{size}\\0{sha256}"
        absent:   f"{path}\\0absent\\0absent"

    lines are joined with "\\n" and a trailing "\\n" is appended. The
    directory digest is the sha256 of that blob encoded UTF-8.

A declared file that is not on disk is written into the manifest as an
absent row and serialised. It is never skipped: "no mask" and "a mask I
did not look at" have to hash differently.
"""

from __future__ import annotations

import hashlib
import json
from enum import Enum
from pathlib import Path

MANIFEST_SCHEMA = "windcheck_mesh_manifest/v1"

# Every file a tifxyz reader consumes semantically. `mask.tif`/`mask.png`
# decide which quads survive; `meta.json` carries the grid metadata. Both
# are optional on disk and both are always declared in the manifest.
COORDINATE_PLANES: tuple[str, ...] = ("x.tif", "y.tif", "z.tif")
SEMANTIC_FILES: tuple[str, ...] = COORDINATE_PLANES + (
    "mask.tif", "mask.png", "meta.json")

ABSENT = "absent"

SERIALISATION_RULE = (
    "manifest rows sorted by path; each row serialised UTF-8 as "
    "\"path\\0size\\0sha256\" when the file is present and "
    "\"path\\0absent\\0absent\" when it is declared but absent; rows joined "
    "by \\n with a trailing \\n; digest = sha256 of that blob. Declared "
    "files that are absent are serialised explicitly, never skipped.")

MANIFEST_RULE = (
    "mesh identity is a content manifest over every file a tifxyz reader "
    "consumes semantically -- x/y/z.tif, mask.tif or mask.png, meta.json -- "
    "not the coordinate planes alone, so a mask change (the carrier of an "
    "excision) cannot leave a recorded mesh hash valid.")


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def file_row(root: str | Path, rel: str) -> dict:
    """One manifest row for `root/rel`, present or explicitly absent."""
    p = Path(root) / rel
    if p.is_file():
        return {"path": rel, "present": True, "size": p.stat().st_size,
                "sha256": sha256_file(p)}
    return {"path": rel, "present": False, "size": None, "sha256": None}


def serialise(files: list[dict]) -> bytes:
    """The frozen serialisation of a manifest's file rows."""
    lines = []
    for row in sorted(files, key=lambda r: str(r["path"])):
        if row.get("present") and row.get("sha256"):
            lines.append(f"{row['path']}\0{int(row['size'])}\0{row['sha256']}")
        else:
            lines.append(f"{row['path']}\0{ABSENT}\0{ABSENT}")
    return ("\n".join(lines) + "\n").encode("utf-8")


def digest(files: list[dict]) -> str:
    """sha256 over the canonical serialisation of the file rows."""
    return hashlib.sha256(serialise(files)).hexdigest()


def mesh_manifest(mesh: str | Path,
                  extra: tuple[str, ...] | list[str] = ()) -> dict:
    """The canonical manifest of a tifxyz directory.

    Carries the file rows, the directory digest, and -- as a compatibility
    view derived from the same rows, never independently computed -- the
    per-plane `x`/`y`/`z`/`mask` sha256 spellings older records used.
    """
    mesh = Path(mesh)
    names = list(SEMANTIC_FILES) + [n for n in extra
                                    if n not in SEMANTIC_FILES]
    rows = sorted((file_row(mesh, n) for n in names),
                  key=lambda r: str(r["path"]))
    by_name = {r["path"]: r["sha256"] for r in rows}
    view = {a: by_name.get(f"{a}.tif") for a in ("x", "y", "z")}
    view["mask"] = by_name.get("mask.tif") or by_name.get("mask.png")
    return {"schema": MANIFEST_SCHEMA, "files": rows,
            "digest": digest(rows), "serialisation": SERIALISATION_RULE,
            **view}


# ------------------------------------------------- comparing two records

def _identity_map(record) -> dict[str, str] | None:
    """Normalise any mesh-identity record to {filename: sha256} over files
    it declares present. Accepts a v1 manifest, a legacy `{"x": sha, ...}`
    or `{"x.tif": sha, ...}` plane dict, and returns None when the record
    carries no usable hash at all.
    """
    if not isinstance(record, dict):
        return None
    out: dict[str, str] = {}
    rows = record.get("files")
    if isinstance(rows, list):
        for row in rows:
            if isinstance(row, dict) and row.get("present") and row.get(
                    "sha256"):
                out[str(row["path"])] = str(row["sha256"])
        return out or None
    for key, val in record.items():
        if not isinstance(val, str) or len(val) != 64:
            continue
        name = str(key)
        if name in ("digest", "schema", "serialisation"):
            continue
        if not (name.endswith(".tif") or name.endswith(".png")
                or name.endswith(".json")):
            name = f"{name}.tif"
        out[name] = val
    return out or None


def same_mesh(a, b) -> bool | None:
    """Semantic equality of two mesh-identity records: identical set of
    present files, identical sha256 for each. None when either record is
    unusable, so an undecidable comparison can never read as agreement.
    """
    ma, mb = _identity_map(a), _identity_map(b)
    if ma is None or mb is None:
        return None
    return ma == mb


# ------------------------------------------------------------- base kind

class BaseKind(str, Enum):
    """What the excision pass was handed. The only two answers there are."""

    ORIGINAL = "original"
    DISPLACEMENT_REPAIRED = "displacement_repaired"


# Spellings older records used for BaseKind.ORIGINAL. Accepted on read,
# never written.
_LEGACY_ALIASES = {
    "original_published": BaseKind.ORIGINAL,
    "original_mesh": BaseKind.ORIGINAL,
    "published": BaseKind.ORIGINAL,
    "repaired": BaseKind.DISPLACEMENT_REPAIRED,
    "displacement": BaseKind.DISPLACEMENT_REPAIRED,
}


def normalise_base_kind(value) -> BaseKind | None:
    """Any historical spelling -> the enum. None if unrecognised."""
    if isinstance(value, BaseKind):
        return value
    if not isinstance(value, str):
        return None
    v = value.strip().lower()
    try:
        return BaseKind(v)
    except ValueError:
        return _LEGACY_ALIASES.get(v)


def base_kind_from_manifests(input_record, original_record) -> BaseKind | None:
    """The hash-verified answer: the base is ORIGINAL exactly when the input
    manifest equals the original manifest. None when undecidable.
    """
    eq = same_mesh(input_record, original_record)
    if eq is None:
        return None
    return BaseKind.ORIGINAL if eq else BaseKind.DISPLACEMENT_REPAIRED


def verify_base_kind(declared, input_record, original_record,
                     *, label: str = "") -> dict:
    """One shared verdict on a declared base_kind.

    Never decides by path. Returns
    `{"declared", "verified", "agrees", "decidable", "reason"}`; `agrees`
    is False both when the declaration contradicts the hashes and when the
    hashes cannot decide, so a caller that requires agreement fails loudly
    on missing evidence too.
    """
    dec = normalise_base_kind(declared)
    ver = base_kind_from_manifests(input_record, original_record)
    where = f"{label}: " if label else ""
    if ver is None:
        return {"declared": dec.value if dec else declared, "verified": None,
                "agrees": False, "decidable": False,
                "reason": (f"{where}base kind is not decidable from content: "
                           "the input and original mesh manifests do not both "
                           "carry usable hashes")}
    if dec is None:
        return {"declared": declared, "verified": ver.value, "agrees": False,
                "decidable": True,
                "reason": (f"{where}unrecognised declared base kind "
                           f"{declared!r}; content says {ver.value}")}
    if dec is not ver:
        return {"declared": dec.value, "verified": ver.value, "agrees": False,
                "decidable": True,
                "reason": (f"{where}declared base kind {dec.value} disagrees "
                           f"with the hash-verified answer {ver.value} "
                           "(input manifest vs original manifest)")}
    return {"declared": dec.value, "verified": ver.value, "agrees": True,
            "decidable": True, "reason": ""}


class BaseKindDisagreement(RuntimeError):
    """A declared base kind that the manifests refuse to confirm."""


def require_base_kind(declared, input_record, original_record,
                      *, label: str = "") -> BaseKind:
    """`verify_base_kind`, raising on anything short of agreement."""
    v = verify_base_kind(declared, input_record, original_record, label=label)
    if not v["agrees"]:
        raise BaseKindDisagreement(v["reason"])
    return BaseKind(v["verified"])


def manifest_json(mesh: str | Path) -> str:
    return json.dumps(mesh_manifest(mesh), indent=2, sort_keys=True)
