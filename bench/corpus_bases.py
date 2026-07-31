"""Round-28 BASE MANIFEST for the 185 pinned segments.

DECISIONS.md 2026-07-31 (Round 28) Q3 requires that, before the corpus
excision pass runs, we write down -- per segment -- exactly which mesh is
about to be cut, why that mesh was selected, and what its bytes hash to:

    "Write a BASE MANIFEST first (per segment: selected input path,
     original | displacement_repaired, input hashes, originating
     repair-certificate hash, original published mesh hashes, voxel scale)."

The base is the displacement-repaired mesh when -- and ONLY when -- the mesh
on disk hash-verifies against the OUTPUT hashes recorded in that segment's
repair certificate (Q3: "the 103 hash-verified displacement-repaired meshes
+ originals elsewhere"). An unverifiable repaired mesh is never silently
accepted: the base falls back to the original published mesh and the reason
is written into the record.

Segment enumeration is imported from bench/repair_corpus.py so the pinned
set here is byte-identical to the set the executor shards over. Partial
workdirs (out/repaired/multi/work_*) are ignored entirely, per Q3.

Exact duplicates (same base geometry about to be cut) are grouped by
geometry_key and aliased to the lexicographically smallest segment name, so
the pass can publish both artifact-count and unique-geometry-weighted
summaries.

    uv run python bench/corpus_bases.py [--out out/corpus_bases.json]
                                        [--segments FILE] [--check]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
# repair_segment.py inserts *relative* "src"/"bench" paths, so make the
# import work from any cwd.
for _p in (REPO_ROOT / "src", REPO_ROOT / "bench"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from repair_segment import CORPORA, RES_UM                      # noqa: E402
from windcheck.manifest import BaseKind, mesh_manifest          # noqa: E402
from windcheck.provenance import release_provenance             # noqa: E402

SCHEMA = "corpus_bases/v1"
AXES = ("x", "y", "z")
DEFAULT_VOXEL_UM = 7.91

REPAIRED_MESH_ROOT = REPO_ROOT / "out" / "repaired" / "multi" / "meshes"
CERT_ROOT = REPO_ROOT / "out" / "repaired" / "multi"
DATA_DIR = REPO_ROOT / "data"
DEFAULT_OUT = REPO_ROOT / "out" / "corpus_bases.json"

# Keys excluded from --check drift comparison: wall-clock is not content,
# and the code provenance block moves whenever the source tree moves,
# without the manifest's own content changing.
VOLATILE_TOP_KEYS = ("generated_utc", "code_commit", "provenance")


# ---------------------------------------------------------------- primitives

def sha_file(p: Path) -> str:
    """sha256 of a file's bytes (same discipline as repair_corpus.sha)."""
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def voxel_scale(name: str, default: float = DEFAULT_VOXEL_UM) -> float:
    """Voxel scale in um parsed off a tifxyz mesh name, exactly as
    bench/excise_shadow.py main() does it: RES_UM match else 7.91."""
    m = RES_UM.search(str(name))
    return float(m.group(1)) if m else default


def mesh_hashes(mesh_dir: Path) -> dict[str, str | None]:
    """The canonical content manifest of a tifxyz directory.

    Delegates to `windcheck.manifest.mesh_manifest`, so the base manifest,
    the pipeline certificates and the release index all pin mesh identity
    the same way: every file a reader consumes semantically -- x/y/z.tif,
    mask.tif or mask.png, meta.json -- with a declared-but-absent file
    written in explicitly, plus one directory digest over the frozen
    serialisation. The {x,y,z,mask} keys older records used are still
    present, derived from the same rows.
    """
    return mesh_manifest(mesh_dir)


def geometry_key(hashes: dict[str, str | None]) -> str:
    """sha256 over the ORDERED x/y/z hash triple. Two segments share a
    geometry_key iff the three coordinate planes are byte-identical."""
    triple = [hashes.get(a) for a in AXES]
    return hashlib.sha256(
        json.dumps(triple, separators=(",", ":")).encode()).hexdigest()


def repo_path(p: Path) -> str:
    """Repo-relative posix path when possible (machine-independent), else
    the absolute path."""
    p = Path(p)
    try:
        return Path(p).resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return p.as_posix()


def head_commit(root: Path = REPO_ROOT) -> str:
    r = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                       capture_output=True, text=True)
    return r.stdout.strip()


# ------------------------------------------------------------ base selection

def classify_base(segment: str, original_mesh: Path,
                  original_hashes: dict[str, str | None],
                  repaired_root: Path = REPAIRED_MESH_ROOT,
                  cert_root: Path = CERT_ROOT) -> dict:
    """Decide the mesh the excision pass will actually cut.

    displacement_repaired REQUIRES both (a) out/repaired/multi/meshes/
    <segment>_repaired.tifxyz on disk with all three planes, and (b) those
    planes hashing equal to hashes.output in <segment>_multi_certificate
    .json. Anything else falls back to the original published mesh with
    the reason recorded. Never silently accept an unverified mesh.
    """
    repaired = Path(repaired_root) / f"{segment}_repaired.tifxyz"
    cert = Path(cert_root) / f"{segment}_multi_certificate.json"
    cert_exists = cert.is_file()
    cert_sha = sha_file(cert) if cert_exists else None

    rec = {
        "base_kind": BaseKind.ORIGINAL.value,
        "base_mesh": Path(original_mesh),
        "base_hashes": dict(original_hashes),
        "repair_certificate": cert if cert_exists else None,
        "repair_certificate_sha256": cert_sha,
        "repair_cert_output_hashes_verified": None,
        "repair_verification_note": "",
    }

    have_mesh = repaired.is_dir() and all(
        (repaired / f"{a}.tif").is_file() for a in AXES)

    if not have_mesh:
        if repaired.is_dir():
            rec["repair_verification_note"] = (
                f"{repo_path(repaired)} exists but is missing coordinate "
                f"planes; base is the original published mesh")
            rec["repair_cert_output_hashes_verified"] = False
        elif cert_exists:
            rec["repair_verification_note"] = (
                "no displacement-repaired mesh on disk; a repair certificate "
                f"exists at {repo_path(cert)} but emitted no mesh (segment "
                "was already clean or the repair was refused), so the base "
                "is the original published mesh")
        else:
            rec["repair_verification_note"] = (
                "no displacement-repaired mesh and no repair certificate; "
                "base is the original published mesh")
        return rec

    if not cert_exists:
        rec["repair_cert_output_hashes_verified"] = False
        rec["repair_verification_note"] = (
            f"repaired mesh {repo_path(repaired)} exists but there is no "
            f"repair certificate at {repo_path(cert)}; refusing an "
            "unverified repaired mesh, base falls back to the original")
        return rec

    try:
        doc = json.loads(cert.read_text())
    except (OSError, ValueError) as exc:
        rec["repair_cert_output_hashes_verified"] = False
        rec["repair_verification_note"] = (
            f"repair certificate {repo_path(cert)} is unreadable ({exc}); "
            "refusing an unverified repaired mesh, base falls back to the "
            "original")
        return rec

    out_h = doc.get("hashes", {})
    out_h = out_h.get("output") if isinstance(out_h, dict) else None
    if not isinstance(out_h, dict):
        rec["repair_cert_output_hashes_verified"] = False
        rec["repair_verification_note"] = (
            f"repair certificate {repo_path(cert)} records no "
            "hashes.output block; refusing an unverified repaired mesh, "
            "base falls back to the original")
        return rec

    disk = mesh_hashes(repaired)
    missing = [a for a in AXES if not out_h.get(f"{a}.tif")]
    if missing:
        rec["repair_cert_output_hashes_verified"] = False
        rec["repair_verification_note"] = (
            f"repair certificate {repo_path(cert)} hashes.output is missing "
            f"{','.join(missing)}; refusing an unverified repaired mesh, "
            "base falls back to the original")
        return rec

    bad = [a for a in AXES if disk[a] != out_h[f"{a}.tif"]]
    if bad:
        rec["repair_cert_output_hashes_verified"] = False
        rec["repair_verification_note"] = (
            f"repaired mesh {repo_path(repaired)} disagrees with certificate "
            f"hashes.output on {','.join(bad)}.tif; refusing an unverified "
            "repaired mesh, base falls back to the original")
        return rec

    rec["base_kind"] = BaseKind.DISPLACEMENT_REPAIRED.value
    rec["base_mesh"] = repaired
    rec["base_hashes"] = disk
    rec["repair_cert_output_hashes_verified"] = True
    rec["repair_verification_note"] = (
        "on-disk x/y/z sha256 match hashes.output in "
        f"{repo_path(cert)}")
    return rec


# -------------------------------------------------------- published manifest

def load_manifests(data_dir: Path = DATA_DIR) -> dict[str, tuple[str, dict]]:
    """volume -> (manifest path str, {relative path: sha256}).

    data/MANIFEST*.json are S3 download manifests, one per sample, each
    carrying its own 'volume' and a flat 'files' list of
    {path, s3_key, bytes, sha256} where path is
    <segment>/mesh/<mesh>.tifxyz/<plane>.
    """
    index: dict[str, tuple[str, dict]] = {}
    for p in sorted(Path(data_dir).glob("MANIFEST*.json")):
        try:
            doc = json.loads(p.read_text())
        except (OSError, ValueError):
            continue
        vol = doc.get("volume")
        files = doc.get("files")
        if not vol or not isinstance(files, list):
            continue
        table = {f["path"]: f.get("sha256") for f in files
                 if isinstance(f, dict) and f.get("path")}
        index[str(vol)] = (repo_path(p), table)
    return index


def published_lookup(index: dict, volume: str, segment: str,
                     mesh_dir_name: str):
    """(hashes|None, source|None) for a segment's published mesh planes."""
    ent = index.get(str(volume))
    if ent is None:
        return None, None
    source, table = ent
    prefix = f"{segment}/mesh/{mesh_dir_name}"
    got = {a: table.get(f"{prefix}/{a}.tif") for a in AXES}
    got["mask"] = table.get(f"{prefix}/mask.tif")
    if all(v is None for v in got.values()):
        return None, None
    return got, source


def manifest_agrees(published: dict | None,
                    disk: dict[str, str | None]) -> bool | None:
    """True iff every hash the manifest actually lists equals disk."""
    if published is None:
        return None
    for k, v in published.items():
        if v is None:
            continue
        if disk.get(k) != v:
            return False
    return True


# --------------------------------------------------------------- duplicates

def assign_duplicates(entries: list[dict]) -> list[dict]:
    """Group by geometry_key; canonical member is the lexicographically
    smallest segment name. Mutates and returns entries."""
    groups: dict[str, list[str]] = {}
    for e in entries:
        groups.setdefault(e["geometry_key"], []).append(e["segment"])
    for k in groups:
        groups[k].sort()
    for e in entries:
        members = groups[e["geometry_key"]]
        canon = members[0]
        e["is_canonical"] = e["segment"] == canon
        e["duplicate_of"] = None if e["is_canonical"] else canon
    return entries


def duplicate_groups(entries: list[dict], key: str = "geometry_key"):
    """Sorted [(key, [segment,...])] for every group of size > 1."""
    groups: dict[str, list[str]] = {}
    for e in entries:
        groups.setdefault(e[key], []).append(e["segment"])
    return sorted((k, sorted(v)) for k, v in groups.items() if len(v) > 1)


# ------------------------------------------------------------------- build

def enumerate_pinned() -> list[dict]:
    """The pinned set, straight from bench/repair_corpus.enumerate_segments
    with cwd-independent corpus roots so the tool runs from anywhere."""
    from repair_corpus import enumerate_segments        # noqa: PLC0415
    corpora = [(corpus, str(REPO_ROOT / root), volume, work)
               for corpus, root, volume, work in CORPORA]
    rows = enumerate_segments(corpora)
    for r in rows:
        r["volume"] = next(v for c, _r, v, _w in CORPORA if c == r["corpus"])
    return rows


def build_entry(row: dict, manifests: dict,
                repaired_root: Path = REPAIRED_MESH_ROOT,
                cert_root: Path = CERT_ROOT) -> dict:
    original = Path(row["mesh"])
    orig_h = mesh_hashes(original)
    published, source = published_lookup(
        manifests, row["volume"], row["segment"], original.name)
    base = classify_base(row["segment"], original, orig_h,
                         repaired_root=repaired_root, cert_root=cert_root)
    return {
        "segment": row["segment"],
        "corpus": row["corpus"],
        "volume": row["volume"],
        "voxel_um": voxel_scale(original.name),
        "original_mesh": repo_path(original),
        "original_hashes": orig_h,
        "published_manifest_hashes": published,
        "published_manifest_source": source,
        "published_manifest_agrees_with_disk": manifest_agrees(published,
                                                               orig_h),
        "base_kind": base["base_kind"],
        "base_mesh": repo_path(base["base_mesh"]),
        "base_hashes": base["base_hashes"],
        "repair_certificate": (repo_path(base["repair_certificate"])
                               if base["repair_certificate"] else None),
        "repair_certificate_sha256": base["repair_certificate_sha256"],
        "repair_cert_output_hashes_verified":
            base["repair_cert_output_hashes_verified"],
        "repair_verification_note": base["repair_verification_note"],
        "geometry_key": geometry_key(base["base_hashes"]),
        "original_geometry_key": geometry_key(orig_h),
        "duplicate_of": None,
        "is_canonical": True,
    }


def build(rows: list[dict], manifests: dict | None = None,
          repaired_root: Path = REPAIRED_MESH_ROOT,
          cert_root: Path = CERT_ROOT,
          code_commit: str | None = None) -> dict:
    """The whole manifest document, deterministic apart from generated_utc."""
    manifests = load_manifests() if manifests is None else manifests
    rows = sorted(rows, key=lambda r: r["segment"])
    entries = [build_entry(r, manifests, repaired_root, cert_root)
               for r in rows]
    assign_duplicates(entries)
    entries.sort(key=lambda e: e["segment"])

    n_canon = sum(1 for e in entries if e["is_canonical"])
    n_rep = sum(1 for e in entries
                if e["base_kind"] == BaseKind.DISPLACEMENT_REPAIRED.value)
    return {
        "schema": SCHEMA,
        "generated_utc": datetime.now(timezone.utc)
                         .strftime("%Y-%m-%dT%H:%M:%SZ"),
        # Provenance a public consumer can verify from the published files
        # alone: code version plus a content digest of the source tree. No
        # commit sha -- a reader of the release cannot resolve one.
        "provenance": release_provenance(),
        "n_segments": len(entries),
        "n_canonical": n_canon,
        "n_duplicate_aliases": len(entries) - n_canon,
        "n_displacement_repaired": n_rep,
        "n_original": len(entries) - n_rep,
        "counts_note": (
            "n_segments counts pinned artifacts, one row per enumerated "
            "segment directory. n_canonical counts UNIQUE BASE GEOMETRIES "
            "(distinct geometry_key); n_duplicate_aliases are byte-identical "
            "aliases of a canonical base and must not be double-counted in "
            "unique-geometry-weighted retention. n_displacement_repaired "
            "counts bases that hash-verified against their repair "
            "certificate's hashes.output; every other base is the original "
            "published mesh. An unverifiable repaired mesh is never used "
            "(see repair_verification_note)."),
        "entries": entries,
    }


def summarize(doc: dict) -> str:
    entries = doc["entries"]
    lines = [f"segments enumerated : {doc['n_segments']}",
             f"base displacement_repaired : {doc['n_displacement_repaired']}",
             f"base original              : {doc['n_original']}"]

    rejected = [e for e in entries
                if e["repair_cert_output_hashes_verified"] is False]
    lines.append(f"repaired meshes rejected (hash/cert failure): "
                 f"{len(rejected)}")
    for e in rejected:
        lines.append(f"  REJECT {e['segment']}: "
                     f"{e['repair_verification_note']}")

    dup = duplicate_groups(entries)
    lines.append(f"duplicate base-geometry groups: {len(dup)} "
                 f"({doc['n_duplicate_aliases']} aliases, "
                 f"{doc['n_canonical']} unique geometries)")
    for k, members in dup:
        lines.append(f"  {k[:12]}  canonical={members[0]}  "
                     f"aliases={', '.join(members[1:])}")

    odup = duplicate_groups(entries, key="original_geometry_key")
    lines.append(f"duplicate ORIGINAL-geometry groups: {len(odup)}")
    for k, members in odup:
        lines.append(f"  {k[:12]}  {', '.join(members)}")

    scales: dict[float, int] = {}
    for e in entries:
        scales[e["voxel_um"]] = scales.get(e["voxel_um"], 0) + 1
    lines.append("voxel scales (um): " + ", ".join(
        f"{v}x{n}" for v, n in sorted(scales.items())))

    bad_pub = [e["segment"] for e in entries
               if e["published_manifest_agrees_with_disk"] is False]
    unlisted = [e["segment"] for e in entries
                if e["published_manifest_hashes"] is None]
    lines.append(f"published-manifest disagreements: {len(bad_pub)}"
                 + (f" ({', '.join(bad_pub)})" if bad_pub else ""))
    lines.append(f"segments absent from data/MANIFEST*.json: {len(unlisted)}"
                 + (f" ({', '.join(unlisted[:8])}"
                    + (" ..." if len(unlisted) > 8 else "") + ")"
                    if unlisted else ""))
    return "\n".join(lines)


def stable(doc: dict) -> dict:
    return {k: v for k, v in doc.items() if k not in VOLATILE_TOP_KEYS}


def diff_docs(old: dict, new: dict) -> list[str]:
    """Human-readable drift report over the non-volatile content."""
    o, n = stable(old), stable(new)
    msgs = []
    for k in sorted(set(o) | set(n)):
        if k == "entries":
            continue
        if o.get(k) != n.get(k):
            msgs.append(f"{k}: {o.get(k)!r} -> {n.get(k)!r}")
    oe = {e["segment"]: e for e in o.get("entries", [])}
    ne = {e["segment"]: e for e in n.get("entries", [])}
    for s in sorted(set(oe) - set(ne)):
        msgs.append(f"entry removed: {s}")
    for s in sorted(set(ne) - set(oe)):
        msgs.append(f"entry added: {s}")
    for s in sorted(set(oe) & set(ne)):
        if oe[s] != ne[s]:
            for f in sorted(set(oe[s]) | set(ne[s])):
                if oe[s].get(f) != ne[s].get(f):
                    msgs.append(f"entry {s}.{f}: "
                                f"{oe[s].get(f)!r} -> {ne[s].get(f)!r}")
    return msgs


def dumps(doc: dict) -> str:
    return json.dumps(doc, indent=1) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--segments", default=None,
                    help="file of segment names, one per line, to restrict to")
    ap.add_argument("--check", action="store_true",
                    help="rebuild and diff against --out; exit 1 on drift, "
                         "write nothing")
    args = ap.parse_args(argv)

    rows = enumerate_pinned()
    if args.segments:
        want = {ln.strip() for ln in Path(args.segments).read_text()
                .splitlines() if ln.strip() and not ln.startswith("#")}
        rows = [r for r in rows if r["segment"] in want]
        missing = want - {r["segment"] for r in rows}
        if missing:
            print(f"WARNING: {len(missing)} requested segment(s) not "
                  f"enumerable: {sorted(missing)}", file=sys.stderr)

    doc = build(rows)
    print(summarize(doc))

    out = Path(args.out)
    if not out.is_absolute():
        out = REPO_ROOT / out

    if args.check:
        if not out.exists():
            print(f"CHECK FAIL: {repo_path(out)} does not exist")
            return 1
        old = json.loads(out.read_text())
        msgs = diff_docs(old, doc)
        od = (old.get("provenance") or {}).get("source_tree_digest")
        nd = (doc.get("provenance") or {}).get("source_tree_digest")
        if od != nd:
            print(f"note: source_tree_digest {od} -> {nd} (not drift)")
        if msgs:
            print(f"CHECK FAIL: {len(msgs)} difference(s) vs "
                  f"{repo_path(out)}")
            for m in msgs[:60]:
                print("  " + m)
            if len(msgs) > 60:
                print(f"  ... {len(msgs) - 60} more")
            return 1
        print(f"CHECK OK: {repo_path(out)} matches a fresh rebuild")
        return 0

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(dumps(doc))
    print(f"wrote {repo_path(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
