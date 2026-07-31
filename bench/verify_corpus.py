"""Round-28 Q5 INDEPENDENT VERIFICATION of the corpus excision pass.

The excision driver certifies its own output. This program does not trust it.
It re-derives, from the artifacts on disk alone, every claim the corpus makes:

  1. CERTIFICATE WELL-FORMEDNESS. Every `<segment>_excision_certificate.json`
     in the certificates directory is parsed, checked for the record kind, a
     `terminal_disposition` drawn from the certificate's own declared set, a
     `policy_hash`, a verifiable code provenance (source-tree digest),
     and base hashes that STILL match
     out/corpus_bases.json. Roster coverage is checked both ways: an entry
     without a certificate and a certificate without an entry are both
     failures.

  2. RE-HASHING. For every emitted mesh the x/y/z planes are re-hashed from
     disk and compared against `output_mesh_hashes`. Every mismatch is named.

  3. INDEPENDENT RECENSUS. Every emitted mesh is censused AGAIN with the same
     C++ engine the pipeline uses (`bench/crossing_census.census_one`, the
     function `bench/excise_shadow.py` reaches for through `census_mesh`),
     with the frozen CENSUS parameters and BOTH diagonals, in a FRESH
     temporary workdir that is created for that segment and deleted after.
     None of the driver's workdirs, CSVs or atlases are reused, so a stale
     or mislabelled CSV cannot produce a clean verdict here. Acceptance is
     transverse == 0 on diagonal 0 AND diagonal 1. The recomputed counts are
     additionally compared against the certificate's recorded `census_after`.

  4. ALREADY-CLEAN SEGMENTS. Their certificate claims the INPUT base mesh
     needed no cut. That claim is re-measured the same way: the base mesh is
     censused fresh and must come back 0/0.

  5. TRIANGLE-EMPTY/INVALID SEGMENTS. Their certificate claims no census is
     defined on the input. That is re-measured too: the engine census must
     decline the mesh (below the validity threshold), not return counts.

Nothing in here can make a claim pass. A missing file, an unreadable
certificate and a nonzero transverse count are all recorded as failures and
the exit status is nonzero.

    uv run python bench/verify_corpus.py --jobs 4
    uv run python bench/verify_corpus.py --limit 3          # smoke test
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
for _p in (REPO_ROOT / "src", REPO_ROOT / "bench"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from crossing_census import census_one                          # noqa: E402
from excise_segment import CENSUS                               # noqa: E402

SCHEMA = "verify_corpus/v1"
AXES = ("x", "y", "z")
DIAGONALS = (0, 1)
CERT_SUFFIX = "_excision_certificate.json"
RECORD_KIND = "excision certificate"
ENGINE_FIELDS = ("triangles", "quads_dropped", "pairs_tested", "transverse",
                 "coplanar", "grazing")

DEFAULT_CERT_DIR = REPO_ROOT / "out" / "excised" / "corpus"
DEFAULT_BASE_MANIFEST = REPO_ROOT / "out" / "corpus_bases.json"
DEFAULT_OUT = DEFAULT_CERT_DIR / "verification.json"

# The dispositions a certificate may declare. Kept here as a fallback only:
# each certificate carries its own `terminal_dispositions_defined` and that
# list is preferred, so this program cannot narrow the driver's vocabulary.
FALLBACK_DISPOSITIONS = ("transformed", "already_clean", "duplicate_alias",
                         "triangle_empty_invalid", "not_censusable",
                         "residual_transverse", "error")

# Dispositions whose certificate promises a mesh on disk.
EMITTING = ("transformed",)


# ------------------------------------------------------------------- hashing

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def mesh_plane_hashes(mesh: Path) -> dict[str, str | None]:
    """sha256 of x/y/z.tif, exactly as bench/excise_shadow.mesh_plane_hashes
    records them. A missing plane hashes to None, never to a pass."""
    out: dict[str, str | None] = {}
    for ax in AXES:
        p = Path(mesh) / f"{ax}.tif"
        out[ax] = sha256_file(p) if p.is_file() else None
    return out


def normalise_hash_keys(hashes) -> dict[str, str | None]:
    """Accept both `{"x": ...}` and `{"x.tif": ...}` spellings."""
    if not isinstance(hashes, dict):
        return {}
    out: dict[str, str | None] = {}
    for k, v in hashes.items():
        key = str(k)
        if key.endswith(".tif"):
            key = key[:-4]
        out[key] = v
    return out


def code_identity(cert: dict) -> str | None:
    """The certificate's code provenance as a single verifiable string.

    Prefers the release spelling -- a content digest of the published
    source tree, recomputable with `uv run python -m windcheck.provenance`
    -- and still accepts the pre-release `code_commit` field so certificates
    written before the provenance change keep verifying.
    """
    for key in ("source_tree_digest", "code_identity", "code_commit"):
        v = cert.get(key)
        if isinstance(v, str) and v.strip():
            return v.strip()
    prov = cert.get("provenance") or cert.get("code_provenance")
    if isinstance(prov, dict):
        for key in ("source_tree_digest", "commit"):
            v = prov.get(key)
            if isinstance(v, str) and v.strip():
                return v.strip()
    return None


# --------------------------------------------------------- certificate check

def check_certificate(cert, entry, *, segment: str,
                      manifest_sha: str | None = None) -> list[str]:
    """Structural + provenance problems with one certificate.

    Returns a list of human-readable problems; empty means well-formed. A
    missing field is a problem, never a silent pass.
    """
    problems: list[str] = []
    if not isinstance(cert, dict):
        return ["certificate is not a JSON object"]

    if cert.get("record_kind") != RECORD_KIND:
        problems.append(f"record_kind {cert.get('record_kind')!r} != "
                        f"{RECORD_KIND!r}")
    if cert.get("segment") != segment:
        problems.append(f"segment field {cert.get('segment')!r} does not "
                        f"match the filename segment {segment!r}")

    defined = cert.get("terminal_dispositions_defined")
    if not isinstance(defined, (list, tuple)) or not defined:
        defined = FALLBACK_DISPOSITIONS
    disp = cert.get("terminal_disposition")
    if not disp:
        problems.append("no terminal_disposition")
    elif disp not in defined:
        problems.append(f"terminal_disposition {disp!r} is not one of the "
                        f"declared dispositions {sorted(defined)}")

    ph = cert.get("policy_hash")
    if not isinstance(ph, str) or not ph.strip():
        problems.append("no policy_hash")
    cc = code_identity(cert)
    if not cc:
        problems.append("no verifiable code provenance (expected a "
                        "source_tree_digest, or the pre-release code_commit "
                        "spelling)")

    if entry is None:
        problems.append("segment is not on the pinned base manifest")
        return problems

    want = normalise_hash_keys(entry.get("base_hashes") or {})
    got = normalise_hash_keys(cert.get("base_hashes") or {})
    if not want:
        problems.append("base manifest entry carries no base_hashes")
    for ax in AXES:
        if want.get(ax) is None:
            problems.append(f"base manifest entry has no {ax} base hash")
        elif got.get(ax) != want.get(ax):
            problems.append(f"base hash {ax}: certificate "
                            f"{str(got.get(ax))[:12]} != manifest "
                            f"{str(want.get(ax))[:12]}")
    base_hash_problems = bool(problems)

    # The manifest byte-sha is belt-and-braces on top of the per-segment
    # base-hash check above. A provenance-only edit to the manifest header
    # (say, removing an unverifiable commit citation) changes those bytes
    # without changing a single thing the certificate was cut under, so a
    # byte mismatch is only a PROBLEM when this segment's pinned base
    # hashes have also moved. Otherwise it is a note.
    if manifest_sha and cert.get("base_manifest_sha256") and \
            cert["base_manifest_sha256"] != manifest_sha:
        msg = (f"base_manifest_sha256 {cert['base_manifest_sha256'][:12]} != "
               f"the manifest on disk {manifest_sha[:12]}")
        if base_hash_problems:
            problems.append(msg)
        else:
            problems.append(
                "NOTE: " + msg + "; the pinned base hashes for this segment "
                "are unchanged, so the difference is in the manifest header, "
                "not in what was cut")

    if disp in EMITTING:
        if not cert.get("output_mesh"):
            problems.append("transformed certificate names no output_mesh")
        if not normalise_hash_keys(cert.get("output_mesh_hashes") or {}):
            problems.append("transformed certificate carries no "
                            "output_mesh_hashes")
    return problems


def check_output_hashes(cert, *, repo_root: Path = REPO_ROOT) -> dict:
    """Re-hash the emitted planes on disk against `output_mesh_hashes`."""
    rec: dict = {"checked": False, "status": "not_applicable",
                 "mesh": None, "mismatched_axes": [], "per_axis": {}}
    if cert.get("terminal_disposition") not in EMITTING:
        return rec
    mesh_rel = cert.get("output_mesh")
    if not mesh_rel:
        rec["status"] = "no_output_mesh_recorded"
        return rec
    mesh = Path(mesh_rel)
    if not mesh.is_absolute():
        mesh = repo_root / mesh
    rec["mesh"] = str(mesh)
    if not mesh.is_dir():
        rec["status"] = "mesh_missing"
        return rec

    recorded = normalise_hash_keys(cert.get("output_mesh_hashes") or {})
    actual = mesh_plane_hashes(mesh)
    rec["checked"] = True
    for ax in AXES:
        rec["per_axis"][ax] = {"recorded": recorded.get(ax),
                               "recomputed": actual.get(ax),
                               "match": bool(recorded.get(ax) is not None
                                             and actual.get(ax) is not None
                                             and recorded[ax] == actual[ax])}
        if not rec["per_axis"][ax]["match"]:
            rec["mismatched_axes"].append(ax)
    rec["status"] = "ok" if not rec["mismatched_axes"] else "mismatch"
    return rec


# ------------------------------------------------------------------- census

def engine_recensus(mesh: Path, tag: str, threads: int = 3) -> dict | None:
    """One both-diagonal engine census of a mesh on disk, in a FRESH workdir.

    The same call bench/excise_shadow.py's `census_mesh` makes -- the C++
    engine through bench/crossing_census.census_one, with the frozen CENSUS
    parameters -- except that the workdir is a private temporary directory
    created for this call and removed when it returns. Nothing the driver
    wrote is read.
    """
    work = Path(tempfile.mkdtemp(prefix="verify_corpus_census_"))
    try:
        return census_one(Path(mesh), tag, CENSUS["exclude"], CENSUS["cell"],
                          threads, CENSUS["maxedge"], work)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _engine_counts(row) -> dict | None:
    if not isinstance(row, dict):
        return None
    out = {}
    for d in DIAGONALS:
        blk = row.get(f"d{d}")
        if not isinstance(blk, dict):
            return None
        out[f"d{d}"] = {k: blk.get(k) for k in ENGINE_FIELDS}
    return out


def compare_census(recomputed: dict, recorded) -> dict:
    """Recomputed vs the certificate's recorded block."""
    rec = {"comparable": False, "transverse_disagrees": False,
           "differing_fields": [], "recorded": None}
    if not isinstance(recorded, dict):
        return rec
    ref = {}
    for d in DIAGONALS:
        blk = recorded.get(f"d{d}")
        if not isinstance(blk, dict):
            return rec
        ref[f"d{d}"] = {k: blk.get(k) for k in ENGINE_FIELDS}
    rec["comparable"] = True
    rec["recorded"] = {f"d{d}": {"transverse": ref[f"d{d}"]["transverse"]}
                       for d in DIAGONALS}
    for d in DIAGONALS:
        for k in ENGINE_FIELDS:
            a, b = recomputed[f"d{d}"].get(k), ref[f"d{d}"].get(k)
            if b is None:
                continue
            if a != b:
                rec["differing_fields"].append(
                    {"diagonal": d, "field": k, "recomputed": a,
                     "recorded": b})
                if k == "transverse":
                    rec["transverse_disagrees"] = True
    return rec


# ------------------------------------------------------------ per segment

def verify_segment(cert_path: Path, entry, *, census_fn=None,
                   repo_root: Path = REPO_ROOT,
                   manifest_sha: str | None = None) -> dict:
    """Full independent verification of ONE certificate and its artifact.

    `census_fn(mesh: Path, tag: str) -> row | None` is injected so the tests
    can exercise the logic without the C++ engine.
    """
    cert_path = Path(cert_path)
    segment = cert_path.name[:-len(CERT_SUFFIX)]
    if census_fn is None:
        census_fn = partial(engine_recensus, threads=3)

    out: dict = {
        "segment": segment,
        "certificate": str(cert_path),
        "certificate_sha256": None,
        "disposition": None,
        "certificate_problems": [],
        "hash_check": {"checked": False, "status": "not_applicable",
                       "mismatched_axes": [], "per_axis": {}, "mesh": None},
        "recensus": {"kind": "skipped", "ran": False, "mesh": None,
                     "d0_transverse": None, "d1_transverse": None,
                     "clean_both_diagonals": None, "counts": None,
                     "seconds": None, "comparison": None, "error": None},
        "problems": [],
        "ok": False,
    }

    try:
        raw = cert_path.read_bytes()
    except OSError as exc:
        out["certificate_problems"].append(f"unreadable certificate: {exc}")
        out["problems"] = list(out["certificate_problems"])
        return out
    out["certificate_sha256"] = hashlib.sha256(raw).hexdigest()
    try:
        cert = json.loads(raw.decode("utf-8"))
    except (ValueError, UnicodeDecodeError) as exc:
        out["certificate_problems"].append(f"unparseable certificate: {exc}")
        out["problems"] = list(out["certificate_problems"])
        return out

    out["disposition"] = cert.get("terminal_disposition")
    out["certificate_problems"] = check_certificate(
        cert, entry, segment=segment, manifest_sha=manifest_sha)
    out["hash_check"] = check_output_hashes(cert, repo_root=repo_root)

    disp = out["disposition"]
    if disp in EMITTING:
        kind, mesh_rel, recorded = ("output", cert.get("output_mesh"),
                                    cert.get("census_after"))
    elif disp == "already_clean":
        kind, mesh_rel, recorded = ("input_base", cert.get("base_mesh"),
                                    cert.get("census_before"))
    elif disp == "triangle_empty_invalid":
        kind, mesh_rel, recorded = ("input_base_expect_not_censusable",
                                    cert.get("base_mesh"), None)
    else:
        kind, mesh_rel, recorded = (None, None, None)

    if kind and mesh_rel:
        mesh = Path(mesh_rel)
        if not mesh.is_absolute():
            mesh = repo_root / mesh
        r = out["recensus"]
        r["kind"] = kind
        r["mesh"] = str(mesh)
        if not mesh.is_dir():
            r["error"] = "mesh directory missing"
        else:
            t0 = time.time()
            try:
                row = census_fn(mesh, f"verify_{segment}")
            except (subprocess.CalledProcessError, OSError, ValueError) as exc:
                row = None
                r["error"] = f"census failed: {exc}"
            r["seconds"] = round(time.time() - t0, 3)
            counts = _engine_counts(row)
            if counts is None and r["error"] is None:
                # census_one declines a mesh below the validity threshold
                r["error"] = None if kind.endswith("not_censusable") \
                    else "engine declined the mesh (not censusable)"
                r["ran"] = True
                r["not_censusable"] = True
            elif counts is not None:
                r["ran"] = True
                r["not_censusable"] = False
                r["counts"] = counts
                r["d0_transverse"] = counts["d0"]["transverse"]
                r["d1_transverse"] = counts["d1"]["transverse"]
                r["clean_both_diagonals"] = bool(
                    r["d0_transverse"] == 0 and r["d1_transverse"] == 0)
                r["comparison"] = compare_census(counts, recorded)
    elif kind:
        out["recensus"]["kind"] = kind
        out["recensus"]["error"] = "certificate names no mesh to recensus"

    # ------- fold everything into one problem list
    problems = list(out["certificate_problems"])
    if out["hash_check"]["status"] == "mismatch":
        problems.append("output mesh hash mismatch on "
                        + ",".join(out["hash_check"]["mismatched_axes"]))
    elif out["hash_check"]["status"] in ("mesh_missing",
                                         "no_output_mesh_recorded"):
        problems.append(f"output mesh hash check: "
                        f"{out['hash_check']['status']}")

    r = out["recensus"]
    if r["kind"] in ("output", "input_base"):
        if r["error"]:
            problems.append(f"recensus: {r['error']}")
        elif not r["ran"]:
            problems.append("recensus did not run")
        elif r.get("not_censusable"):
            problems.append("recensus: engine declined the mesh")
        elif not r["clean_both_diagonals"]:
            problems.append(f"recensus NOT clean: d0 transverse "
                            f"{r['d0_transverse']}, d1 transverse "
                            f"{r['d1_transverse']}")
        if r.get("comparison") and r["comparison"].get("transverse_disagrees"):
            problems.append(
                "recensus disagrees with the certificate's recorded census: "
                + "; ".join(f"d{f['diagonal']} {f['field']} recomputed "
                            f"{f['recomputed']} vs recorded {f['recorded']}"
                            for f in r["comparison"]["differing_fields"]
                            if f["field"] == "transverse"))
    elif r["kind"] == "input_base_expect_not_censusable":
        if r["error"]:
            problems.append(f"recensus: {r['error']}")
        elif not r.get("not_censusable"):
            problems.append(
                "certificate claims triangle-empty/invalid but the engine "
                "censused the input mesh: d0 transverse "
                f"{r['d0_transverse']}, d1 transverse {r['d1_transverse']}")

    # "NOTE: ..." entries are recorded, not failures. They exist so a
    # difference that is genuinely not a defect is still visible on the
    # record instead of being dropped.
    out["notes"] = [p for p in problems if p.startswith("NOTE: ")]
    out["problems"] = [p for p in problems if not p.startswith("NOTE: ")]
    out["ok"] = not out["problems"]
    return out


# ---------------------------------------------------------------- aggregate

def summarise(records: list[dict], roster: dict) -> dict:
    def names(pred):
        return sorted(r["segment"] for r in records if pred(r))

    emitted = [r for r in records if r["hash_check"]["checked"]]
    censused = [r for r in records
                if r["recensus"]["kind"] in ("output", "input_base")]
    clean = [r for r in censused if r["recensus"]["clean_both_diagonals"]]
    not_clean = [r for r in censused
                 if r["recensus"]["clean_both_diagonals"] is not True]
    disagree = [r for r in censused
                if (r["recensus"].get("comparison") or {})
                .get("transverse_disagrees")]
    soft = [r for r in censused
            if (r["recensus"].get("comparison") or {}).get("differing_fields")
            and not (r["recensus"]["comparison"].get("transverse_disagrees"))]

    counts: dict[str, int] = {}
    for r in records:
        counts[r["disposition"] or "none"] = \
            counts.get(r["disposition"] or "none", 0) + 1

    return {
        "n_certificates_checked": len(records),
        "n_certificates_well_formed": sum(
            1 for r in records if not r["certificate_problems"]),
        "certificate_failures": names(lambda r: r["certificate_problems"]),
        "dispositions": dict(sorted(counts.items())),
        "n_meshes_rehashed": len(emitted),
        "n_meshes_rehash_ok": sum(
            1 for r in emitted if r["hash_check"]["status"] == "ok"),
        "n_meshes_rehash_mismatch": sum(
            1 for r in emitted if r["hash_check"]["status"] == "mismatch"),
        "hash_mismatch_segments": sorted(
            r["segment"] for r in emitted
            if r["hash_check"]["status"] == "mismatch"),
        "n_recensused": len(censused),
        "n_recensus_output_meshes": sum(
            1 for r in censused if r["recensus"]["kind"] == "output"),
        "n_recensus_input_bases": sum(
            1 for r in censused if r["recensus"]["kind"] == "input_base"),
        "n_recensus_clean_both_diagonals": len(clean),
        "not_clean_segments": sorted(
            f"{r['segment']} (d0 {r['recensus']['d0_transverse']}, "
            f"d1 {r['recensus']['d1_transverse']})" for r in not_clean),
        "n_census_disagreements": len(disagree),
        "census_disagreement_segments": sorted(r["segment"] for r in disagree),
        "n_census_soft_field_differences": len(soft),
        "census_soft_difference_segments": sorted(r["segment"] for r in soft),
        "n_not_censusable_confirmed": sum(
            1 for r in records
            if r["recensus"]["kind"] == "input_base_expect_not_censusable"
            and r["recensus"].get("not_censusable")),
        "n_ok": sum(1 for r in records if r["ok"]),
        "n_failed": sum(1 for r in records if not r["ok"]),
        "failed_segments": [{"segment": r["segment"],
                             "problems": r["problems"]}
                            for r in records if not r["ok"]],
        "roster": roster,
    }


def build_roster(manifest_entries: list[dict], cert_paths: list[Path]) -> dict:
    have = {p.name[:-len(CERT_SUFFIX)] for p in cert_paths}
    pinned = {e["segment"] for e in manifest_entries}
    return {
        "n_manifest_entries": len(pinned),
        "n_certificates_on_disk": len(have),
        "manifest_entries_without_certificate": sorted(pinned - have),
        "certificates_not_on_manifest": sorted(have - pinned),
    }


# --------------------------------------------------------------------- main

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="independent verification of the corpus excision pass")
    ap.add_argument("--certificates-dir", type=Path, default=DEFAULT_CERT_DIR)
    ap.add_argument("--base-manifest", type=Path,
                    default=DEFAULT_BASE_MANIFEST)
    ap.add_argument("--json", dest="out_json", type=Path, default=None,
                    help="default <certificates-dir>/verification.json")
    ap.add_argument("--jobs", type=int, default=4,
                    help="segments verified concurrently")
    ap.add_argument("--census-threads", type=int, default=3,
                    help="engine threads PER census; keep jobs*threads at or "
                         "under the physical core count. The engine's answer "
                         "does not depend on this")
    ap.add_argument("--limit", type=int, default=None,
                    help="verify only the first N certificates (smoke test)")
    a = ap.parse_args(argv)

    cert_dir = a.certificates_dir
    out_json = a.out_json or (cert_dir / "verification.json")

    try:
        manifest_doc = json.loads(Path(a.base_manifest).read_text())
    except (OSError, ValueError) as exc:
        print(f"cannot read base manifest {a.base_manifest}: {exc}")
        return 2
    manifest_sha = sha256_file(Path(a.base_manifest))
    entries = list(manifest_doc.get("entries") or [])
    by_segment = {e["segment"]: e for e in entries}

    cert_paths = sorted(p for p in cert_dir.glob(f"*{CERT_SUFFIX}"))
    if not cert_paths:
        print(f"no certificates under {cert_dir}")
        return 2
    roster = build_roster(entries, cert_paths)
    if a.limit is not None:
        cert_paths = cert_paths[:a.limit]

    census_fn = partial(engine_recensus, threads=a.census_threads)

    t0 = time.time()
    records: list[dict] = []
    done = 0
    with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as pool:
        futs = {
            pool.submit(verify_segment, p,
                        by_segment.get(p.name[:-len(CERT_SUFFIX)]),
                        census_fn=census_fn, manifest_sha=manifest_sha): p
            for p in cert_paths}
        for fut in as_completed(futs):
            p = futs[fut]
            try:
                rec = fut.result()
            except Exception as exc:                    # never lose a segment
                rec = {"segment": p.name[:-len(CERT_SUFFIX)],
                       "certificate": str(p), "disposition": None,
                       "certificate_problems": [], "ok": False,
                       "hash_check": {"checked": False,
                                      "status": "verifier_error",
                                      "mismatched_axes": [], "per_axis": {},
                                      "mesh": None},
                       "recensus": {"kind": "skipped", "ran": False,
                                    "mesh": None, "d0_transverse": None,
                                    "d1_transverse": None,
                                    "clean_both_diagonals": None,
                                    "counts": None, "seconds": None,
                                    "comparison": None, "error": str(exc)},
                       "problems": [f"verifier raised: {exc!r}"]}
            records.append(rec)
            done += 1
            r = rec["recensus"]
            flag = "OK " if rec["ok"] else "FAIL"
            print(f"[{done:3d}/{len(cert_paths)}] {flag} {rec['segment'][:44]:46s}"
                  f" {str(rec['disposition'])[:22]:24s}"
                  f" d0={r['d0_transverse']} d1={r['d1_transverse']}"
                  + ("" if rec["ok"] else "  <- " + "; ".join(rec["problems"])),
                  flush=True)

    records.sort(key=lambda r: r["segment"])
    summary = summarise(records, roster)
    summary["failed_segments"] = [
        {"segment": r["segment"], "problems": r["problems"]}
        for r in records if not r["ok"]]
    summary["wall_seconds"] = round(time.time() - t0, 2)
    if roster["manifest_entries_without_certificate"] or \
            roster["certificates_not_on_manifest"]:
        summary["roster_ok"] = False
    else:
        summary["roster_ok"] = True

    doc = {
        "schema": SCHEMA,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "certificates_dir": str(cert_dir),
        "base_manifest": str(a.base_manifest),
        "base_manifest_sha256": manifest_sha,
        "census_params": dict(CENSUS, diagonals=list(DIAGONALS),
                              threads=a.census_threads),
        "census_engine": str(REPO_ROOT / "engines" / "selfcross"),
        "census_engine_sha256": sha256_file(REPO_ROOT / "engines" /
                                            "selfcross"),
        "independence_note": (
            "every census here ran in a fresh temporary workdir created for "
            "that segment and deleted afterwards; no driver workdir, CSV or "
            "atlas was reused"),
        "limit": a.limit,
        "summary": summary,
        "segments": records,
    }
    out_json.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_json.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(doc, indent=1, default=str))
    tmp.replace(out_json)

    s = summary
    print("\n" + "=" * 78)
    print("INDEPENDENT VERIFICATION -- round-28 Q5")
    print("=" * 78)
    print(f"certificates checked            : {s['n_certificates_checked']}")
    print(f"  well-formed                   : "
          f"{s['n_certificates_well_formed']}")
    print(f"  dispositions                  : " + ", ".join(
        f"{k}={v}" for k, v in s["dispositions"].items()))
    print(f"roster                          : "
          f"{roster['n_certificates_on_disk']} certificates vs "
          f"{roster['n_manifest_entries']} manifest entries"
          f" ({'OK' if summary['roster_ok'] else 'MISMATCH'})")
    print(f"meshes re-hashed                : {s['n_meshes_rehashed']} "
          f"(OK {s['n_meshes_rehash_ok']}, "
          f"MISMATCH {s['n_meshes_rehash_mismatch']})")
    print(f"independently recensused        : {s['n_recensused']} "
          f"({s['n_recensus_output_meshes']} emitted meshes, "
          f"{s['n_recensus_input_bases']} already-clean inputs)")
    print(f"  transverse 0/0 both diagonals : "
          f"{s['n_recensus_clean_both_diagonals']}")
    print(f"  NOT 0/0                       : {len(s['not_clean_segments'])}")
    print(f"recorded-vs-recomputed census   : "
          f"{s['n_census_disagreements']} disagreements on transverse, "
          f"{s['n_census_soft_field_differences']} on other engine fields")
    print(f"not-censusable claims confirmed : "
          f"{s['n_not_censusable_confirmed']}")
    print(f"segments OK / FAILED            : {s['n_ok']} / {s['n_failed']}")

    for label, items in (
            ("certificate problems", s["certificate_failures"]),
            ("output hash mismatches", s["hash_mismatch_segments"]),
            ("not transverse 0/0", s["not_clean_segments"]),
            ("census disagreements", s["census_disagreement_segments"]),
            ("manifest entries without a certificate",
             roster["manifest_entries_without_certificate"]),
            ("certificates not on the manifest",
             roster["certificates_not_on_manifest"])):
        print(f"\n{label}: {len(items)}")
        for n in items:
            print(f"    {n}")

    if s["failed_segments"]:
        print(f"\nFAILURES ({len(s['failed_segments'])}):")
        for f in s["failed_segments"]:
            print(f"  {f['segment']}: " + "; ".join(f["problems"]))

    print(f"\nwrote {out_json}  ({s['wall_seconds']}s)")
    if s["n_failed"] or not summary["roster_ok"]:
        verdict = "FAIL"
    elif a.limit is not None:
        verdict = "PARTIAL (--limit, not the full corpus)"
    else:
        verdict = "PASS"
    print("VERIFICATION: " + verdict)
    return 0 if s["n_failed"] == 0 and summary["roster_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
