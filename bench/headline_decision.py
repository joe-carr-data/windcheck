"""Round-28 Q4 HEADLINE DECISION RULE, evaluated automatically from the index.

DECISIONS.md 2026-07-31 (Round 28) Q4 froze eight conditions and the two
candidate headline sentences BEFORE the corpus pass ran.  Q5 then required
that the decision be computed "AUTOMATICALLY from the index (no editable
prose arithmetic)".  This module is that computation: it reads the base
manifest, the corpus summary and every excision certificate, evaluates the
eight frozen conditions, and PRINTS the sentence with X and Y filled in from
the data.  Nothing here may be edited to make a number nicer; the only way to
change the sentence is to change the corpus.

    uv run python bench/headline_decision.py
        [--base-manifest out/corpus_bases.json]
        [--summary out/excised/corpus/corpus_summary.jsonl]
        [--certificates-dir out/excised/corpus]
        [--expected-total 185] [--json OUT] [--no-strict-artifacts]
        [--area-cache out/headline_original_areas.json] [--no-area-recompute]

Exit status is the gate: 0 == STRONG (all eight conditions pass), 1 ==
FALLBACK, 2 == the inputs could not be read at all.

TWO RETENTION FIGURES (Q3: "publish artifact-count AND unique-geometry-
weighted summaries"):

  * UNIQUE-GEOMETRY-WEIGHTED (this is X, and this is what condition 3 gates)

        1 - sum_canonical A_removed_priced_on_original
            / sum_canonical A_original_canonical

    Duplicate aliases contribute NOTHING to either sum -- a geometry is
    counted exactly once, through its canonical segment.  An already-clean
    canonical segment contributes its full original area to the denominator
    and zero to the numerator.

    An already-clean certificate records no cut, so it carries no
    `headline_area` block.  That is a BOOKKEEPING gap, not a missing
    measurement: the segment was not modified, so nothing was removed and
    its retained fraction is exactly 1.0.  This module therefore PRICES
    those segments itself, using the same definition the cut segments use
    (see `original_canonical_area`), so that every pinned entry appears in
    both retention figures instead of silently dropping out of the
    denominator.

  * ARTIFACT-COUNT-WEIGHTED: the same ratio over all pinned artifacts, with
    each alias priced on the geometry it aliases.  Reported, never gated.

Y is the MINIMUM per-segment headline retained fraction over the canonical
segments.

ROUNDING.  Q4 says, verbatim: "Do not round a failing 98.997% to 99%."
Condition 3 compares the exact float against 0.99.  Every percentage printed
anywhere in this program is TRUNCATED toward zero at three decimals, and a
retention that fails condition 3 is additionally floored below 99.000 before
it is formatted, so no failing number can ever be displayed as 99%.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from decimal import Decimal, ROUND_DOWN
from math import fsum
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))
from windcheck.manifest import (BaseKind,                         # noqa: E402
                                base_kind_from_manifests,
                                normalise_base_kind)


def _resolved_base_kind(entry: dict) -> str | None:
    """The hash-verified base kind of one base-manifest entry, falling back
    to the declared spelling when content cannot decide. One shared rule,
    shared with the pipeline and the release index."""
    kind = base_kind_from_manifests(entry.get("base_hashes"),
                                    entry.get("original_hashes"))
    if kind is None:
        kind = normalise_base_kind(entry.get("base_kind"))
    return kind.value if kind is not None else entry.get("base_kind")


SCHEMA = "headline_decision/v1"
DEFAULT_BASE_MANIFEST = REPO_ROOT / "out" / "corpus_bases.json"
DEFAULT_CERT_DIR = REPO_ROOT / "out" / "excised" / "corpus"
DEFAULT_SUMMARY = DEFAULT_CERT_DIR / "corpus_summary.jsonl"
CERT_SUFFIX = "_excision_certificate.json"

# Original-area cache.  Re-deriving an unmodified segment's original
# canonical area means re-reading its published mesh; the answer is a pure
# function of that mesh, so it is memoised next to the corpus and re-used on
# every later run.  A cached row is trusted only while the segment still
# names the same original mesh with the same hashes.
DEFAULT_AREA_CACHE = REPO_ROOT / "out" / "headline_original_areas.json"
AREA_CACHE_SCHEMA = "headline_original_areas/v1"

# bench/excise_segment.py CENSUS["maxedge"]: the retained-quad set is the one
# the census used.  Only a fallback -- the value is read off the certificate's
# own census_params whenever it is recorded there.
CENSUS_MAXEDGE = 60.0

EXPECTED_TOTAL = 185
RETENTION_GATE = 0.99          # condition 3: >= 99.000%
PER_SEGMENT_GATE = 0.95        # condition 4: >= 95.000% of ORIGINAL area
WALL_LIMIT_S = 600.0           # condition 7: ten minutes
MAX_LISTED = 20                # offender lists are capped, with a count

# Q4 condition 2: the four dispositions that COUNT as terminal.  "error",
# "residual_transverse" and "not_censusable" are records, but they are not on
# that list, so they fail condition 2 by name.
TERMINAL_DISPOSITIONS = ("transformed", "already_clean", "duplicate_alias",
                         "triangle_empty_invalid")
CLEAN_DISPOSITIONS = ("transformed", "already_clean")

STRONG_SENTENCE = (
    "windcheck produced a reload-verified, transverse-self-intersection-free "
    "tifxyz version of all {total} pinned trace artifacts from five scrolls, "
    "retaining {X}% of original represented surface area overall and at "
    "least {Y}% per trace.")

# The fallback sentence states the CENSUS BOOKKEEPING and nothing else.  The
# artifacts that are not censusable are triangle-empty or invalid INPUTS --
# they were never audited and carry no cleanliness claim either way -- so the
# sentence names them for what they are.  Area retention and fragmentation
# are NOT part of this sentence: they are printed as separate qualifying
# statements underneath it (see `qualifications`), never folded in as a
# clause, because they qualify a different thing from the census.
FALLBACK_SENTENCE = (
    "Of {total} pinned trace artifacts, {censusable} were censusable. All "
    "{censusable} now have reload-verified tifxyz outputs with zero "
    "non-adjacent transverse contacts under both canonical triangulations: "
    "{transformed} were transformed and {unchanged} required no change. "
    "{word} triangle-empty or invalid {input_s} {have} {article}explicit "
    "terminal {record_s}.")

NUMBER_WORDS = ("Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven",
                "Eight", "Nine", "Ten", "Eleven", "Twelve", "Thirteen",
                "Fourteen", "Fifteen", "Sixteen", "Seventeen", "Eighteen",
                "Nineteen", "Twenty")


def count_word(n: int) -> str:
    """A small count spelled out, for the head of a sentence."""
    return NUMBER_WORDS[n] if 0 <= n < len(NUMBER_WORDS) else str(n)


# --------------------------------------------------------------- formatting

def truncate_pct(fraction: float | None, decimals: int = 3) -> str:
    """A fraction as a percentage string, TRUNCATED toward zero.

    Truncation, never rounding: 0.98997 -> "98.997" and 0.98999 -> "98.999",
    never "99.000".  The float is taken through its shortest round-tripping
    decimal representation so that the ordinary binary noise of 1 - a/b does
    not silently shave a digit off a number the corpus really achieved.
    """
    if fraction is None:
        return "n/a"
    q = Decimal(1).scaleb(-decimals)
    pct = Decimal(repr(float(fraction))) * 100
    return str(pct.quantize(q, rounding=ROUND_DOWN))


def fmt_retention(fraction: float | None, passed: bool) -> str:
    """Format a retention figure, with the anti-round-up floor of Q4.

    If the figure FAILED the 99.000% gate it may not be printed as 99% -- so
    a failing value whose truncation still lands on or above 99.000 (possible
    only for a float within an ulp of the gate) is floored to 98.999.
    """
    s = truncate_pct(fraction)
    if not passed and s != "n/a" and Decimal(s) >= Decimal("99.000"):
        return "98.999"
    return s


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _get(d, *path, default=None):
    """Nested lookup that tolerates missing/oddly-typed intermediate nodes."""
    cur = d
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return default if cur is None else cur


def _num(value):
    return value if isinstance(value, (int, float)) and not isinstance(
        value, bool) else None


# ------------------------------------------------------- original-area pricing

def original_canonical_area(mesh: Path, maxedge: float) -> dict:
    """Canonical area of an ORIGINAL published mesh's retained quads.

    This is the denominator `bench/excise_shadow.headline_denominator` builds
    for a segment that WAS cut, computed here by the identical route so the
    two are commensurable: the retained-quad set is
    `windcheck.intrinsic.retained_quads(P, V, maxedge)` -- four valid corners
    and all six edges within the census `maxedge` -- and the per-quad area is
    the canonical one from `bench/excise_segment.quad_area_grids`, the mean of
    the d0 and d1 tessellations.  `math.fsum` makes the total independent of
    iteration order.

    Nothing about the excision enters: with no cut there is no removed-quad
    set to price, so the numerator is exactly zero and only this denominator
    has to be recovered.
    """
    import numpy as np                                     # heavy; on demand
    for extra in (REPO_ROOT / "src", REPO_ROOT / "bench"):
        if str(extra) not in sys.path:
            sys.path.insert(0, str(extra))
    from windcheck import tifxyz
    from windcheck.intrinsic import retained_quads
    from excise_segment import quad_area_grids

    surf = tifxyz.read(Path(mesh))
    P0 = np.asarray(surf.points, np.float64)
    V0 = np.asarray(surf.valid, bool)
    Q0 = retained_quads(P0, V0, maxedge)
    A0 = quad_area_grids(P0)["canonical"]
    return {
        "A_original_canonical": fsum(A0[Q0].tolist()),
        "original_mesh": str(mesh),
        "grid_shape": [int(s) for s in surf.shape],
        "n_original_retained_quads": int(Q0.sum()),
        "maxedge": float(maxedge),
    }


def load_area_cache(path: Path | None) -> dict:
    """The memoised original areas, keyed by segment. Missing file -> {}."""
    if path is None:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    try:
        doc = json.loads(p.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    areas = doc.get("areas") if isinstance(doc, dict) else None
    return areas if isinstance(areas, dict) else {}


def write_area_cache(path: Path, areas: dict) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps({
        "schema": AREA_CACHE_SCHEMA,
        "generated_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "note": ("canonical area of each ORIGINAL published mesh's retained "
                 "quads, for pinned segments whose certificate records no "
                 "cut and therefore no headline_area block; a row is used "
                 "again only while the segment still names the same original "
                 "mesh with the same hashes"),
        "areas": dict(sorted(areas.items())),
    }, indent=1, default=str))


def price_unmodified_segments(index: dict, area_cache: dict,
                              recompute: bool = True) -> dict:
    """Give every UNMODIFIED segment its missing area bookkeeping.

    An `already_clean` certificate is the record of a segment that was left
    alone: no cut, no output mesh, and so no `headline_area` block.  Its two
    headline area numbers are not unknown, they are determined:

        A_removed_priced_on_original = 0.0        (nothing was removed)
        A_original_canonical         = the full original canonical area

    The denominator is recovered exactly as the cut segments recover theirs.
    When the certificate's input mesh IS the original published mesh, the
    certificate already carries that number as `input_area_canonical` -- the
    same shortcut `excise_shadow` takes when it sets `input_is_original` --
    and no mesh has to be re-read.  Otherwise (a displacement-repaired base)
    the original mesh is re-read once and memoised.

    Returns a report; mutates the records in place.
    """
    segs = index["segments"]
    derived: list[dict] = []
    fresh: dict[str, dict] = {}
    failures: list[str] = []

    for seg in index["order"]:
        rec = segs[seg]
        if rec["disposition"] != "already_clean":
            continue
        if rec["A_original_canonical"] is not None:
            continue                       # already priced by the certificate
        cert = rec["certificate"] or {}
        maxedge = _num(_get(cert, "census_params", "maxedge"))
        if maxedge is None:
            maxedge = CENSUS_MAXEDGE
        original = cert.get("original_mesh") or rec.get("original_mesh")
        input_mesh = cert.get("input_mesh") or rec.get("base_mesh")
        hashes = rec.get("original_hashes")

        row = None
        # "The input IS the original" is a statement about CONTENT, decided
        # by the one shared verifier against the recorded manifests. It was
        # once a path comparison, which is wrong: paths differ between a
        # downloaded archive and a fresh workdir while the bytes do not, and
        # two different meshes can sit at the same relative path.
        # `rec["base_kind"]` was resolved by `_resolved_base_kind`: the
        # manifests decide, and the declared spelling is used only where no
        # hashes were recorded at all. Either way it is the same answer the
        # release index and the pipeline give.
        kind = normalise_base_kind(rec.get("base_kind"))
        if kind is None:
            kind = normalise_base_kind(cert.get("base_kind"))
        input_is_original = kind is BaseKind.ORIGINAL
        if original and input_mesh and input_is_original:
            a_in = _num(cert.get("input_area_canonical"))
            if a_in is not None:
                row = {"A_original_canonical": a_in,
                       "original_mesh": str(original),
                       "grid_shape": cert.get("grid_shape"),
                       "n_original_retained_quads": _num(
                           cert.get("n_retained_quads")),
                       "maxedge": float(maxedge),
                       "source": "certificate input_area_canonical "
                                 "(the input IS the original published mesh)"}
        if row is None and original:
            cached = area_cache.get(seg)
            if (isinstance(cached, dict)
                    and _num(cached.get("A_original_canonical")) is not None
                    and cached.get("original_mesh") == str(original)
                    and cached.get("original_hashes") == hashes
                    and _num(cached.get("maxedge")) == float(maxedge)):
                row = dict(cached)
                row["source"] = "cached recomputation from the original mesh"
            elif recompute:
                try:
                    row = original_canonical_area(Path(original), maxedge)
                except Exception as exc:                   # unreadable mesh
                    failures.append(f"{seg} ({exc!r})")
                    row = None
                if row is not None:
                    row["original_hashes"] = hashes
                    fresh[seg] = dict(row)
                    row["source"] = "recomputed from the original mesh"

        if row is None:
            if not failures or not failures[-1].startswith(seg):
                failures.append(f"{seg} (no original mesh to price it on)")
            continue

        rec["A_original_canonical"] = row["A_original_canonical"]
        rec["A_removed_priced_on_original"] = 0.0
        rec["area_source"] = row["source"]
        # Unmodified means unmodified: the retained fraction is exactly 1.0.
        if rec["headline_retained_fraction"] is None:
            rec["headline_retained_fraction"] = 1.0
        if rec["operational_retained_fraction"] is None:
            rec["operational_retained_fraction"] = 1.0
        derived.append({"segment": seg, **row})

    return {"derived": derived, "fresh": fresh, "failures": failures}


# ------------------------------------------------------------------ loading

def load_base_manifest(path: Path) -> list[dict]:
    """Per-segment base records; accepts the list or dict entry shapes."""
    doc = json.loads(Path(path).read_text())
    entries = doc.get("entries", doc) if isinstance(doc, dict) else doc
    if isinstance(entries, dict):
        rows = []
        for seg, rec in entries.items():
            rec = dict(rec or {})
            rec.setdefault("segment", seg)
            rows.append(rec)
        entries = rows
    if not isinstance(entries, list):
        raise ValueError(f"{path}: unrecognised base-manifest shape")
    out = []
    for rec in entries:
        if not isinstance(rec, dict) or not rec.get("segment"):
            raise ValueError(f"{path}: entry without a segment name")
        out.append(rec)
    return out


def load_summary(path: Path) -> dict[str, dict]:
    """The corpus summary index, keyed by segment. Missing file -> {}."""
    p = Path(path)
    if not p.exists():
        return {}
    rows: dict[str, dict] = {}
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        seg = str(rec.get("segment") or rec.get("seg") or "")
        if seg:
            rows[seg] = rec
    return rows


def certificate_path(segment: str, row: dict, cert_dir: Path) -> Path:
    """Where this segment's certificate is claimed to live."""
    for key in ("certificate", "certificate_path", "cert", "cert_path",
                "record", "record_path"):
        val = row.get(key)
        if isinstance(val, str) and val:
            p = Path(val)
            return p if p.is_absolute() else (REPO_ROOT / p)
    return Path(cert_dir) / f"{segment}{CERT_SUFFIX}"


def summary_hash(row: dict) -> str | None:
    for key in ("certificate_sha256", "cert_sha256", "sha256",
                "certificate_hash"):
        val = row.get(key)
        if isinstance(val, str) and val:
            return val.lower()
    return None


# ------------------------------------------------------------- per-segment

def collect(base_manifest: Path, summary_path: Path, cert_dir: Path,
            strict_artifacts: bool, area_cache: Path | None = None,
            recompute_areas: bool = True) -> dict:
    """One normalised record per pinned segment, plus artifact problems."""
    bases = load_base_manifest(base_manifest)
    summary = load_summary(summary_path)
    cert_dir = Path(cert_dir)

    segs: dict[str, dict] = {}
    order: list[str] = []
    for base in bases:
        seg = str(base["segment"])
        row = summary.get(seg, {})
        cert_p = certificate_path(seg, row, cert_dir)
        rec: dict = {
            "segment": seg,
            "is_canonical": bool(base.get("is_canonical", True)),
            "duplicate_of": base.get("duplicate_of"),
            "geometry_key": base.get("geometry_key"),
            # base_kind is resolved by content, through the one shared
            # verifier, not by comparing paths.
            "base_kind": _resolved_base_kind(base),
            "base_kind_declared": base.get("base_kind"),
            "original_mesh": base.get("original_mesh"),
            "original_hashes": base.get("original_hashes"),
            "base_mesh": base.get("base_mesh"),
            "input_hashes": base.get("base_hashes"),
            "in_summary": seg in summary,
            "certificate_path": str(cert_p),
            "certificate_present": cert_p.exists(),
            "artifact_problems": [],
            "certificate": None,
        }
        if not rec["certificate_present"]:
            rec["artifact_problems"].append("missing certificate")
        else:
            try:
                rec["certificate"] = json.loads(cert_p.read_text())
            except Exception as exc:                       # unreadable == missing
                rec["certificate_present"] = False
                rec["artifact_problems"].append(f"unreadable certificate: {exc!r}")
            want = summary_hash(row)
            if strict_artifacts and rec["certificate_present"]:
                got = sha256_file(cert_p)
                rec["certificate_sha256"] = got
                if want and got != want:
                    rec["artifact_problems"].append(
                        f"certificate hash mismatch (summary {want[:12]}..., "
                        f"file {got[:12]}...)")
        segs[seg] = rec
        order.append(seg)

    for seg in order:
        normalise(segs[seg])

    # A duplicate alias inherits the geometry -- and therefore the
    # cleanliness and the areas -- of the canonical it points at.
    for seg in order:
        rec = segs[seg]
        canon = rec.get("duplicate_of")
        rec["canonical_segment"] = seg if rec["is_canonical"] else canon
        rec["canonical_record"] = segs.get(canon) if canon in segs else None
        if canon and canon not in segs:
            rec["artifact_problems"].append(
                f"duplicate_of names an unknown segment {canon!r}")

    index = {"segments": segs, "order": order, "summary_extras": extras_of(
        summary, order), "n_summary_rows": len(summary)}

    # An unmodified segment's certificate records no cut and therefore no
    # headline_area.  Price it now, so that condition 3 and condition 4 see
    # EVERY pinned entry rather than only the ones that were cut.
    cache = load_area_cache(area_cache)
    report = price_unmodified_segments(index, cache, recompute=recompute_areas)
    if report["fresh"] and area_cache is not None:
        cache.update(report["fresh"])
        write_area_cache(area_cache, cache)
    index["area_pricing"] = {
        "derived": report["derived"],
        "n_derived": len(report["derived"]),
        "failures": report["failures"],
        "cache": str(area_cache) if area_cache is not None else None,
        "n_recomputed_from_mesh": len(report["fresh"]),
    }
    return index


def extras_of(summary: dict, order: list[str]) -> list[str]:
    """Segments present in the summary but absent from the pinned roster."""
    return sorted(set(summary) - set(order))


def normalise(rec: dict) -> None:
    """Pull the fields the eight conditions need out of one certificate."""
    cert = rec["certificate"] or {}
    rec["disposition"] = cert.get("terminal_disposition")
    rec["status"] = cert.get("status")
    rec["record_kind"] = cert.get("record_kind")
    rec["policy_hash"] = cert.get("policy_hash")
    rec["wall_seconds"] = _num(cert.get("wall_seconds"))

    # -- census (condition 1 / 8)
    d0 = _num(_get(cert, "census_after", "d0", "transverse"))
    d1 = _num(_get(cert, "census_after", "d1", "transverse"))
    rec["census_after"] = {"d0": d0, "d1": d1}
    rec["census_after_present"] = d0 is not None and d1 is not None
    total = _num(cert.get("output_transverse_total"))
    if total is None and rec["census_after_present"]:
        total = d0 + d1
    rec["output_transverse_total"] = total

    # -- reload checks (condition 1 / 8): EVERY boolean must be true
    checks = cert.get("reload_checks")
    rec["reload_checks_present"] = isinstance(checks, dict) and any(
        isinstance(v, bool) for v in checks.values())
    rec["reload_failures"] = sorted(
        k for k, v in (checks or {}).items() if isinstance(v, bool) and not v
    ) if isinstance(checks, dict) else []

    # -- areas (conditions 3 / 4)
    a_orig = _num(_get(cert, "headline_area", "A_original_canonical"))
    a_rem = _num(_get(cert, "headline_area", "A_removed_priced_on_original"))
    if a_orig is None:
        a_orig = _num(_get(cert, "area", "canonical", "A_input"))
        if a_orig is not None and a_rem is None:
            a_rem = _num(_get(cert, "area", "canonical", "A_excised"))
    if a_rem is None and a_orig is not None and \
            rec["disposition"] == "already_clean":
        a_rem = 0.0                     # nothing was cut, by definition
    rec["A_original_canonical"] = a_orig
    rec["A_removed_priced_on_original"] = a_rem

    hrf = _num(cert.get("headline_retained_fraction"))
    if hrf is None and a_orig:
        hrf = 1.0 - (a_rem / a_orig) if a_rem is not None else None
    if hrf is None and rec["disposition"] == "already_clean":
        hrf = _num(_get(cert, "area", "canonical", "retained_fraction"))
    rec["headline_retained_fraction"] = hrf
    rec["operational_retained_fraction"] = _num(
        cert.get("operational_retained_fraction"))

    # -- core gate (condition 5)
    gate = cert.get("core_gate_pass")
    if not isinstance(gate, bool):
        gate = _get(cert, "component_recovery", "core_gate", "core_gate_pass")
    rec["core_gate_pass"] = gate if isinstance(gate, bool) else None
    rec["min_R_main_core"] = _num(
        _get(cert, "component_recovery", "core_gate", "min_R_main_core"))

    # -- emptiness guard (condition 6)
    empty = _get(cert, "emptiness_guard", "clean_by_emptiness")
    rec["clean_by_emptiness"] = empty if isinstance(empty, bool) else None

    # A transformed segment is a production artifact: the evidence for every
    # gate must be ON the certificate.  Silence is not a pass.
    if rec["disposition"] == "transformed":
        for field, label in (("census_after_present", "census_after"),
                             ("reload_checks_present", "reload_checks")):
            if not rec[field]:
                rec["artifact_problems"].append(f"{label} missing")
        for key, label in (("A_original_canonical", "headline_area"),
                           ("core_gate_pass", "core_gate_pass"),
                           ("clean_by_emptiness", "emptiness_guard"),
                           ("wall_seconds", "wall_seconds")):
            if rec[key] is None:
                rec["artifact_problems"].append(f"{label} missing")


# ------------------------------------------------------------- the numbers

def retention_figures(index: dict) -> dict:
    """Unique-geometry-weighted and artifact-count-weighted retention."""
    segs, order = index["segments"], index["order"]
    uni_num = uni_den = 0.0
    art_num = art_den = 0.0
    gaps_unique: list[str] = []
    gaps_artifact: list[str] = []
    n_unique = n_artifact = 0

    def priced(rec):
        """(A_original, A_removed) for this artifact, or None if unknown."""
        a_o, a_r = rec["A_original_canonical"], rec["A_removed_priced_on_original"]
        if a_o is None and rec["disposition"] == "duplicate_alias":
            canon = rec.get("canonical_record")
            if canon is not None:
                a_o = canon["A_original_canonical"]
                a_r = canon["A_removed_priced_on_original"]
        if a_o is None:
            return None
        return a_o, (a_r or 0.0)

    for seg in order:
        rec = segs[seg]
        pair = priced(rec)
        # An explicitly triangle-empty/invalid record carries no surface: it
        # is neither numerator nor denominator, and that is not a data gap.
        skip_ok = rec["disposition"] in ("triangle_empty_invalid",)
        if rec["is_canonical"]:
            n_unique += 1
            if pair is None:
                if not skip_ok:
                    gaps_unique.append(seg)
            else:
                uni_den += pair[0]
                uni_num += pair[1]
        n_artifact += 1
        if pair is None:
            if not skip_ok:
                gaps_artifact.append(seg)
        else:
            art_den += pair[0]
            art_num += pair[1]

    uni = (1.0 - uni_num / uni_den) if uni_den > 0 else None
    art = (1.0 - art_num / art_den) if art_den > 0 else None

    per_seg = {seg: segs[seg]["headline_retained_fraction"] for seg in order
               if segs[seg]["is_canonical"]
               and segs[seg]["headline_retained_fraction"] is not None}
    min_seg = min(per_seg, key=per_seg.get) if per_seg else None

    return {
        "unique_geometry_weighted": uni,
        "unique_A_original_total": uni_den,
        "unique_A_removed_total": uni_num,
        "unique_segments_counted": n_unique - len(gaps_unique),
        "unique_segments": n_unique,
        "unique_missing_area": gaps_unique,
        "artifact_count_weighted": art,
        "artifact_A_original_total": art_den,
        "artifact_A_removed_total": art_num,
        "artifact_segments": n_artifact,
        "artifact_missing_area": gaps_artifact,
        "min_headline_retained_fraction": per_seg.get(min_seg),
        "min_headline_segment": min_seg,
    }


def census_counts(index: dict) -> dict:
    """The census bookkeeping the fallback sentence states, from the data.

    An alias is counted under the disposition of the canonical geometry it
    aliases, because that is the artifact it is a copy of.  A triangle-empty
    or invalid input is NOT censusable: it carries no triangles, so it was
    never audited and no cleanliness claim is defined on it either way.
    """
    segs = index["segments"]
    n = {"transformed": 0, "already_clean": 0, "not_censusable": 0,
         "other": 0}
    for seg in index["order"]:
        rec = segs[seg]
        disp = rec["disposition"]
        if disp == "duplicate_alias":
            canon = rec.get("canonical_record")
            disp = canon["disposition"] if canon is not None else None
        if disp == "transformed":
            n["transformed"] += 1
        elif disp == "already_clean":
            n["already_clean"] += 1
        elif disp == "triangle_empty_invalid":
            n["not_censusable"] += 1
        else:
            n["other"] += 1
    n["censusable"] = n["transformed"] + n["already_clean"]
    return n


def clean_artifact_count(index: dict) -> int:
    """N: pinned artifacts delivered transverse-clean (aliases included)."""
    segs, n = index["segments"], 0
    for seg in index["order"]:
        rec = segs[seg]
        if rec["disposition"] in CLEAN_DISPOSITIONS:
            n += 1
        elif rec["disposition"] == "duplicate_alias":
            canon = rec.get("canonical_record")
            if canon is not None and canon["disposition"] in CLEAN_DISPOSITIONS:
                n += 1
    return n


# ------------------------------------------------------- the eight conditions

def condition_text(expected_total: int) -> dict[int, str]:
    """The frozen Q4 text, transcribed verbatim from DECISIONS.md."""
    return {
        1: ("every unique censusable mesh emits a reload-verified aggregate "
            "with transverse 0/0 both diagonals"),
        2: (f"all {expected_total} entries have a terminal disposition "
            "(transformed | already-clean | duplicate alias | explicitly "
            "triangle-empty/invalid)"),
        3: "unique-geometry-weighted corpus retention >= 99.000%",
        4: ("every claimed-clean segment >= 95.000% of ORIGINAL canonical "
            "area"),
        5: ("every component in every segment's input 99.9%-area core has "
            "R_main >= 0.90"),
        6: "no segment clean by empty/near-empty output",
        7: ("every production artifact completes within ten minutes under "
            "the frozen policy"),
        8: ("no missing certificate, failed reload, unresolved transverse "
            "contact, timeout or unaccounted hash mismatch"),
    }


def evaluate(index: dict, retention: dict, expected_total: int) -> list[dict]:
    segs, order = index["segments"], index["order"]
    text = condition_text(expected_total)
    conds: list[dict] = []

    def add(n, passed, evidence, offenders=()):
        conds.append({"n": n, "text": text[n], "pass": bool(passed),
                      "evidence": evidence,
                      "offenders": [str(o) for o in offenders],
                      "n_offenders": len(list(offenders))})

    # -- 1 -----------------------------------------------------------------
    off1, checked1 = [], 0
    for seg in order:
        rec = segs[seg]
        if not (rec["is_canonical"] and rec["disposition"] == "transformed"):
            continue
        checked1 += 1
        why = []
        if not rec["census_after_present"]:
            why.append("no census_after")
        else:
            if rec["census_after"]["d0"] != 0:
                why.append(f"d0 transverse {rec['census_after']['d0']}")
            if rec["census_after"]["d1"] != 0:
                why.append(f"d1 transverse {rec['census_after']['d1']}")
        if not rec["reload_checks_present"]:
            why.append("no reload_checks")
        if rec["reload_failures"]:
            why.append("reload " + ",".join(rec["reload_failures"]))
        if why:
            off1.append(f"{seg} ({'; '.join(why)})")
    add(1, not off1,
        f"{checked1 - len(off1)}/{checked1} transformed unique meshes are "
        f"reload-verified with transverse 0/0 on both diagonals", off1)

    # -- 2 -----------------------------------------------------------------
    off2 = []
    for seg in order:
        d = segs[seg]["disposition"]
        if d not in TERMINAL_DISPOSITIONS:
            off2.append(f"{seg} ({d if d else 'no certificate/disposition'})")
    roster_ok = len(order) == expected_total
    if not roster_ok:
        off2.append(f"roster has {len(order)} entries, expected "
                    f"{expected_total}")
    if index["summary_extras"]:
        off2.extend(f"{s} (in summary but not on the pinned roster)"
                    for s in index["summary_extras"])
    counts = {}
    for seg in order:
        counts[segs[seg]["disposition"] or "none"] = counts.get(
            segs[seg]["disposition"] or "none", 0) + 1
    add(2, not off2,
        f"{len(order)} entries; dispositions " + ", ".join(
            f"{k}={v}" for k, v in sorted(counts.items())), off2)

    # -- 3 -----------------------------------------------------------------
    uni = retention["unique_geometry_weighted"]
    gaps = retention["unique_missing_area"]
    passed3 = uni is not None and not gaps and uni >= RETENTION_GATE
    ev3 = (f"unique-geometry-weighted retention = "
           f"{fmt_retention(uni, passed3)}% "
           f"(1 - {retention['unique_A_removed_total']!r} / "
           f"{retention['unique_A_original_total']!r} over "
           f"{retention['unique_segments_counted']} canonical segments); "
           f"gate 99.000%")
    if gaps:
        ev3 += f"; {len(gaps)} canonical segments have no priced area"
    add(3, passed3, ev3, gaps)

    # -- 4 -----------------------------------------------------------------
    off4, checked4 = [], 0
    for seg in order:
        rec = segs[seg]
        frac = rec["headline_retained_fraction"]
        if rec["disposition"] not in CLEAN_DISPOSITIONS:
            if rec["disposition"] == "duplicate_alias":
                canon = rec.get("canonical_record")
                if canon is None or canon["disposition"] not in CLEAN_DISPOSITIONS:
                    continue
                frac = canon["headline_retained_fraction"]
            else:
                continue
        checked4 += 1
        if frac is None:
            off4.append(f"{seg} (no headline retained fraction)")
        elif frac < PER_SEGMENT_GATE:
            off4.append(f"{seg} ({truncate_pct(frac)}%)")
    worst = retention["min_headline_retained_fraction"]
    add(4, not off4,
        f"{checked4 - len(off4)}/{checked4} claimed-clean artifacts keep "
        f">= 95.000% of ORIGINAL canonical area; minimum over canonical "
        f"segments {truncate_pct(worst)}% "
        f"({retention['min_headline_segment'] or 'n/a'})", off4)

    # -- 5 -----------------------------------------------------------------
    off5, checked5 = [], 0
    for seg in order:
        rec = segs[seg]
        gate = rec["core_gate_pass"]
        if gate is None:
            if rec["disposition"] == "transformed":
                checked5 += 1
                off5.append(f"{seg} (core gate not reported)")
            continue
        checked5 += 1
        if not gate:
            off5.append(f"{seg} (min R_main core "
                        f"{rec['min_R_main_core']})")
    add(5, not off5,
        f"{checked5 - len(off5)}/{checked5} segments where the 99.9%-area "
        f"core gate applies pass it", off5)

    # -- 6 -----------------------------------------------------------------
    off6, checked6 = [], 0
    for seg in order:
        rec = segs[seg]
        flag = rec["clean_by_emptiness"]
        if flag is None:
            if rec["disposition"] == "transformed":
                checked6 += 1
                off6.append(f"{seg} (emptiness guard not reported)")
            continue
        checked6 += 1
        if flag:
            off6.append(f"{seg} (clean_by_emptiness, operational retained "
                        f"{truncate_pct(rec['operational_retained_fraction'])}%)")
    add(6, not off6,
        f"{checked6 - len(off6)}/{checked6} emitted segments are clean on "
        f"their own geometry, none by empty/near-empty output", off6)

    # -- 7 -----------------------------------------------------------------
    off7, checked7, slowest, slowseg = [], 0, None, None
    for seg in order:
        rec = segs[seg]
        wall = rec["wall_seconds"]
        if wall is None:
            if rec["disposition"] == "transformed":
                checked7 += 1
                off7.append(f"{seg} (no wall_seconds)")
            continue
        checked7 += 1
        if slowest is None or wall > slowest:
            slowest, slowseg = wall, seg
        if wall > WALL_LIMIT_S:
            off7.append(f"{seg} ({wall}s)")
    add(7, not off7,
        f"{checked7 - len(off7)}/{checked7} production artifacts finish "
        f"within {WALL_LIMIT_S:.0f}s; slowest "
        f"{slowest if slowest is not None else 'n/a'}s ({slowseg or 'n/a'})",
        off7)

    # -- 8 -----------------------------------------------------------------
    off8 = []
    for seg in order:
        rec = segs[seg]
        why = list(rec["artifact_problems"])
        if rec["reload_failures"]:
            why.append("failed reload: " + ",".join(rec["reload_failures"]))
        if rec["output_transverse_total"]:
            why.append(f"unresolved transverse contact "
                       f"({rec['output_transverse_total']})")
        if rec["disposition"] == "residual_transverse":
            why.append("residual transverse disposition")
        if rec["wall_seconds"] is not None and rec["wall_seconds"] > WALL_LIMIT_S:
            why.append(f"timeout ({rec['wall_seconds']}s)")
        if rec["disposition"] == "error":
            why.append("error disposition")
        if why:
            off8.append(f"{seg} ({'; '.join(why)})")
    if index["summary_extras"]:
        off8.extend(f"{s} (summary row with no pinned base entry)"
                    for s in index["summary_extras"])
    n_present = sum(1 for s in order if segs[s]["certificate_present"])
    add(8, not off8,
        f"{n_present}/{len(order)} certificates present and hash-checked; "
        f"{len(off8)} artifact problems", off8)

    return conds


# ---------------------------------------------------------------- reporting

def build_qualifications(index: dict, retention: dict, conds: list[dict],
                         X: str, Y: str) -> list[dict]:
    """The qualifying statements that sit BELOW the headline sentence.

    Area retention and fragmentation qualify something different from the
    census the sentence states, so they are separate statements and are never
    folded into it -- not as a clause, not as a trailing "while ...".  A
    reader who quotes the sentence quotes the census; a reader who quotes a
    qualification quotes the qualification.
    """
    segs = index["segments"]
    cond3 = next(c for c in conds if c["n"] == 3)
    cond5 = next(c for c in conds if c["n"] == 5)

    counted = retention["unique_segments_counted"]
    area = (f"Area retention is a separate figure: unique-geometry-weighted "
            f"retention of original represented surface area over the "
            f"{counted} priced canonical segments is {X}%")
    if not cond3["pass"]:
        area += f", below the gate {truncate_pct(RETENTION_GATE)}%"
    worst_seg = retention["min_headline_segment"] or "n/a"
    area += (f", and the lowest per-segment retained fraction is {Y}% "
             f"({worst_seg}).")

    if cond5["pass"]:
        frag = ("Fragmentation is a separate figure: every component in "
                "every segment's input 99.9%-area core meets the R_main "
                ">= 0.90 recovery gate.")
    else:
        worst = None
        for entry in cond5["offenders"]:
            seg = entry.split(" (")[0]
            val = _num((segs.get(seg) or {}).get("min_R_main_core"))
            if val is not None and (worst is None or val < worst[1]):
                worst = (seg, val)
        n = cond5["n_offenders"]
        frag = (f"Fragmentation is a separate figure: {n} "
                f"{'segment does' if n == 1 else 'segments do'} not meet the "
                f"R_main >= 0.90 recovery gate on every component of the "
                f"input 99.9%-area core")
        frag += (f"; the lowest reported core R_main is {worst[1]!r} "
                 f"({worst[0]})." if worst else ".")

    return [{"kind": "area_retention", "text": area},
            {"kind": "fragmentation", "text": frag}]


def decide(index: dict, expected_total: int) -> dict:
    retention = retention_figures(index)
    conds = evaluate(index, retention, expected_total)
    strong = all(c["pass"] for c in conds)
    cond3 = next(c for c in conds if c["n"] == 3)

    X = fmt_retention(retention["unique_geometry_weighted"], cond3["pass"])
    Y = truncate_pct(retention["min_headline_retained_fraction"])
    N = clean_artifact_count(index)
    M = expected_total - N
    census = census_counts(index)

    if strong:
        sentence = STRONG_SENTENCE.format(total=expected_total, X=X, Y=Y)
    else:
        nc = census["not_censusable"]
        sentence = FALLBACK_SENTENCE.format(
            total=expected_total,
            censusable=census["censusable"],
            transformed=census["transformed"],
            unchanged=census["already_clean"],
            word=count_word(nc),
            input_s="input" if nc == 1 else "inputs",
            have="has" if nc == 1 else "have",
            article="an " if nc == 1 else "",
            record_s="record" if nc == 1 else "records")
        if census["other"]:
            sentence += (f" {count_word(census['other'])} further "
                         f"{'entry has' if census['other'] == 1 else 'entries have'} "
                         f"no terminal disposition.")
    return {
        "schema": SCHEMA,
        "generated_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "decision_rule": ("DECISIONS.md 2026-07-31 Round 28 Q4, frozen "
                          "before the corpus pass"),
        "verdict": "STRONG" if strong else "FALLBACK",
        "sentence": sentence,
        "sentence_states": ("the census bookkeeping only; area retention and "
                            "fragmentation are separate qualifying "
                            "statements and are never folded into it"),
        "qualifications": build_qualifications(index, retention, conds, X, Y),
        "conditions": conds,
        "n_conditions_passed": sum(1 for c in conds if c["pass"]),
        "retention": {
            **retention,
            "X_percent": X,
            "Y_percent": Y,
            "artifact_count_weighted_percent": truncate_pct(
                retention["artifact_count_weighted"]),
            "rounding": ("percentages are TRUNCATED toward zero at three "
                         "decimals; a retention that fails condition 3 is "
                         "never displayed as 99%"),
        },
        "counts": {
            "expected_total": expected_total,
            "roster_entries": len(index["order"]),
            "summary_rows": index["n_summary_rows"],
            "canonical": sum(1 for s in index["order"]
                             if index["segments"][s]["is_canonical"]),
            "duplicate_aliases": sum(1 for s in index["order"]
                                     if not index["segments"][s]["is_canonical"]),
            "transverse_clean_artifacts_N": N,
            "other_dispositions_M": M,
            "censusable": census["censusable"],
            "transformed": census["transformed"],
            "already_clean_unchanged": census["already_clean"],
            "not_censusable_triangle_empty_invalid": census["not_censusable"],
            "no_terminal_disposition": census["other"],
        },
        "area_pricing": index.get("area_pricing", {}),
    }


def render(result: dict, out=None) -> None:
    # Resolved at call time, not import time, so a redirected stdout is used.
    stream = sys.stdout if out is None else out
    p = lambda s="": print(s, file=stream)                     # noqa: E731
    p("=" * 78)
    p("ROUND-28 Q4 HEADLINE DECISION RULE (frozen 2026-07-31)")
    p("=" * 78)
    for c in result["conditions"]:
        p(f"COND {c['n']} [{'PASS' if c['pass'] else 'FAIL'}] {c['text']}")
        p(f"    evidence: {c['evidence']}")
        if c["offenders"]:
            shown = c["offenders"][:MAX_LISTED]
            p(f"    offending segments ({c['n_offenders']}):")
            for o in shown:
                p(f"      - {o}")
            if c["n_offenders"] > len(shown):
                p(f"      ... and {c['n_offenders'] - len(shown)} more")
        p()
    r, k = result["retention"], result["counts"]
    pct = lambda s: s if s == "n/a" else f"{s}%"                # noqa: E731
    p(f"unique-geometry-weighted retention : {pct(r['X_percent'])}  "
      f"({r['unique_segments_counted']}/{k['canonical']} canonical segments, "
      f"aliases counted once)")
    p(f"artifact-count-weighted retention  : "
      f"{pct(r['artifact_count_weighted_percent'])}  "
      f"({k['roster_entries']} pinned artifacts)")
    p(f"minimum per-segment headline       : {pct(r['Y_percent'])}  "
      f"({r['min_headline_segment'] or 'n/a'})")
    p(f"transverse-clean artifacts         : "
      f"{k['transverse_clean_artifacts_N']} of {k['expected_total']} "
      f"(other dispositions {k['other_dispositions_M']})")
    p(f"conditions passed                  : "
      f"{result['n_conditions_passed']}/8")
    ap = result.get("area_pricing") or {}
    if ap.get("n_derived"):
        p(f"unmodified segments priced at 1.0  : {ap['n_derived']}  "
          f"({ap.get('n_recomputed_from_mesh', 0)} recomputed from the "
          f"original mesh this run)")
    if ap.get("failures"):
        p(f"unmodified segments still unpriced : {len(ap['failures'])}")
    p()
    p(f"VERDICT: {result['verdict']}")
    p("SENTENCE:")
    p(result["sentence"])
    p()
    p("QUALIFICATIONS (separate statements, NOT part of the sentence above):")
    for q in result.get("qualifications", []):
        p(f"  - {q['text']}")
    p("=" * 78)


# --------------------------------------------------------------------- main

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=("Evaluate the frozen round-28 Q4 headline decision rule "
                     "automatically from the corpus index."))
    ap.add_argument("--base-manifest", type=Path, default=DEFAULT_BASE_MANIFEST)
    ap.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    ap.add_argument("--certificates-dir", type=Path, default=DEFAULT_CERT_DIR)
    ap.add_argument("--expected-total", type=int, default=EXPECTED_TOTAL,
                    help="pinned artifact count named in the frozen sentence")
    ap.add_argument("--json", type=Path, default=None,
                    help="write the machine-readable verdict here")
    ap.add_argument("--strict-artifacts", dest="strict_artifacts",
                    action="store_true", default=True,
                    help="(default) every referenced certificate must exist "
                         "and hash-match")
    ap.add_argument("--no-strict-artifacts", dest="strict_artifacts",
                    action="store_false")
    ap.add_argument("--area-cache", type=Path, default=DEFAULT_AREA_CACHE,
                    help="memoised original canonical areas for unmodified "
                         "segments whose certificate records no cut")
    ap.add_argument("--no-area-cache", dest="area_cache",
                    action="store_const", const=None,
                    help="neither read nor write the area cache")
    ap.add_argument("--recompute-areas", dest="recompute_areas",
                    action="store_true", default=True,
                    help="(default) re-read an original mesh when an "
                         "unmodified segment's area is not already known")
    ap.add_argument("--no-area-recompute", dest="recompute_areas",
                    action="store_false")
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        index = collect(args.base_manifest, args.summary,
                        args.certificates_dir, args.strict_artifacts,
                        area_cache=args.area_cache,
                        recompute_areas=args.recompute_areas)
    except FileNotFoundError as exc:
        print(f"[input error] {exc}", file=sys.stderr)
        return 2
    except (ValueError, json.JSONDecodeError) as exc:
        print(f"[input error] {exc}", file=sys.stderr)
        return 2

    result = decide(index, args.expected_total)
    result["inputs"] = {
        "base_manifest": str(args.base_manifest),
        "summary": str(args.summary),
        "certificates_dir": str(args.certificates_dir),
        "strict_artifacts": bool(args.strict_artifacts),
        "area_cache": str(args.area_cache) if args.area_cache else None,
        "recompute_areas": bool(args.recompute_areas),
    }
    render(result)
    if args.json:
        Path(args.json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json).write_text(json.dumps(result, indent=1, default=str))
        print(f"wrote {args.json}")
    return 0 if result["verdict"] == "STRONG" else 1


if __name__ == "__main__":
    sys.exit(main())
