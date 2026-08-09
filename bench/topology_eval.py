"""Score a candidate derivative against a benchmark trace (PREREG B4/B5).

The benchmark ships paired inputs and reference derivatives. This is the
part that lets someone else's method be measured the same way ours was.

Design decisions, all fixed in the prereg before any of this ran:

  clean is primary and binary. A candidate that still self-intersects is
  NOT scored on cost. Reporting "retained 99.9%" for a surface that never
  satisfied the invariant would rank a method that did nothing at the top.

  no single scalar. A leaderboard number invites trading fragmentation
  for retention. The evaluator emits the vector and refuses to collapse
  it; `--json` gives the whole thing.

  the reference is one admissible answer, not the target. A candidate
  that keeps more area while staying clean and unfragmented is BETTER,
  and the comparison says so rather than scoring distance-to-reference.

  same input or nothing. The candidate is checked against the input by
  hash. A method evaluated on different bytes is not evaluated.

  moving geometry is not cheating, but it is not excision either. A
  candidate whose retained cells differ from the input is reported with
  coordinate_fidelity=false and its retention is NOT comparable
  like-for-like; the report says so in words, not just a flag.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import tifffile

CENSUS = {"exclude": 1, "maxedge": 60.0, "cell": 40.0,
          "touch_tolerance": 0.001}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def load_surface(d: Path):
    """Load a tifxyz under CERTIFICATE semantics, and say so.

    Certificate semantics -- what every published retention, excision and
    fragmentation figure in this project was measured under -- are: a cell
    is valid unless it carries the `x = y = z = -1` marker, unless a
    coordinate is non-finite, or unless `mask.tif` says below 255.

    This is NOT identical to the official loader, which additionally
    invalidates `z <= 0`. On the traces that carry retained geometry at
    `z <= 0` our measurement covers a SUPERSET of what the official
    validator examines. That divergence is reported per surface as
    `n_retained_cells_at_z_le_0` rather than quietly reconciled, because
    the two loaders are both right about different things and picking one
    silently would misstate whichever figure the reader assumed.

    A mask whose shape does not match the coordinate planes used to be
    ignored. That is how a malformed artifact gets scored as if it were
    well-formed, so it now raises.
    """
    arrs = []
    for a in "xyz":
        raw = np.asarray(tifffile.imread(d / f"{a}.tif"))
        if raw.dtype != np.float32:
            # Recorded, not coerced: "bit-identical" downstream must mean
            # the decoded values, not values after a silent conversion.
            raise ValueError(f"{d}/{a}.tif has dtype {raw.dtype}, expected "
                             "float32; comparisons in this tool assume the "
                             "native decoded dtype")
        arrs.append(raw)
    x, y, z = arrs
    valid = ~((x == -1) & (y == -1) & (z == -1))
    valid &= np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    mp = d / "mask.tif"
    if mp.is_file():
        m = np.asarray(tifffile.imread(mp))
        if m.ndim == 3:
            m = m[..., 0]
        if m.shape != valid.shape:
            raise ValueError(f"{mp} has shape {m.shape}, coordinate planes "
                             f"are {valid.shape}")
        valid &= m >= 255
    return x, y, z, valid


def quad_areas(x, y, z, valid):
    """Canonical quad area on the (v,u) lattice, both triangulations
    averaged -- the same pricing the corpus used.

    A quad is retained on the SAME criterion the census applies: four
    valid corners AND all six pairwise corner distances within `maxedge`.
    The edge-length rule is not decoration. `vc_tifxyz_selfcross` drops
    long quads before it builds triangles, so a quad over the limit is not
    part of the censused surface at all. Pricing it would credit a
    candidate with area the invariant never examined, and would let a
    dropped quad bridge two components that the census sees as separate.
    Omitting this rule put the evaluator 1e-5 away from the certificates
    on 87 of 274 traces, which is how it was caught.
    """
    P = np.stack([x, y, z], axis=-1).astype(np.float64)
    ok = valid[:-1, :-1] & valid[:-1, 1:] & valid[1:, :-1] & valid[1:, 1:]
    p00, p10 = P[:-1, :-1], P[1:, :-1]
    p01, p11 = P[:-1, 1:], P[1:, 1:]
    maxedge = float(CENSUS["maxedge"])
    if maxedge > 0:
        e = np.zeros(ok.shape)
        for u, v in ((p00, p01), (p01, p11), (p11, p10), (p00, p10),
                     (p00, p11), (p10, p01)):
            e = np.maximum(e, np.linalg.norm(u - v, axis=-1))
        ok = ok & (e <= maxedge)
    a, b = P[:-1, :-1], P[:-1, 1:]
    c, d = P[1:, :-1], P[1:, 1:]
    def tri(p, q, r):
        return 0.5 * np.linalg.norm(np.cross(q - p, r - p), axis=-1)
    d0 = tri(a, b, c) + tri(b, d, c)
    d1 = tri(a, b, d) + tri(a, d, c)
    area = 0.5 * (d0 + d1)
    area[~ok] = 0.0
    return area, ok


def components(ok: np.ndarray):
    from scipy import ndimage
    lab, n = ndimage.label(ok, structure=np.ones((3, 3), bool))
    return lab, n


def _int(v):
    """An integer, and not a bool. `True == 1` in Python, so a report
    carrying `"transverse": false` would otherwise read as zero."""
    return None if isinstance(v, bool) or not isinstance(v, int) else v


def read_census_report(r: dict) -> dict:
    """Parse a census report STRICTLY. Nothing is coerced.

    Every relaxation here is a way for a malformed or hostile report to
    be read as a clean verdict: a single diagonal passing as "both", a
    missing count read as zero, a string coerced to an integer, or -- the
    shape the mask defect took -- a clean verdict over zero triangles,
    which is true of any empty surface.
    """
    problems: list[str] = []
    entries = r.get("census")
    if not isinstance(entries, list) or len(entries) != 2:
        problems.append(
            f"expected exactly 2 census entries, got "
            f"{len(entries) if isinstance(entries, list) else type(entries).__name__}")
        entries = entries if isinstance(entries, list) else []

    per: dict = {}
    tri: dict = {}
    diags = [e.get("diagonal") for e in entries if isinstance(e, dict)]
    if sorted(d for d in diags if _int(d) is not None) != [0, 1] \
            or len(set(diags)) != len(diags):
        problems.append(f"diagonals are {diags!r}; require exactly one d0 "
                        "and one d1")
    for e in entries:
        if not isinstance(e, dict):
            problems.append(f"census entry is {type(e).__name__}, not an object")
            continue
        d = e.get("diagonal")
        t, n = _int(e.get("transverse")), _int(e.get("triangles"))
        if t is None or t < 0:
            problems.append(f"d{d}: transverse is {e.get('transverse')!r}; "
                            "require a non-negative integer")
        if n is None or n <= 0:
            # Zero triangles means the loader saw no surface. A clean
            # verdict over nothing is vacuous, and is how this project's
            # own published defect stayed invisible.
            problems.append(f"d{d}: triangles is {e.get('triangles')!r}; "
                            "require a positive integer -- a census over "
                            "zero triangles is vacuous, not clean")
        per[f"d{d}"] = t
        tri[f"d{d}"] = n
    for k, v in CENSUS.items():
        params = r.get("parameters") or {}
        if k not in params:
            problems.append(f"parameter {k} absent from the report")
        else:
            try:
                if abs(float(params[k]) - v) > 1e-12:
                    problems.append(f"parameter {k} is {params[k]!r}, frozen "
                                    f"value is {v}")
            except (TypeError, ValueError):
                problems.append(f"parameter {k} is {params[k]!r}, not numeric")
    for gk in ("grid_rows", "grid_cols"):
        if _int(r.get(gk)) is None or _int(r.get(gk)) <= 0:
            problems.append(f"{gk} is {r.get(gk)!r}")

    return {"ran": True, "transverse": per, "triangles": tri,
            "clean": (not problems
                      and all(v == 0 for v in per.values())),
            "grid": [r.get("grid_rows"), r.get("grid_cols")],
            "parameters_ok": not problems,
            "parameter_problems": problems}


def census(surface: Path, runner: list[str]) -> dict:
    # The report is written BESIDE the surface, not into the system temp
    # directory. A containerised validator only sees the paths its
    # launcher mounts, and on macOS $TMPDIR is not one of them: the
    # census then succeeds while the report lands nowhere the caller can
    # read, which reads as "could not be censused".
    # ABSOLUTE paths, always. The validator is frequently a container
    # launcher whose working directory is not the caller's, so a relative
    # path it is handed resolves to nothing and the census fails with
    # "cannot open meta.json" -- which reads as a broken surface rather
    # than a broken invocation. Found by cold-running the published
    # package from a fresh extraction.
    surface = surface.resolve()
    with tempfile.TemporaryDirectory(prefix=".topoeval-",
                                     dir=surface.parent) as td:
        rep = Path(td).resolve() / "r.json"
        cmd = runner + [str(surface), "-o", str(rep)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0 or not rep.is_file():
            return {"ran": False,
                    "detail": (proc.stderr or proc.stdout)[-400:]}
        out = read_census_report(json.loads(rep.read_text()))
        out["source"] = "run by this tool on the candidate"
        return out


def plane_hashes(d: Path) -> dict:
    out = {}
    for name in ("x.tif", "y.tif", "z.tif", "mask.tif", "meta.json"):
        p = d / name
        out[name] = sha256(p) if p.is_file() else None
    return out


def check_census_binding(binding: Path, report: Path,
                         candidate: Path) -> list[str]:
    """A supplied census must be tied to the candidate it is supposed to
    describe, or it certifies nothing.

    Without this, a caller hands the evaluator a dirty candidate and a
    clean report belonging to any other surface of the same shape, and the
    evaluator reports the dirty candidate clean. Grid and parameter checks
    do not catch it: the two surfaces agree on both.

    A binding cannot make a self-reported census trustworthy -- whoever
    produced the report could produce the binding too. What it does is
    make an ACCIDENTAL or CASUAL swap impossible, and it is why a census
    this tool did not run is never marked authoritative.
    """
    problems: list[str] = []
    try:
        b = json.loads(Path(binding).read_text())
    except Exception as exc:                       # noqa: BLE001
        return [f"census binding unreadable: {exc}"]
    if sha256(Path(report)) != b.get("report_sha256"):
        problems.append("the census report's sha256 is not the one the "
                        "binding records")
    want = b.get("operand_sha256") or {}
    got = plane_hashes(Path(candidate))
    if not want:
        problems.append("the binding records no operand hashes")
    for k, v in want.items():
        if got.get(k) != v:
            problems.append(f"candidate {k} is {str(got.get(k))[:12]}, the "
                            f"binding says {str(v)[:12]}")
    if not b.get("validator"):
        problems.append("the binding records no validator identity")
    return problems


def evaluate(input_dir: Path, candidate: Path, runner: list[str],
             reference: Path | None = None,
             census_report: Path | None = None,
             census_binding: Path | None = None,
             provenance: Path | None = None) -> dict:
    xi, yi, zi, vi = load_surface(input_dir)
    xc, yc, zc, vc = load_surface(candidate)
    out: dict = {"input": str(input_dir), "candidate": str(candidate),
                 "census_parameters": CENSUS}

    # ---- same input, or refuse (B5) ----------------------------------
    if vc.shape != vi.shape:
        out["comparable"] = False
        out["refused"] = (f"candidate grid {list(vc.shape)} != input grid "
                          f"{list(vi.shape)}")
        return out
    out["comparable"] = True
    out["grid"] = list(vi.shape)

    # The benchmark's rule is "same input, BY HASH". Shape agreement is
    # not identity: two different traces can share a grid. When the
    # trace's provenance.json is supplied the input is verified against
    # it; when it is not, the report says so instead of implying a check
    # that never ran.
    if provenance is not None:
        p = json.loads(Path(provenance).read_text())
        want = ((p.get("input") or {}).get("hashes") or {})
        got = {a: sha256(input_dir / f"{a}.tif") for a in "xyz"
               if (input_dir / f"{a}.tif").is_file()}
        wrong = {a: [got.get(a), want.get(a)] for a in want
                 if want.get(a) and got.get(a) != want.get(a)}
        out["input_identity"] = {
            "verified": not wrong and bool(want),
            "against": str(provenance),
            "mismatches": wrong}
        if wrong:
            out["refused"] = (f"the input is not the benchmark's input: "
                              f"{sorted(wrong)} differ by hash")
            return out
        if not want:
            out["input_identity"]["note"] = (
                "the provenance file records no input hashes")
    else:
        out["input_identity"] = {
            "verified": False,
            "note": ("no provenance file supplied: the input was matched "
                     "by GRID SHAPE ONLY, which is not identity. Pass "
                     "--provenance <trace>/provenance.json for a "
                     "benchmark-grade result.")}

    # ---- the invariant, primary and binary (B4.1) ---------------------
    if census_report is not None:
        # A caller may supply a census it already ran -- a corpus-wide
        # sweep does one container for hundreds of surfaces instead of one
        # each. The report is NOT taken on trust: its grid must match the
        # candidate, its parameters must be the frozen ones, and the
        # result records that this tool did not run it.
        cen = read_census_report(json.loads(Path(census_report).read_text()))
        cen["source"] = (f"PRECOMPUTED report supplied by the caller "
                         f"({census_report}); NOT run by this tool")
        cen["authoritative"] = False
        if cen["grid"] != list(vc.shape):
            out["census"] = cen
            out["refused"] = (f"the supplied census report is for grid "
                              f"{cen['grid']}, the candidate grid is "
                              f"{list(vc.shape)}")
            return out
        if census_binding is None:
            out["census"] = cen
            out["refused"] = (
                "a census this tool did not run must come with "
                "--census-binding tying the report to THIS candidate's "
                "bytes. Without it, a clean report from any other surface "
                "of the same shape would certify this one.")
            return out
        bad = check_census_binding(census_binding, census_report, candidate)
        cen["binding"] = {"file": str(census_binding), "problems": bad,
                          "ok": not bad}
        if bad:
            out["census"] = cen
            out["refused"] = f"census binding rejected: {bad[0]}"
            return out
    else:
        cen = census(candidate, runner)
        cen["authoritative"] = True
    out["census"] = cen
    if not cen.get("ran"):
        out["refused"] = "the candidate could not be censused"
        return out
    if not cen.get("parameters_ok"):
        out["refused"] = f"census parameters differ: {cen['parameter_problems']}"
        return out
    out["clean"] = bool(cen["clean"])

    ai, oki = quad_areas(xi, yi, zi, vi)
    ac, okc = quad_areas(xc, yc, zc, vc)
    A_in = float(ai.sum())
    out["input_valid_cells"] = int(vi.sum())
    out["input_area_canonical"] = A_in

    if not out["clean"]:
        # B4.1: cost is not reported for a candidate that never satisfied
        # the invariant. Ranking it on retention would put "do nothing"
        # first.
        out["scored"] = False
        out["verdict"] = ("FAILS THE INVARIANT: still self-intersecting "
                          f"({cen['transverse']}); cost metrics are "
                          "deliberately not reported")
        return out

    out["scored"] = True

    # ---- cost (B4.2-4) -------------------------------------------------
    retained = vi & vc
    removed = vi & ~vc
    added = ~vi & vc
    out["removed_cells"] = int(removed.sum())
    out["removed_cell_fraction"] = (float(removed.sum() / vi.sum())
                                    if vi.sum() else None)
    out["added_cells"] = int(added.sum())

    # ---- retention, priced on the INPUT (B4.2) --------------------------
    # The comparable number is how much of the INPUT's surface the
    # candidate still covers, priced with the INPUT's own quad areas.
    #
    # Pricing the candidate's own geometry instead lets a method score
    # above 1.0 by inventing surface: add cells outside the input, or move
    # retained coordinates apart, and "retained_fraction" rises while less
    # of the original is actually kept. Both are reported separately below
    # so a method that legitimately does something other than excision is
    # visible rather than penalised silently.
    covered = oki & okc
    out["retained_fraction"] = (float(ai[covered].sum() / A_in)
                                if A_in else None)
    out["retained_fraction_definition"] = (
        "input quad area covered by the candidate, divided by input quad "
        "area. Priced on the INPUT, so geometry the candidate adds or "
        "stretches cannot inflate it.")
    out["candidate_geometric_area_fraction"] = (
        float(ac.sum() / A_in) if A_in else None)
    out["candidate_geometric_area_note"] = (
        "the candidate's OWN area over the input's. Equal to "
        "retained_fraction for a pure excision; above it means the "
        "candidate added or enlarged geometry, which is not excision and "
        "is not comparable like-for-like.")

    # ---- coordinate fidelity (B4.5) ------------------------------------
    same = (np.array_equal(xi[retained], xc[retained])
            and np.array_equal(yi[retained], yc[retained])
            and np.array_equal(zi[retained], zc[retained]))
    out["coordinate_fidelity"] = bool(same)
    if not same:
        moved = int((~np.isclose(xi[retained], xc[retained])
                     | ~np.isclose(yi[retained], yc[retained])
                     | ~np.isclose(zi[retained], zc[retained])).sum())
        out["moved_cells"] = moved
        out["coordinate_note"] = (
            "retained cells are NOT bit-identical to the input: this "
            "candidate MOVES geometry rather than only removing it. That "
            "is not disqualifying, but its retained_fraction is not "
            "comparable like-for-like with an excision method, and it is "
            "not evidence about the scroll.")

    # ---- fragmentation (B4.4) ------------------------------------------
    lab_i, n_i = components(oki)
    lab_c, n_c = components(okc)
    area_in = np.bincount(lab_i[oki], weights=ai[oki], minlength=n_i + 1)
    area_c = np.bincount(lab_c[okc], weights=ac[okc], minlength=n_c + 1)
    # A component of exclusively degenerate quads has zero canonical area
    # and is not surface. The published certificates count only
    # positive-area components; counting the rest would report
    # fragmentation that no area measurement can see.
    out["components_input"] = int((area_in[1:] > 0).sum())
    out["components_candidate"] = int((area_c[1:] > 0).sum())
    out["components_input_including_zero_area"] = int(n_i)
    out["components_candidate_including_zero_area"] = int(n_c)
    # The parent map below assumes every candidate quad is also an input
    # quad, so each candidate component lies inside exactly one input
    # component. A candidate that ADDS quads breaks that: one of its
    # components can bridge two input components, and R_main would then be
    # computed against an arbitrary parent. Refuse rather than report a
    # number whose meaning has quietly changed.
    outside = okc & ~oki
    out["candidate_quads_outside_input"] = int(outside.sum())
    r_main = {}
    if outside.any():
        out["fragmentation_comparable"] = False
        out["fragmentation_note"] = (
            f"{int(outside.sum())} candidate quads are not input quads, so "
            "a candidate component can span more than one input component "
            "and R_main has no well-defined parent. Fragmentation is not "
            "reported for this candidate.")
        out["R_main_min"] = None
        out["R_main_area_weighted"] = None
    elif n_i and n_c:
        out["fragmentation_comparable"] = True
        pi = lab_i[okc]                       # parent input component
        po = lab_c[okc]
        area_by_out = area_c
        parent = np.zeros(n_c + 1, np.int64)
        parent[po] = pi
        largest = np.zeros(n_i + 1)
        for oc in range(1, n_c + 1):
            p = int(parent[oc])
            largest[p] = max(largest[p], area_by_out[oc])
        for ic in range(1, n_i + 1):
            if area_in[ic] > 0:
                r_main[ic] = float(largest[ic] / area_in[ic])
    out.setdefault("fragmentation_comparable", bool(r_main))
    out["R_main_min"] = min(r_main.values()) if r_main else None
    out["R_main_area_weighted"] = (
        float(sum(r_main[i] * area_in[i] for i in r_main)
              / sum(area_in[i] for i in r_main)) if r_main else None)

    # ---- the 99.9%-area core gate, same definition as the certificates --
    # Without this the evaluator reported only the aggregate R_main, and a
    # method could destroy a small-but-real component while the weighted
    # average stayed near 1. The core gate is the check that notices.
    core = {"core_fraction_target": 0.999, "gate_R_main": 0.9}
    if r_main:
        live = sorted(r_main, key=lambda i: (-area_in[i], i))
        total = float(sum(area_in[i] for i in live))
        cum, core_ids = 0.0, []
        for i in live:
            if total > 0 and cum >= 0.999 * total:
                break
            core_ids.append(i)
            cum += float(area_in[i])
        core_r = [r_main[i] for i in core_ids]
        core.update(
            n_input_components=len(live),
            n_core_components=len(core_ids),
            n_tail_components=len(live) - len(core_ids),
            core_area_fraction=(cum / total) if total else None,
            min_R_main_core=(min(core_r) if core_r else None),
            n_core_components_below_gate=sum(1 for v in core_r if v < 0.9),
            core_gate_pass=bool(core_r and min(core_r) >= 0.9) or not core_r)
    else:
        core.update(n_input_components=None, n_core_components=None,
                    n_tail_components=None, core_area_fraction=None,
                    min_R_main_core=None, n_core_components_below_gate=None,
                    core_gate_pass=None)
    out["core_gate"] = core

    # ---- the loader divergence, reported rather than reconciled --------
    out["n_retained_cells_at_z_le_0"] = int((vc & (zc <= 0)).sum())
    out["loader_semantics"] = (
        "CERTIFICATE semantics: mask >= 255, finite coordinates, not the "
        "-1 marker. The official loader ALSO invalidates z <= 0, so where "
        "n_retained_cells_at_z_le_0 is non-zero this measurement covers a "
        "superset of what the official validator examines.")

    # ---- comparison against the reference, if given (B5) --------------
    if reference is not None and reference.is_dir():
        xr, yr, zr, vr = load_surface(reference)
        ar, okr = quad_areas(xr, yr, zr, vr)
        ref_frac = float(ar.sum() / A_in) if A_in else None
        out["reference"] = {
            "retained_fraction": ref_frac,
            "removed_cells": int((vi & ~vr).sum()),
            "note": ("the reference is ONE admissible answer, not the "
                     "target: a clean candidate retaining more area is "
                     "better, and this comparison says so"),
        }
        if ref_frac is not None and out["retained_fraction"] is not None:
            delta = out["retained_fraction"] - ref_frac
            out["reference"]["retention_delta"] = delta
            out["reference"]["verdict"] = (
                "retains MORE than the reference" if delta > 1e-9 else
                "retains LESS than the reference" if delta < -1e-9 else
                "retains the same area as the reference")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Score a candidate derivative against a benchmark trace")
    ap.add_argument("--input", required=True, help="the input tifxyz")
    ap.add_argument("--candidate", required=True,
                    help="the derivative to score")
    ap.add_argument("--reference", default=None,
                    help="the benchmark's own derivative, for comparison")
    ap.add_argument("--validator", default="",
                    help="vc_tifxyz_selfcross, or a launcher for it "
                         "(shell-split). Not needed with --census-report.")
    ap.add_argument("--census-report", default=None,
                    help="a census of the CANDIDATE this tool did not run. "
                         "Requires --census-binding, and the result is "
                         "never marked authoritative.")
    ap.add_argument("--census-binding", default=None,
                    help="manifest tying --census-report to this "
                         "candidate's bytes and to a validator identity")
    ap.add_argument("--provenance", default=None,
                    help="<trace>/provenance.json; verifies that --input "
                         "IS the benchmark's input, by hash")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    if not a.validator and not a.census_report:
        ap.error("give --validator, or --census-report for a census you "
                 "already ran")
    import shlex
    r = evaluate(Path(a.input), Path(a.candidate), shlex.split(a.validator),
                 Path(a.reference) if a.reference else None,
                 Path(a.census_report) if a.census_report else None,
                 Path(a.census_binding) if a.census_binding else None,
                 Path(a.provenance) if a.provenance else None)
    if a.json:
        Path(a.json).write_text(json.dumps(r, indent=1) + "\n")

    if not r.get("comparable"):
        print(f"REFUSED: {r.get('refused')}")
        return 2
    if r.get("refused"):
        print(f"REFUSED: {r['refused']}")
        return 2
    if not r.get("scored"):
        print(r["verdict"])
        return 3
    print(f"clean                : yes ({r['census']['transverse']})")
    print(f"  census             : {r['census'].get('source')}")
    if not r["census"].get("authoritative"):
        print("  NOT AUTHORITATIVE  : this tool did not run the census. A "
              "benchmark-grade result requires --validator.")
    ii = r.get("input_identity") or {}
    print(f"input identity       : "
          f"{'verified by hash' if ii.get('verified') else 'NOT VERIFIED'}")
    if ii.get("note"):
        print(f"  {ii['note']}")
    print(f"retained_fraction    : {r['retained_fraction']:.6f}  "
          f"(input-priced)")
    if r.get("candidate_geometric_area_fraction") is not None and \
            abs(r["candidate_geometric_area_fraction"]
                - r["retained_fraction"]) > 1e-12:
        print(f"  candidate own area : "
              f"{r['candidate_geometric_area_fraction']:.6f} -- differs, so "
              f"this candidate is not a pure excision")
    print(f"removed_cells        : {r['removed_cells']} "
          f"({(r['removed_cell_fraction'] or 0) * 100:.4f}% of valid)")
    print(f"components           : {r['components_input']} -> "
          f"{r['components_candidate']}")
    if r.get("R_main_min") is not None:
        print(f"R_main min / weighted: {r['R_main_min']:.6f} / "
              f"{r['R_main_area_weighted']:.6f}")
        cg = r.get("core_gate") or {}
        print(f"99.9%-area core gate : "
              f"{'pass' if cg.get('core_gate_pass') else 'FAIL'} "
              f"({cg.get('n_core_components')} core / "
              f"{cg.get('n_tail_components')} tail, min R_main "
              f"{cg.get('min_R_main_core')})")
    elif r.get("fragmentation_comparable") is False:
        print(f"fragmentation        : NOT COMPARABLE -- "
              f"{r.get('fragmentation_note')}")
    print(f"coordinate_fidelity  : "
          f"{'bit-identical' if r['coordinate_fidelity'] else 'MOVED'}")
    if not r["coordinate_fidelity"]:
        print(f"  {r['coordinate_note']}")
    if "reference" in r:
        print(f"vs reference         : {r['reference']['verdict']} "
              f"({r['reference']['retention_delta']:+.6f})")
    print("\nNo single score is emitted by design: retention, "
          "fragmentation and fidelity trade against each other.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
