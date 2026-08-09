"""Decide whether a quickstart demonstration actually demonstrated anything.

The quickstart runs the evaluator twice: once on the ORIGINAL trace,
which must be rejected, and once on the reconstructed REFERENCE, which
must be clean and scored. The first of those is the dangerous one. Its
process exit status is non-zero *because it worked* -- and non-zero is
also what a dead validator, an unmounted path, a malformed report, a
stale file from a previous run and a census over an empty surface all
produce. Reading "non-zero" as the success condition made every one of
those failures print as a passing demonstration.

So the demonstration is judged on its REPORT, not its exit status alone,
and the report must satisfy every one of:

  the census RAN, and this tool ran it (a precomputed census is not
  authoritative here); the frozen parameters were used; both diagonals
  carry a POSITIVE triangle count, because a clean verdict over an empty
  surface is vacuous; the report describes THE OPERANDS WE INTENDED, not
  whatever a previous invocation left behind; the input was verified by
  hash against provenance.json; the verdict matches the expected
  polarity; and the process exit status agrees with the verdict.

For the clean case there is a further binding. "Some clean surface"
proves nothing about this package: a corrupted `masks.npz` that removed
far more geometry would still be clean. The reconstruction is therefore
required to BE the packaged reference -- its coordinate planes must hash
to the values `provenance.json` records for the reference derivative --
and its cost vector must match the certificate figures the package
ships in `fields.json`.

This lives in its own module rather than inside the shell script so the
failure modes above can be regression-tested directly, by feeding it
reports instead of by arranging for a validator to break.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")

# The same relative tolerance the acceptance test uses to compare the
# evaluator against a certificate. Reusing it keeps one definition of
# "agrees" in the package.
RETAINED_FRACTION_RTOL = 1e-9

# What an inaccessible path looks like coming out of a containerised
# validator. Matching on this decides whether the mount advice is
# relevant -- printing it after every failure trains the reader to
# ignore it, which is worse than not printing it.
_MOUNT_HINT = re.compile(
    r"cannot open|no such file|not a directory|meta\.json|permission denied",
    re.IGNORECASE)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def _same_path(a, b) -> bool:
    try:
        return Path(str(a)).resolve() == Path(str(b)).resolve()
    except OSError:
        return False


def check(report: Path, expect: str, exit_status: int,
          input_dir: Path, candidate: Path,
          provenance: Path | None = None,
          fields: Path | None = None) -> list[str]:
    """Return the reasons this demonstration is not evidence. Empty is a
    pass."""
    if expect not in ("dirty", "clean"):
        return [f"unknown expectation {expect!r}"]

    if not report.is_file():
        return [f"the evaluator wrote no report to {report}"]
    try:
        d = json.loads(report.read_text())
    except Exception as exc:                                  # noqa: BLE001
        return [f"the report is not readable JSON: {exc}"]

    cen = d.get("census") or {}
    ii = d.get("input_identity") or {}

    # Reported alone. Everything below is derived from a census that did
    # not happen, so listing ten consequences of one cause buries it.
    if cen.get("ran") is not True:
        return [f"the census did not run: {cen.get('detail')}"]
    if d.get("comparable") is not True:
        return [f"candidate not comparable to the input: {d.get('refused')}"]
    if d.get("refused"):
        return [f"the evaluator refused: {d['refused']}"]

    bad: list[str] = []

    # The report must be about the operands we asked for. Without this a
    # stale file from an earlier run -- or a report for another surface
    # entirely -- is read as this run's result.
    if not _same_path(d.get("input"), input_dir):
        bad.append(f"the report is for input {d.get('input')!r}, not "
                   f"{str(input_dir)!r}")
    if not _same_path(d.get("candidate"), candidate):
        bad.append(f"the report is for candidate {d.get('candidate')!r}, not "
                   f"{str(candidate)!r}")

    if ii.get("verified") is not True:
        bad.append("the input was NOT verified by hash against "
                   "provenance.json, so this is not the benchmark's trace")
    if cen.get("authoritative") is not True:
        bad.append("the census was not run by the evaluator, so it is not "
                   "authoritative")
    if cen.get("parameters_ok") is not True:
        bad.append(f"census parameters/report rejected: "
                   f"{cen.get('parameter_problems')}")

    tri = cen.get("triangles") or {}
    tr = cen.get("transverse") or {}
    for k in ("d0", "d1"):
        v = tri.get(k)
        if not isinstance(v, int) or isinstance(v, bool) or v <= 0:
            bad.append(f"{k}: {v!r} triangles -- a verdict over an empty "
                       "surface is vacuous, not a result")
        t = tr.get(k)
        if not isinstance(t, int) or isinstance(t, bool) or t < 0:
            bad.append(f"{k}: transverse count is {t!r}")
    live = [v for v in tr.values()
            if isinstance(v, int) and not isinstance(v, bool)]

    if expect == "dirty":
        if d.get("clean") is not False:
            bad.append(f"expected an unclean verdict, got "
                       f"clean={d.get('clean')!r}")
        if d.get("scored") is not False:
            bad.append("an unclean candidate must NOT be scored on cost, "
                       f"and scored={d.get('scored')!r}")
        if len(live) == 2 and not any(v > 0 for v in live):
            bad.append("no transverse crossings were found, so this trace "
                       "does not demonstrate rejection")
        if exit_status != 3:
            bad.append(f"exit status {exit_status}; an unclean candidate "
                       "must exit 3 (2 is a refusal, 0 is a pass)")
        return bad

    # ---- clean ---------------------------------------------------------
    if d.get("clean") is not True:
        bad.append(f"expected a clean verdict, got clean={d.get('clean')!r}")
    if d.get("scored") is not True:
        bad.append("the reference must be scored on cost")
    if any(v != 0 for v in live):
        bad.append(f"transverse crossings remain: {tr}")
    if exit_status != 0:
        bad.append(f"exit status {exit_status}; a clean scored candidate "
                   "exits 0")

    # The bindings below are what make this "the packaged reference"
    # rather than "a clean surface". Without them the clean check is a
    # different, weaker claim than this module says it makes, so their
    # absence is a failure rather than a reduced mode.
    if provenance is None:
        bad.append("no provenance supplied: the reconstruction cannot be "
                   "bound to the packaged reference derivative")
    if fields is None:
        bad.append("no fields.json supplied: the cost vector cannot be "
                   "compared to the certificate figures")

    rf = d.get("retained_fraction")
    # NaN is a float, and every comparison against it is false -- so a
    # NaN retention would satisfy both "is a number" and any tolerance
    # test written as `abs(a - b) > tol`.
    if not isinstance(rf, float) or not math.isfinite(rf):
        bad.append(f"retained_fraction is {rf!r}, not a finite number")
        rf = None
    elif not 0.0 <= rf <= 1.0:
        bad.append(f"retained_fraction is {rf!r}, outside [0, 1]")
        rf = None
    if d.get("coordinate_fidelity") is not True:
        bad.append("the reconstruction MOVED geometry; the reference is a "
                   "pure excision and its retained cells are bit-identical "
                   "to the input")
    if d.get("added_cells") != 0:
        bad.append(f"the reconstruction added {d.get('added_cells')!r} cells; "
                   "the reference removes geometry and adds none")
    if d.get("fragmentation_comparable") is not True:
        bad.append("fragmentation is not comparable, so this is not the "
                   "packaged reference")

    # The strongest binding: this IS the packaged reference, not merely
    # some clean derivative of the same input.
    if provenance is not None:
        try:
            p = json.loads(Path(provenance).read_text())
        except Exception as exc:                              # noqa: BLE001
            bad.append(f"provenance.json unreadable: {exc}")
        else:
            want = ((p.get("reference_derivative") or {}).get("hashes") or {})
            # ALL THREE, or the binding does not hold. Skipping an axis
            # whose hash is absent would let a provenance record carrying
            # only `x` satisfy a check this module describes as x/y/z.
            for a in ("x", "y", "z"):
                h = want.get(a)
                if not isinstance(h, str) or not _SHA256.match(h):
                    bad.append(f"provenance.json records no valid {a} hash "
                               "for the reference derivative, so the "
                               "reconstruction cannot be bound")
                    continue
                f = Path(candidate) / f"{a}.tif"
                if not f.is_file():
                    bad.append(f"the reconstruction has no {a}.tif")
                elif sha256(f) != h:
                    bad.append(
                        f"the reconstructed {a}.tif does not hash to the "
                        "packaged reference derivative: this is a clean "
                        "surface, but it is not THE reference")

    if fields is not None and rf is not None:
        try:
            fj = json.loads(Path(fields).read_text())
        except Exception as exc:                              # noqa: BLE001
            bad.append(f"fields.json unreadable: {exc}")
        else:
            ref = (((fj.get("retention") or {}).get("area") or {})
                   .get("canonical") or {}).get("retained_fraction")
            if not isinstance(ref, float) or not math.isfinite(ref) \
                    or not 0.0 <= ref <= 1.0:
                # Range-checked explicitly. It cannot currently fall
                # outside [0,1] and still match a range-constrained report
                # value, but a contract the code only implies is a
                # contract the next edit can drop.
                bad.append("fields.json records no finite canonical retained "
                           "fraction in [0, 1] to compare against")
            elif abs(rf - ref) > RETAINED_FRACTION_RTOL * max(abs(ref), 1.0):
                bad.append(f"retained_fraction {rf!r} does not match the "
                           f"certificate figure {ref!r} shipped in "
                           "fields.json")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--report", required=True)
    ap.add_argument("--expect", required=True, choices=("dirty", "clean"))
    ap.add_argument("--exit-status", required=True, type=int)
    ap.add_argument("--input", required=True)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--provenance", default=None)
    ap.add_argument("--fields", default=None)
    # Only used to make the mount advice concrete, and only printed when
    # the recorded failure actually looks like an inaccessible path.
    ap.add_argument("--package", default=None)
    ap.add_argument("--workdir", default=None)
    a = ap.parse_args()

    report = Path(a.report)
    bad = check(report, a.expect, a.exit_status, Path(a.input),
                Path(a.candidate),
                Path(a.provenance) if a.provenance else None,
                Path(a.fields) if a.fields else None)
    if not bad:
        print(f"    gate: the {a.expect} demonstration behaved exactly as "
              "specified")
        return 0

    print(f"\n    GATE FAILED for the {a.expect.upper()} demonstration:")
    for b in bad:
        print(f"      - {b}")

    detail = ""
    if report.is_file():
        try:
            detail = str(((json.loads(report.read_text()).get("census") or {})
                          .get("detail")) or "")
        except Exception:                                     # noqa: BLE001
            detail = ""
    blob = detail + " " + " ".join(bad)
    if detail:
        print(f"\n    what the validator actually said:\n      {detail.strip()}")
    if _MOUNT_HINT.search(blob):
        print(f"""
    That is a path the validator could not open, which usually means the
    surface is fine and the VALIDATOR cannot SEE it: a containerised
    vc_tifxyz_selfcross only reads directories its launcher mounts.

    Either extract this package somewhere the launcher mounts, or set
    WORK=<a directory inside the mount> and re-run.

    Currently: package  {a.package}
               workdir  {a.workdir}""")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
