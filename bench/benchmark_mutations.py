"""Negative controls for the evaluator: every mutation MUST be caught.

A benchmark that only ever reports "agrees" is not evidence. The
acceptance test scores reference derivatives, which are all clean, all
pure excisions and all correctly paired -- so on its own it cannot show
that the evaluator would reject anything.

Each case below breaks exactly one property and states what must happen.
The suite fails if any mutation is accepted.

  dirty            the input itself, which self-intersects
  empty            a candidate whose mask retains nothing
  moved            a retained coordinate perturbed
  added            a removed cell reinstated, so the candidate has
                   geometry the reference does not
  wrong_input      scored against a DIFFERENT trace of the same grid
  swapped_census   a clean report belonging to another surface
  unbound_census   a report with no binding at all

`wrong_input` and `swapped_census` are the two that matter most: they are
how a sceptical entrant would obtain a pass without doing the work.

One case, `relative_paths`, is a POSITIVE control and not a mutation: it
requires the evaluator to SUCCEED. Reporting the suite as "8/8 caught"
therefore misdescribes it -- seven things must be rejected and one must be
accepted, and a reader counting rejections gets the wrong number. Each
case records its own `polarity`, and the summary counts the two
populations separately rather than leaving a consumer to infer them from
the prose in `must`.

Every case is also DECLARED before the run. A case whose census fails is
otherwise skipped silently and the suite prints "5/5 caught -- PASS": a
missing control reads exactly like a passing one. The declared list is
reconciled at the end and any absence is a FAIL.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parent))
from topology_eval import evaluate  # noqa: E402


def copy_surface(src: Path, dst: Path) -> Path:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    return dst


def census(surface: Path, validator: str, out: Path) -> Path | None:
    r = subprocess.run([validator, str(surface.resolve()), "-o",
                        str(out.resolve())],
                       capture_output=True, text=True)
    return out if (r.returncode == 0 and out.is_file()) else None


def bind(report: Path, surface: Path, out: Path) -> Path:
    subprocess.run([sys.executable,
                    str(Path(__file__).parent / "bind_census_reports.py"),
                    "--manifest", "/dev/stdin", "--reports", str(out.parent)],
                   input=f"{out.stem.replace('.binding','')}\t{surface.resolve()}\n",
                   text=True, capture_output=True, check=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", default="out/benchmark")
    ap.add_argument("--work", default="out/benchmark-evidence/mutations")
    ap.add_argument("--validator",
                    default=os.environ.get("VC_SELFCROSS",
                                           "vc_tifxyz_selfcross"),
                    help="vc_tifxyz_selfcross, or a launcher for it; "
                         "defaults to $VC_SELFCROSS")
    ap.add_argument("--indices", nargs="+",
                    default=["out/release/index.json",
                             "out/release/expansion_index.json"])
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    entries: dict = {}
    for p in a.indices:
        entries.update(json.loads(Path(p).read_text())["segments"])
    bench = json.loads((Path(a.benchmark) / "index.json").read_text())
    scored = [t for t in bench["traces"] if t.get("scored")]

    # A transformed trace: one with a real cut to mutate.
    cut = [t for t in scored if t.get("n_removed_cells")]
    if not cut:
        print("REFUSING: no packaged trace has a cut to mutate")
        return 2
    victim = max(cut, key=lambda t: t["n_removed_cells"])
    seg = victim["segment"]
    ev = entries[seg]
    inp = Path(ev["input_mesh"])
    ref = Path(ev["output_mesh"] or ev["input_mesh"])
    prov = Path(a.benchmark) / "traces" / seg / "provenance.json"

    # No two packaged traces share a grid, so the swap attacks are built
    # from SYNTHESISED same-shaped surfaces instead. That is the realistic
    # threat anyway: an attacker does not need a second published trace,
    # only a second array of the same shape.
    oseg = "(synthesised same-grid decoy)"

    work = Path(a.work)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    print(f"victim  {seg}  grid {victim['grid']}\nother   {oseg}\n")

    # The honest baseline for this trace: reference + its own bound census.
    base_rep = census(ref, a.validator, work / "ref.json")
    if base_rep is None:
        print("REFUSING: the reference did not census")
        return 2
    bind(base_rep, ref, work / "ref.binding.json")

    # Declared before anything runs. A case that never executes -- because
    # its census failed, or because an edit dropped it -- must not simply
    # vanish from the denominator.
    DECLARED = {"dirty": "negative", "empty": "negative", "moved": "negative",
                "added": "negative", "wrong_input": "negative",
                "swapped_census": "negative", "relative_paths": "positive",
                "unbound_census": "negative"}

    def refused_because(res, *needles) -> bool:
        """The control must be caught for ITS OWN reason.

        Accepting any refusal at all means a globally broken evaluator --
        one that refuses everything it is handed -- passes four of the
        seven negative controls, and the suite reports PASS while proving
        nothing about the specific defence each case exists to test.
        """
        r = res.get("refused")
        if not isinstance(r, str):
            return False
        low = r.lower()
        return any(n.lower() in low for n in needles)

    cases: list[dict] = []

    def record(name, must, res, ok):
        polarity = DECLARED.get(name)
        if polarity is None:
            raise SystemExit(f"REFUSING: case {name!r} was not declared, so "
                             "the suite cannot say whether it is a mutation "
                             "to reject or a control to accept")
        cases.append({"case": name, "polarity": polarity,
                      "must": must, "caught": bool(ok),
                      "refused": res.get("refused"),
                      "clean": res.get("clean"),
                      "scored": res.get("scored"),
                      "retained_fraction": res.get("retained_fraction"),
                      "coordinate_fidelity": res.get("coordinate_fidelity"),
                      "added_cells": res.get("added_cells"),
                      "fragmentation_comparable":
                          res.get("fragmentation_comparable")})
        print(f"{'OK  ' if ok else 'MISS'}  {name:16s} {must}")

    # 1. dirty: the input self-intersects, so it must not be scored on cost
    rep = census(inp, a.validator, work / "dirty.json")
    if rep:
        bind(rep, inp, work / "dirty.binding.json")
        r = evaluate(inp, inp, [], census_report=rep,
                     census_binding=work / "dirty.binding.json",
                     provenance=prov)
        record("dirty", "not clean, and NOT scored on cost", r,
               r.get("clean") is False and not r.get("scored"))

    # 2. empty: retains nothing. A clean verdict here would be vacuous.
    empty = copy_surface(ref, work / "empty.tifxyz")
    m = np.asarray(tifffile.imread(empty / "mask.tif"))
    tifffile.imwrite(empty / "mask.tif", np.zeros_like(m))
    rep = census(empty, a.validator, work / "empty.json")
    if rep:
        bind(rep, empty, work / "empty.binding.json")
        r = evaluate(inp, empty, [], census_report=rep,
                     census_binding=work / "empty.binding.json",
                     provenance=prov)
        record("empty", "REFUSED as a census over zero triangles", r,
               refused_because(r, "zero triangles", "vacuous"))

    # 3. moved: a retained coordinate is perturbed
    moved = copy_surface(ref, work / "moved.tifxyz")
    x = np.asarray(tifffile.imread(moved / "x.tif"))
    mm = np.asarray(tifffile.imread(moved / "mask.tif"))
    idx = np.argwhere(mm >= 255)
    v, u = idx[len(idx) // 2]
    x[v, u] = np.float32(x[v, u] + 5.0)
    tifffile.imwrite(moved / "x.tif", x)
    rep = census(moved, a.validator, work / "moved.json")
    if rep:
        bind(rep, moved, work / "moved.binding.json")
        r = evaluate(inp, moved, [], census_report=rep,
                     census_binding=work / "moved.binding.json",
                     provenance=prov)
        record("moved", "coordinate_fidelity FALSE", r,
               r.get("coordinate_fidelity") is False)

    # 4. added: reinstate a removed cell -- geometry the reference lacks
    added = copy_surface(ref, work / "added.tifxyz")
    masks = np.load(Path(a.benchmark) / "traces" / seg / "masks.npz")
    rm = masks["removed"].astype(bool)
    am = np.asarray(tifffile.imread(added / "mask.tif"))
    xi = np.asarray(tifffile.imread(inp / "x.tif"))
    yi = np.asarray(tifffile.imread(inp / "y.tif"))
    zi = np.asarray(tifffile.imread(inp / "z.tif"))
    am[rm] = 255
    tifffile.imwrite(added / "mask.tif", am)
    for nm, arr in (("x.tif", xi), ("y.tif", yi), ("z.tif", zi)):
        cur = np.asarray(tifffile.imread(added / nm))
        cur[rm] = arr[rm]
        tifffile.imwrite(added / nm, cur)
    rep = census(added, a.validator, work / "added.json")
    if rep:
        bind(rep, added, work / "added.binding.json")
        r = evaluate(inp, added, [], census_report=rep,
                     census_binding=work / "added.binding.json",
                     provenance=prov)
        # Re-adding the excised cells restores the crossings, so this is
        # caught as not-clean. What must NEVER happen is a scored result
        # whose retention exceeds 1.
        # Every condition, not the weakest one. The old predicate was
        # satisfied by `retained_fraction <= 1`, which is true of a clean
        # SCORED result at 0.9 -- so the control passed without the
        # crossings having come back at all.
        record("added",
               "not clean (the crossings come back), and NOT scored, so no "
               "retention is reported at all",
               r,
               r.get("clean") is False
               and r.get("scored") is not True
               and r.get("retained_fraction") is None)

    # 5. wrong input: a same-grid input that is NOT this trace's
    wrong = copy_surface(inp, work / "wrong_input.tifxyz")
    wx = np.asarray(tifffile.imread(wrong / "x.tif"))
    wx[0, 0] = np.float32(wx[0, 0] + 1.0)
    tifffile.imwrite(wrong / "x.tif", wx)
    r = evaluate(wrong, ref, [], census_report=base_rep,
                 census_binding=work / "ref.binding.json", provenance=prov)
    record("wrong_input", "REFUSED: the input is not this trace's, by hash",
           r, refused_because(r, "not the benchmark's input", "differ by hash"))

    # 6. swapped census: a clean report belonging to another surface of
    # the same shape. The decoy is a copy of the reference with one
    # coordinate nudged -- still clean, still the same grid, different
    # bytes. This is the attack the binding exists to stop.
    oref = copy_surface(ref, work / "decoy.tifxyz")
    dx = np.asarray(tifffile.imread(oref / "x.tif"))
    dm = np.asarray(tifffile.imread(oref / "mask.tif"))
    di = np.argwhere(dm >= 255)
    dv, du = di[len(di) // 3]
    dx[dv, du] = np.float32(dx[dv, du] + 0.01)
    tifffile.imwrite(oref / "x.tif", dx)
    orep = census(oref, a.validator, work / "other.json")
    if orep:
        bind(orep, oref, work / "other.binding.json")
        r = evaluate(inp, ref, [], census_report=orep,
                     census_binding=work / "other.binding.json",
                     provenance=prov)
        record("swapped_census",
               "REFUSED: the binding is for another surface's bytes", r,
               refused_because(r, "census binding rejected"))

    # 6b. POSITIVE CONTROL, not a mutation: the reference scored through a
    # RELATIVE candidate path must succeed. The published v0.3.0 evaluator
    # passed the path to the validator verbatim, so a container launcher
    # resolved it against its own working directory and every relative
    # invocation failed with "cannot open meta.json" -- which reads as a
    # broken surface, not a broken call. Cold-running the published
    # package found it; this keeps it found.
    import os as _os
    rel_in = _os.path.relpath(inp, _os.getcwd())
    rel_ref = _os.path.relpath(ref, _os.getcwd())
    r = evaluate(Path(rel_in), Path(rel_ref), [a.validator], provenance=prov)
    record("relative_paths",
           "POSITIVE CONTROL: a relative candidate path must still census",
           r, bool(r.get("clean")) and bool(r.get("scored")))

    # 7. unbound census: no binding supplied at all
    r = evaluate(inp, ref, [], census_report=base_rep, provenance=prov)
    record("unbound_census", "REFUSED: a census we did not run needs a "
           "binding", r, refused_because(r, "--census-binding"))

    missed = [c for c in cases if not c["caught"]]
    # Every declared case must appear EXACTLY once. Set difference alone
    # catches an absent case but not a duplicated one, and a case recorded
    # twice would let one pass cover one failure.
    from collections import Counter
    seen = Counter(c["case"] for c in cases)
    absent = sorted(k for k in DECLARED if seen[k] == 0)
    repeated = sorted(k for k, n in seen.items() if n > 1)
    neg = [c for c in cases if c["polarity"] == "negative"]
    pos = [c for c in cases if c["polarity"] == "positive"]
    out = {"test": "evaluator controls",
           "victim": seg, "other_same_grid": oseg,
           "n_declared": len(DECLARED),
           "n_cases": len(cases),
           # Deliberately NOT "n_caught"/"n_missed": one of the eight is a
           # positive control that must be ACCEPTED, so a combined count
           # of things "caught" misdescribes the suite in the raw record
           # even when the console and card are right.
           "n_cases_satisfying_expected_outcome": len(cases) - len(missed),
           "n_cases_failing_expected_outcome": len(missed),
           "cases_not_run": absent,
           "cases_recorded_more_than_once": repeated,
           # Counted by polarity because they are not the same claim: seven
           # mutations must be REJECTED and one control must be ACCEPTED.
           # A combined "8/8 caught" reads as eight rejections.
           "n_negative": len(neg),
           "n_negative_rejected": sum(1 for c in neg if c["caught"]),
           "n_positive": len(pos),
           "n_positive_accepted": sum(1 for c in pos if c["caught"]),
           "cases": cases,
           "verdict": ("PASS" if not missed and not absent and not repeated
                       else "FAIL")}
    if a.out:
        Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"\n{out['n_negative_rejected']}/{out['n_negative']} mutations "
          f"rejected; {out['n_positive_accepted']}/{out['n_positive']} "
          f"positive control accepted -- {out['verdict']}")
    for c in missed:
        print(f"  {'MISSED' if c['polarity'] == 'negative' else 'FAILED'} "
              f"{c['case']}: {c['must']}")
    for name in absent:
        print(f"  NEVER RAN {name}: declared but produced no record")
    for name in repeated:
        print(f"  DUPLICATED {name}: recorded {seen[name]} times")
    return 0 if out["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
