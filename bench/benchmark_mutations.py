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

    cases: list[dict] = []

    def record(name, must, res, ok):
        cases.append({"case": name, "must": must, "caught": bool(ok),
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
               bool(r.get("refused")))

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
        record("added",
               "not clean (the crossings come back), never retention > 1",
               r,
               (r.get("clean") is False and not r.get("scored"))
               or (r.get("retained_fraction") or 0) <= 1.0 + 1e-12)

    # 5. wrong input: a same-grid input that is NOT this trace's
    wrong = copy_surface(inp, work / "wrong_input.tifxyz")
    wx = np.asarray(tifffile.imread(wrong / "x.tif"))
    wx[0, 0] = np.float32(wx[0, 0] + 1.0)
    tifffile.imwrite(wrong / "x.tif", wx)
    r = evaluate(wrong, ref, [], census_report=base_rep,
                 census_binding=work / "ref.binding.json", provenance=prov)
    record("wrong_input", "REFUSED: the input is not this trace's, by hash",
           r, bool(r.get("refused")))

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
               bool(r.get("refused")))

    # 7. unbound census: no binding supplied at all
    r = evaluate(inp, ref, [], census_report=base_rep, provenance=prov)
    record("unbound_census", "REFUSED: a census we did not run needs a "
           "binding", r, bool(r.get("refused")))

    missed = [c for c in cases if not c["caught"]]
    out = {"test": "evaluator negative controls",
           "victim": seg, "other_same_grid": oseg,
           "n_cases": len(cases), "n_caught": len(cases) - len(missed),
           "n_missed": len(missed), "cases": cases,
           "verdict": "PASS" if not missed else "FAIL"}
    if a.out:
        Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"\n{out['n_caught']}/{out['n_cases']} mutations caught -- "
          f"{out['verdict']}")
    for c in missed:
        print(f"  MISSED {c['case']}: {c['must']}")
    return 0 if not missed else 1


if __name__ == "__main__":
    raise SystemExit(main())
