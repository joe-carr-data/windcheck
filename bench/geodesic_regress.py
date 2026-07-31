"""Regression harness: engines/geodesic vs the golden Python spectra.

Zero regression is the objective. The contract is the golden files
out/spectrum_final_d{0,1}.json (185 segments per diagonal): for every
segment the C++ engine must reproduce every event's separation_mm (exact
equality of the 4-decimal rounded value), same_component, ambiguous,
self_touching, n_pairs, du_max, endpoint_exact, and the median
intersection length to 1e-9.

Events are matched by REGION SIGNATURE, not order: the golden files do
not store regions, so each segment's Python spectrum is rerun through
bench/spectrum_one.py (which mirrors segment_spectrum and emits
regions), verified field-for-field against the golden row, and then
matched to the C++ events by the unordered pair of region quad sets.
Rerunning Python also yields the Python-vs-C++ timing table for free.

    uv run python bench/geodesic_regress.py analytic
    uv run python bench/geodesic_regress.py run --quick 5
    uv run python bench/geodesic_regress.py run            # full 185 x 2
    uv run python bench/geodesic_regress.py report

`run` checkpoints one JSON line per (segment, diagonal) into the
workdir, so a killed run resumes with --resume at no cost.
"""
from __future__ import annotations

import argparse
import json
import math
import struct
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import spectrum_one  # noqa: E402

GOLDEN = {0: Path("out/spectrum_final_d0.json"),
          1: Path("out/spectrum_final_d1.json")}
ENGINE = Path("engines/geodesic")
WORKDIR = Path("out/geodesic_regress")
MEDIAN_TOL = 1e-9
COMPARE_FIELDS = ("separation_mm", "same_component", "ambiguous",
                  "self_touching", "n_pairs", "du_max", "endpoint_exact",
                  "median_intersection_length_vx")


# ------------------------------------------------------------- fixtures
def write_fixture_atlas(path: Path, P: np.ndarray, V: np.ndarray) -> None:
    """A one-surface WCAT atlas, byte-identical to windcheck.atlas layout."""
    P = np.ascontiguousarray(P, dtype="<f4")
    V = np.ascontiguousarray(V, dtype=np.uint8)
    nv, nu = V.shape
    assert P.shape == (nv, nu, 3)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        fh.write(b"WCAT")
        fh.write(struct.pack("<II", 1, 1))
        fh.write(struct.pack("<iII", -1, nv, nu))
        fh.write(P.tobytes())
        fh.write(V.tobytes())


def write_fixture_pairs(path: Path, quads) -> None:
    with path.open("w") as fh:
        fh.write("v1,u1,v2,u2,verdict,penetration,angle_deg\n")
        for v1, u1, v2, u2 in quads:
            fh.write(f"{v1},{u1},{v2},{u2},transverse,1,10\n")


def plane(n: int = 7, spacing: float = 1.0):
    v, u = np.mgrid[0:n, 0:n].astype(np.float64)
    P = np.stack([u * spacing, v * spacing, np.zeros_like(u)], axis=-1)
    return P, np.ones((n, n), dtype=bool)


# ----------------------------------------------------------- comparison
def _num_eq(a, b, tol=0.0):
    if a is None or b is None:
        return a is None and b is None
    fa, fb = float(a), float(b)
    if math.isnan(fa) or math.isnan(fb):
        return math.isnan(fa) and math.isnan(fb)
    return fa == fb if tol == 0.0 else abs(fa - fb) <= tol


def event_diffs(ref: dict, got: dict) -> list[str]:
    """Field-level mismatches between a reference event and a C++ event."""
    out = []
    for f in COMPARE_FIELDS:
        a, b = ref.get(f), got.get(f)
        if f == "median_intersection_length_vx":
            ok = _num_eq(a, b, MEDIAN_TOL)
        elif f == "separation_mm":
            ok = _num_eq(a, b)          # exact equality of the rounded value
        else:
            ok = a == b
        if not ok:
            out.append(f"{f}: py={a!r} cpp={b!r}")
    return out


def signature(ev: dict):
    a = frozenset(map(tuple, ev["region_a"]))
    b = frozenset(map(tuple, ev["region_b"]))
    return frozenset((a, b))


def sort_key(ev: dict):
    m = ev.get("median_intersection_length_vx")
    s = ev.get("separation_mm")
    return (ev.get("n_pairs"), ev.get("du_max"),
            -1e30 if m is None or (isinstance(m, float) and math.isnan(m))
            else float(m),
            -1e30 if s is None else float(s))


def match_and_compare(py_events: list[dict], cpp_events: list[dict]):
    """Match by region signature (multisets within a signature are paired
    after sorting by their scalar fields). Returns (n_exact, mismatches)."""
    mism: list[str] = []
    by_sig_py: dict = {}
    by_sig_cpp: dict = {}
    for e in py_events:
        by_sig_py.setdefault(signature(e), []).append(e)
    for e in cpp_events:
        by_sig_cpp.setdefault(signature(e), []).append(e)
    n_exact = 0
    for sig in set(by_sig_py) | set(by_sig_cpp):
        ps = sorted(by_sig_py.get(sig, []), key=sort_key)
        cs = sorted(by_sig_cpp.get(sig, []), key=sort_key)
        if len(ps) != len(cs):
            mism.append(f"signature {sorted(map(sorted, sig))}: "
                        f"{len(ps)} python events vs {len(cs)} c++ events")
            continue
        for p, c in zip(ps, cs):
            d = event_diffs(p, c)
            if d:
                mism.append(f"signature n_pairs={p['n_pairs']} "
                            f"du_max={p['du_max']}: " + "; ".join(d))
            else:
                n_exact += 1
    return n_exact, mism


def golden_check(golden_events: list[dict], py_events: list[dict]) -> list[str]:
    """The spectrum_one mirror must reproduce the golden rows in order."""
    if len(golden_events) != len(py_events):
        return [f"golden has {len(golden_events)} events, "
                f"rerun produced {len(py_events)}"]
    out = []
    for k, (a, b) in enumerate(zip(golden_events, py_events)):
        bad = []
        for f in COMPARE_FIELDS + ("distance_exact",):
            x, y = a.get(f), b.get(f)
            ok = _num_eq(x, y) if isinstance(x, float) or isinstance(y, float) \
                else x == y
            if not ok:
                bad.append(f"{f}: golden={x!r} rerun={y!r}")
        if bad:
            out.append(f"event[{k}]: " + "; ".join(bad))
    return out


# -------------------------------------------------------------- running
def run_engine(engine: Path, atlas: Path, csv: Path, out: Path,
               diagonal: int, voxel_um: float, threads: int = 0,
               maxedge: float = 60.0):
    t = time.time()
    r = subprocess.run(
        [str(engine), str(atlas), str(csv), str(out), str(threads),
         str(diagonal), str(maxedge), repr(voxel_um)],
        capture_output=True, text=True)
    wall = time.time() - t
    if r.returncode != 0:
        raise RuntimeError(f"geodesic failed ({r.returncode}): {r.stderr}")
    stats = json.loads(r.stdout) if r.stdout.strip() else {}
    events = json.loads(out.read_text())
    return events, stats, wall


def build_atlas(mesh: Path, out: Path) -> None:
    from windcheck import atlas as watlas

    class _E:
        def __init__(self, p):
            self.path, self.winding = p, None
    watlas.write_atlas([_E(mesh)], out)


def segment_list(diagonals):
    """(corpus, segment, voxel_um, pairs, events) in golden order + goldens."""
    golden = {}
    for d in diagonals:
        golden[d] = {(r["corpus"], r["segment"]): r
                     for r in json.loads(GOLDEN[d].read_text())}
    first = json.loads(GOLDEN[diagonals[0]].read_text())
    segs = [(r["corpus"], r["segment"], r["voxel_um"], r["pairs"],
             r["events_measured"]) for r in first]
    return segs, golden


def checkpoint_rows(workdir: Path) -> list[dict]:
    """All rows from every regress*.jsonl in the workdir (parallel runs
    write tagged sidecars; they all count)."""
    rows = []
    for f in sorted(workdir.glob("regress*.jsonl")):
        rows += [json.loads(x) for x in f.read_text().splitlines() if x]
    return rows


def cmd_run(a) -> int:
    workdir = a.workdir
    workdir.mkdir(parents=True, exist_ok=True)
    tag = f"-{a.tag}" if a.tag else ""
    jsonl = workdir / f"regress{tag}.jsonl"
    done = set()
    if a.resume:
        done = {(r["corpus"], r["segment"], r["diagonal"])
                for r in checkpoint_rows(workdir)}
    elif jsonl.exists():
        jsonl.unlink()

    segs, golden = segment_list(a.diagonal)
    if a.corpus:
        segs = [s for s in segs if s[0] == a.corpus]
    if a.segment:
        segs = [s for s in segs if a.segment in s[1]]
    if a.quick:
        # smallest segment per corpus that actually HAS events to compare
        by_corpus: dict = {}
        for s in sorted(segs, key=lambda s: s[3]):
            if s[4] > 0:
                by_corpus.setdefault(s[0], s)
        segs = sorted(by_corpus.values(), key=lambda s: s[3])[:a.quick]

    total_mism = 0
    for corpus, segment, vx, _pairs, _nev in segs:
        todo = [d for d in a.diagonal if (corpus, segment, d) not in done]
        if not todo:
            continue
        mesh, _, _ = spectrum_one.find_segment(corpus, segment, a.diagonal[0])
        atlas_bin = workdir / f"atlas_tmp{tag}.bin"
        build_atlas(mesh, atlas_bin)
        for d in todo:
            _, csv, _ = spectrum_one.find_segment(corpus, segment, d)
            out_json = workdir / f"cpp_events{tag}.json"
            cpp_events, stats, wall = run_engine(
                a.engine, atlas_bin, csv, out_json, d, vx, a.threads)
            py = spectrum_one.run_one(corpus, segment, d, emit_regions=True)
            gold = golden[d][(corpus, segment)]["events"]
            g_mism = golden_check(gold, py["events"])
            n_exact, mism = match_and_compare(py["events"], cpp_events)
            row = {"corpus": corpus, "segment": segment, "diagonal": d,
                   "voxel_um": vx, "pairs": py["pairs"],
                   "events_golden": len(gold), "events_cpp": len(cpp_events),
                   "exact": n_exact, "mismatched": len(mism),
                   "golden_rerun_mismatches": g_mism,
                   "mismatches": mism,
                   "t_python_s": py["elapsed_s"],
                   "t_cpp_s": stats.get("seconds"), "t_cpp_wall_s": round(wall, 3)}
            with jsonl.open("a") as fh:
                fh.write(json.dumps(row) + "\n")
            total_mism += len(mism) + len(g_mism)
            flag = "OK " if not (mism or g_mism) else "MISMATCH"
            print(f"{flag} {corpus:10s} {segment[:36]:38s} d{d} "
                  f"ev {len(cpp_events):5d}  exact {n_exact:5d}  "
                  f"py {py['elapsed_s']:8.1f}s  cpp {wall:7.1f}s", flush=True)
            for m in g_mism:
                print(f"    GOLDEN-RERUN {m}", flush=True)
            for m in mism:
                print(f"    {m}", flush=True)
        atlas_bin.unlink(missing_ok=True)
    print(f"\nrun complete: {'ZERO mismatches' if total_mism == 0 else f'{total_mism} MISMATCHES'}")
    return 0 if total_mism == 0 else 1


def cmd_report(a) -> int:
    rows = checkpoint_rows(a.workdir)
    agg: dict = {}
    monster = None
    for r in rows:
        k = (r["corpus"], r["diagonal"])
        g = agg.setdefault(k, {"segs": 0, "events": 0, "exact": 0,
                               "mism": 0, "t_py": 0.0, "t_cpp": 0.0})
        g["segs"] += 1
        g["events"] += r["events_golden"]
        g["exact"] += r["exact"]
        g["mism"] += r["mismatched"] + len(r["golden_rerun_mismatches"])
        g["t_py"] += r["t_python_s"] or 0.0
        g["t_cpp"] += r["t_cpp_wall_s"] or 0.0
        if monster is None or r["pairs"] > monster["pairs"]:
            monster = r
    print(f"{'corpus':12s} d  segs  events   exact  mism   t_python    t_cpp   speedup")
    tot = {"segs": 0, "events": 0, "exact": 0, "mism": 0, "t_py": 0.0, "t_cpp": 0.0}
    for (corpus, d) in sorted(agg):
        g = agg[(corpus, d)]
        for k in tot:
            tot[k] += g[k]
        sp = g["t_py"] / g["t_cpp"] if g["t_cpp"] else float("nan")
        print(f"{corpus:12s} {d}  {g['segs']:4d}  {g['events']:6d}  "
              f"{g['exact']:6d}  {g['mism']:4d}  {g['t_py']:8.1f}s "
              f"{g['t_cpp']:7.1f}s  {sp:7.1f}x")
    sp = tot["t_py"] / tot["t_cpp"] if tot["t_cpp"] else float("nan")
    print(f"{'TOTAL':12s} -  {tot['segs']:4d}  {tot['events']:6d}  "
          f"{tot['exact']:6d}  {tot['mism']:4d}  {tot['t_py']:8.1f}s "
          f"{tot['t_cpp']:7.1f}s  {sp:7.1f}x")
    if monster:
        print(f"monster: {monster['corpus']} {monster['segment'][:40]} "
              f"d{monster['diagonal']} pairs {monster['pairs']} "
              f"events {monster['events_golden']} "
              f"py {monster['t_python_s']:.1f}s cpp {monster['t_cpp_wall_s']:.1f}s")
    bad = [r for r in rows
           if r["mismatched"] or r["golden_rerun_mismatches"]]
    for r in bad:
        print(f"\nMISMATCHES {r['corpus']} {r['segment']} d{r['diagonal']}:")
        for m in r["golden_rerun_mismatches"]:
            print(f"  GOLDEN-RERUN {m}")
        for m in r["mismatches"]:
            print(f"  {m}")
    return 1 if bad else 0


# ------------------------------------------------------------- analytic
def cmd_analytic(a) -> int:
    """The 14 analytic cases of tests/test_intrinsic.py, driven through the
    engine on tiny synthetic atlases (plane / slit / strip constructions)."""
    fx = a.workdir / "fixtures"
    fx.mkdir(parents=True, exist_ok=True)
    fails = []

    def check(name, ok, detail=""):
        print(f"  {name:58s} {'ok' if ok else 'FAIL ' + detail}")
        if not ok:
            fails.append(name)

    def tq(atlas, diag, maxedge, *args):
        r = subprocess.run([str(a.engine), "--test", str(atlas), str(diag),
                            str(maxedge)] + [str(x) for x in args],
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"--test failed: {r.stderr}")
        return r.stdout.strip()

    def full(atlas, quads, diag=0, voxel_um=1000.0):
        csv = fx / "pairs.csv"
        out = fx / "events.json"
        write_fixture_pairs(csv, quads)
        events, _, _ = run_engine(a.engine, atlas, csv, out, diag, voxel_um)
        return events

    SQ2 = math.sqrt(2.0)
    apx = lambda x, y, tol=1e-5: abs(float(x) - y) <= tol  # noqa: E731

    # 1-2: exact row/diagonal distances (f32 edge weights: 1e-5 tolerance)
    P, V = plane(7)
    a7 = fx / "plane7.bin"
    write_fixture_atlas(a7, P, V)
    check("row distance (0,0)-(0,6) = 6", apx(tq(a7, 0, 60, "dist", 0, 0, 0, 6), 6.0))
    check("col distance (3,1)-(5,1) = 2", apx(tq(a7, 0, 60, "dist", 3, 1, 5, 1), 2.0))
    check("d0 (0,0)-(5,5) = 5*sqrt2", apx(tq(a7, 0, 60, "dist", 0, 0, 5, 5), 5 * SQ2))
    check("d1 (0,0)-(5,5) = 10", apx(tq(a7, 1, 60, "dist", 0, 0, 5, 5), 10.0))
    check("d1 (5,0)-(0,5) = 5*sqrt2", apx(tq(a7, 1, 60, "dist", 5, 0, 0, 5), 5 * SQ2))
    check("d0 (5,0)-(0,5) = 10", apx(tq(a7, 0, 60, "dist", 5, 0, 0, 5), 10.0))

    # 3: one-row slit disconnects (both vertices exist, no path)
    P, V = plane(7)
    V[0:6, 3] = False
    slit1 = fx / "slit1.bin"
    write_fixture_atlas(slit1, P, V)
    for d in (0, 1):
        check(f"one-row slit d{d}: vertices exist",
              int(tq(slit1, d, 60, "idx", 0, 0)) >= 0
              and int(tq(slit1, d, 60, "idx", 0, 6)) >= 0)
        check(f"one-row slit d{d}: disconnected",
              tq(slit1, d, 60, "dist", 0, 0, 0, 6) == "inf")

    # 4: two-row passage forces the exact 12+2*sqrt2 detour
    P, V = plane(7)
    V[0:5, 3] = False
    slit2 = fx / "slit2.bin"
    write_fixture_atlas(slit2, P, V)
    for d in (0, 1):
        check(f"two-row passage d{d}: 12+2*sqrt2",
              apx(tq(slit2, d, 60, "dist", 0, 0, 0, 6), 12 + 2 * SQ2))

    # 5: a single-vertex chain is not surface
    P, V = plane(5)
    V[:, 2] = False
    V[2, 2] = True
    chain = fx / "chain.bin"
    write_fixture_atlas(chain, P, V)
    check("lone bridge vertex not in graph", tq(chain, 0, 60, "idx", 2, 2) == "-1")
    check("chain leaves two components", tq(chain, 0, 60, "ncomp") == "2")

    # 6: maxedge drops the quad exactly as selfcross drops it
    P, V = plane(4)
    P[1, 1, 2] = 5.0
    spike = fx / "spike.bin"
    write_fixture_atlas(spike, P, V)
    check("maxedge: spiked vertex belongs to no quad",
          tq(spike, 0, 2.0, "idx", 1, 1) == "-1")
    check("maxedge: corner quad dropped too", tq(spike, 0, 2.0, "idx", 0, 0) == "-1")
    check("maxedge: rest stays one sheet", tq(spike, 0, 2.0, "ncomp") == "1")

    # 7: oriented grouping, straight and swapped
    P, V = plane(22)
    p22 = fx / "plane22.bin"
    write_fixture_atlas(p22, P, V)
    evs = full(p22, [(0, 0, 10, 0), (1, 0, 11, 0), (12, 0, 2, 0), (0, 20, 10, 20)])
    big = evs[0] if evs else {}
    regs = ({tuple(map(tuple, sorted(big.get("region_a", [])))),
             tuple(map(tuple, sorted(big.get("region_b", []))))} if evs else set())
    check("grouping: two events", len(evs) == 2)
    check("grouping: big event has 3 pairs, unambiguous",
          bool(evs) and big["n_pairs"] == 3 and not big["ambiguous"]
          and not big["self_touching"])
    check("grouping: parity-normalised regions",
          regs == {((0, 0), (1, 0), (2, 0)), ((10, 0), (11, 0), (12, 0))})

    # 8-9: orientation conflict -> ambiguous, reported not measured
    P, V = plane(8)
    p8 = fx / "plane8.bin"
    write_fixture_atlas(p8, P, V)
    evs = full(p8, [(0, 0, 2, 0), (1, 0, 1, 1)])
    check("conflict: one ambiguous self-touching event",
          len(evs) == 1 and evs[0]["ambiguous"] and evs[0]["self_touching"])
    check("conflict: never measured",
          bool(evs) and evs[0]["separation_mm"] is None
          and evs[0]["same_component"] is None)

    # 10: tri-tri interval method (hard-coded values live in --selftest)
    st = subprocess.run([str(a.engine), "--selftest"], capture_output=True)
    check("selftest: tri-tri segment + nanmedian", st.returncode == 0)

    # 11: planted transverse crossing, measured barycentrically
    yz = [(0, -0.5), (1, -0.5), (2, -0.5), (3, -0.5),
          (3.5, 2.0), (1.5, 2.0), (1.5, -2.0)]
    P = np.array([[[v, y, z] for (y, z) in yz] for v in (0.0, 1.0)])
    V = np.ones((2, 7), dtype=bool)
    strip = fx / "strip.bin"
    write_fixture_atlas(strip, P, V)
    evs = full(strip, [(0, 1, 0, 5)])
    check("planted crossing: one exact-endpoint event",
          len(evs) == 1 and evs[0]["endpoint_exact"] is True
          and evs[0]["same_component"] is True)
    check("planted crossing: separation 6+sqrt(6.5)",
          bool(evs) and evs[0]["separation_mm"] is not None
          and apx(evs[0]["separation_mm"], 6 + math.sqrt(6.5), 1e-3))

    # 12: seeded distance with unequal offsets never mixes them
    want = 12 + 2 * SQ2
    check("seeded: plain detour", apx(tq(slit2, 0, 60, "seeded",
          1, 0, 0, 0.0, 1, 0, 6, 0.0), want))
    check("seeded: decoy with offset 100 does not shortcut",
          apx(tq(slit2, 0, 60, "seeded",
                 2, 0, 0, 0.0, 0, 5, 100.0, 1, 0, 6, 0.0), min(want, 101.0)))

    # 13: disconnection is an event property; corner fallback distance exact
    P, V = plane(10)
    V[:, 4:6] = False
    isl = fx / "islands.bin"
    write_fixture_atlas(isl, P, V)
    check("islands: two components", tq(isl, 0, 60, "ncomp") == "2")
    evs = full(isl, [(0, 0, 0, 7)])
    check("islands: cross-island event has null separation",
          len(evs) == 1 and evs[0]["separation_mm"] is None
          and evs[0]["same_component"] is False)
    evs = full(isl, [(0, 0, 0, 2)])
    check("islands: same-island corner fallback = 1.0",
          len(evs) == 1 and evs[0]["same_component"] is True
          and evs[0]["endpoint_exact"] is False
          and _num_eq(evs[0]["separation_mm"], 1.0))

    # 14: events emitted largest-first (the Python cap is caller-side)
    P, V = plane(14)
    p14 = fx / "plane14.bin"
    write_fixture_atlas(p14, P, V)
    evs = full(p14, [(0, 0, 0, 8), (5, 0, 5, 8), (6, 1, 6, 9)])
    check("largest event first", len(evs) == 2 and evs[0]["n_pairs"] == 2
          and evs[1]["n_pairs"] == 1)

    print(f"\nanalytic: {'ALL PASS' if not fails else f'{len(fails)} FAILURES'}")
    return 0 if not fails else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("analytic", "run", "report"))
    ap.add_argument("--engine", type=Path, default=ENGINE)
    ap.add_argument("--workdir", type=Path, default=WORKDIR)
    ap.add_argument("--diagonal", type=int, nargs="+", default=[0, 1],
                    choices=(0, 1))
    ap.add_argument("--corpus", default="")
    ap.add_argument("--segment", default="",
                    help="substring filter on segment names")
    ap.add_argument("--quick", type=int, default=0,
                    help="N smallest segments, at most one per corpus")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--tag", default="",
                    help="suffix for this run's checkpoint/temp files so "
                         "several corpora can run in parallel")
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    if a.mode == "analytic":
        return cmd_analytic(a)
    if a.mode == "run":
        return cmd_run(a)
    return cmd_report(a)


if __name__ == "__main__":
    sys.exit(main())
