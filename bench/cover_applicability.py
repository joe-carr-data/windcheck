"""Alternate-cover applicability gates (PREREG-COVER-APPLICABILITY.md,
Addenda 1-2). Engine-authoritative; Python never decides a verdict.

G1: fraction of events touching >= 1 ACTIONABLE contested pixel
    (>= 2 DISTINCT candidate XYZ tuples).
G2: fraction SOLVED by the corrected greedy. Outcomes are three-state
    (solved / unsolved_exhausted / censored_timeout_or_cap) with
    lower/upper bounds per Addendum 2.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import struct
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import tifffile

REPO_ROOT = Path(__file__).resolve().parents[1]
for _p in (REPO_ROOT / "src", REPO_ROOT / "bench"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from crossing_analyse import components, events as region_events  # noqa: E402
from excise_segment import parse_census_csv                       # noqa: E402
from windcheck.pipeline import ENGINE                             # noqa: E402

PROCEDURE_VERSION = "addendum3"
CAND_DTYPE = np.dtype([("pi", "<u8"), ("face", "<i4"),
                       ("x", "<f4"), ("y", "<f4"), ("z", "<f4")])
MAXEDGE = 60.0
TIME_BUDGET_S = 900.0
SUBPROC_TIMEOUT_S = 120.0


class EngineError(RuntimeError):
    pass


def load_rows_events(csv_dir: Path, diag: int, H: int, W: int):
    rec = parse_census_csv(csv_dir / f"s12345_d{diag}.csv", diag, H, W)
    rows = [r for r in rec["rows"] if r["verdict"] == "transverse"]
    arr = np.array([(r["q1"][0], r["q1"][1], r["q2"][0], r["q2"][1])
                    for r in rows],
                   dtype=[("v1", "i4"), ("u1", "i4"),
                          ("v2", "i4"), ("u2", "i4")])
    label, _ = components(arr)
    by_event = defaultdict(list)
    for r in rows:
        a = label[(r["q1"][0], r["q1"][1])]
        b = label[(r["q2"][0], r["q2"][1])]
        by_event[(a, b) if a <= b else (b, a)].append(r)
    return by_event


def corner_pixels(cell):
    v, u = cell
    return [(v, u), (v + 1, u), (v, u + 1), (v + 1, u + 1)]


def incident_quads(pixel, H, W):
    pv, pu = pixel
    out = []
    for cv in (pv - 1, pv):
        for cu in (pu - 1, pu):
            if 0 <= cv < H - 1 and 0 <= cu < W - 1:
                out.append((cv, cu))
    return out


def patch_of(pixel, H, W):
    cells = set()
    for q in incident_quads(pixel, H, W):
        for dv in (-1, 0, 1):
            for du in (-1, 0, 1):
                v, u = q[0] + dv, q[1] + du
                if 0 <= v < H - 1 and 0 <= u < W - 1:
                    cells.add((v, u))
    return cells


def row_key(r, diag):
    c1, c2 = tuple(r["q1"]), tuple(r["q2"])
    t1, t2 = r["t1"], r["t2"]
    if (c1, t1) <= (c2, t2):
        return (diag, c1[0], c1[1], t1, c2[0], c2[1], t2, "transverse")
    return (diag, c2[0], c2[1], t2, c1[0], c1[1], t1, "transverse")


class PatchTester:
    """Both-diagonal engine census of a cell set with overrides.

    Returns (contact multiset over schema-v2 triangle keys of ALL verdict
    classes for in-set pairs, retained-quad set under the engine's own
    six-edge maxedge criterion). Raises EngineError on failure."""

    def __init__(self, bands, workdir: Path):
        self.x, self.y, self.z = bands
        self.H, self.W = self.x.shape
        self.work = workdir
        self.n_calls = 0

    def _crop(self, cells, overrides):
        vs = [c[0] for c in cells]
        us = [c[1] for c in cells]
        v0, v1 = min(vs), max(vs) + 2
        u0, u1 = min(us), max(us) + 2
        cx = self.x[v0:v1, u0:u1].copy()
        cy = self.y[v0:v1, u0:u1].copy()
        cz = self.z[v0:v1, u0:u1].copy()
        for (pv, pu), (nx, ny, nz) in overrides.items():
            if v0 <= pv < v1 and u0 <= pu < u1:
                cx[pv - v0, pu - u0] = nx
                cy[pv - v0, pu - u0] = ny
                cz[pv - v0, pu - u0] = nz
        return v0, u0, cx, cy, cz

    def retained_quads(self, cells, overrides):
        """Engine criterion mirrored exactly: 4 valid corners AND max of
        the six pairwise corner distances <= MAXEDGE."""
        v0, u0, cx, cy, cz = self._crop(cells, overrides)
        P = np.stack([cx, cy, cz], axis=-1).astype(np.float64)
        V = ~((cx == -1) & (cy == -1) & (cz == -1))
        out = set()
        for (cv, cu) in cells:
            i, j = cv - v0, cu - u0
            if not (V[i, j] and V[i + 1, j] and V[i, j + 1]
                    and V[i + 1, j + 1]):
                continue
            c00, c10 = P[i, j], P[i + 1, j]
            c01, c11 = P[i, j + 1], P[i + 1, j + 1]
            e = max(np.linalg.norm(c01 - c00), np.linalg.norm(c11 - c01),
                    np.linalg.norm(c10 - c11), np.linalg.norm(c00 - c10),
                    np.linalg.norm(c11 - c00), np.linalg.norm(c10 - c01))
            if e <= MAXEDGE:
                out.add((cv, cu))
        return out

    def census(self, cells: set, overrides: dict) -> Counter:
        v0, u0, cx, cy, cz = self._crop(cells, overrides)
        h, w = cx.shape
        with tempfile.NamedTemporaryFile(dir=self.work, suffix=".bin",
                                         delete=False) as fh:
            abin = Path(fh.name)
            fh.write(b"WCAT" + struct.pack("<II", 1, 1))
            fh.write(struct.pack("<iII", -1, h, w))
            fh.write(np.stack([cx, cy, cz], axis=-1)
                     .astype("<f4").tobytes())
            fh.write((~((cx == -1) & (cy == -1) & (cz == -1)))
                     .astype(np.uint8).tobytes())
        found: Counter = Counter()
        try:
            for diag in (0, 1):
                out = abin.with_suffix(f".d{diag}.csv")
                try:
                    r = subprocess.run(
                        [str(ENGINE), str(abin), str(out), "1", "40.0",
                         "1", str(diag), str(MAXEDGE)],
                        capture_output=True, text=True,
                        timeout=SUBPROC_TIMEOUT_S)
                except subprocess.TimeoutExpired:
                    raise EngineError("engine subprocess timeout")
                self.n_calls += 1
                if r.returncode != 0:
                    raise EngineError(f"engine rc={r.returncode}")
                with open(out) as fh2:
                    fh2.readline()
                    for line in fh2:
                        p = line.split(",")
                        c1 = (int(p[0]) + v0, int(p[1]) + u0)
                        c2 = (int(p[2]) + v0, int(p[3]) + u0)
                        if c1 not in cells or c2 not in cells:
                            continue
                        t1, t2 = int(p[7]), int(p[8])
                        verdict = p[4]
                        a, b = ((c1, t1), (c2, t2))
                        if a <= b:
                            key = (diag, c1[0], c1[1], t1,
                                   c2[0], c2[1], t2, verdict)
                        else:
                            key = (diag, c2[0], c2[1], t2,
                                   c1[0], c1[1], t1, verdict)
                        found[key] += 1
                out.unlink(missing_ok=True)
        finally:
            abin.unlink(missing_ok=True)
        return found


ALLOW_NEW_GRAZING = False


def accept_ok(before: Counter, after: Counter, target_in: set) -> bool:
    """Addendum-2 acceptance on one union: target rows strictly decrease;
    no key of any class new or increased."""
    tb = sum(c for k, c in before.items() if k in target_in)
    ta = sum(c for k, c in after.items() if k in target_in)
    if not ta < tb:
        return False
    # Addendum 3: no-increase applies to EVERY key, targets included --
    # intermediate resurrection is not forgiven.
    for k, c in after.items():
        if c > before.get(k, 0):
            if ALLOW_NEW_GRAZING and k[7] == "grazing":
                continue
            return False
    return True


def greedy_solve(rows, tester, idx_of, cand, max_switches, rec, diag):
    H, W = tester.H, tester.W
    t0 = time.time()
    targets = {row_key(r, diag): r for r in rows}
    unresolved = set(targets)
    overrides: dict = {}
    switched: set = set()

    def local_union(extra_pixel=None):
        cells = set()
        for p in switched:
            cells |= patch_of(p, H, W)
        if extra_pixel is not None:
            cells |= patch_of(extra_pixel, H, W)
        return cells

    for _ in range(max_switches):
        if time.time() - t0 > TIME_BUDGET_S:
            rec["state"] = "censored_timeout_or_cap"
            return False
        improved = False
        for r in rows:
            k = row_key(r, diag)
            if k not in unresolved:
                continue
            if time.time() - t0 > TIME_BUDGET_S:
                rec["state"] = "censored_timeout_or_cap"
                return False
            for cell in (tuple(r["q1"]), tuple(r["q2"])):
                for p in corner_pixels(cell):
                    # Addendum-2 skip: incidence with unresolved targets
                    inc = set(incident_quads(p, H, W))
                    if not any((tuple(targets[u]["q1"]) in inc
                                or tuple(targets[u]["q2"]) in inc)
                               for u in unresolved):
                        continue
                    pi = p[0] * W + p[1]
                    if pi not in idx_of:
                        continue
                    s, c = idx_of[pi]
                    if c < 2:
                        continue
                    union = local_union(p)
                    before = tester.census(union, overrides)
                    rq_before = tester.retained_quads(union, overrides)
                    tin = {t for t in targets
                           if (tuple(targets[t]["q1"]) in union
                               and tuple(targets[t]["q2"]) in union)}
                    for kk in range(s, s + c):
                        alt = (float(cand["x"][kk]), float(cand["y"][kk]),
                               float(cand["z"][kk]))
                        old = overrides.get(p)
                        base = (float(tester.x[p]), float(tester.y[p]),
                                float(tester.z[p]))
                        if (old or base) == alt:
                            continue
                        trial = dict(overrides)
                        trial[p] = alt
                        after = tester.census(union, trial)
                        if not accept_ok(before, after, tin):
                            continue
                        rq_after = tester.retained_quads(union, trial)
                        if not rq_before <= rq_after:
                            continue
                        overrides.clear()
                        overrides.update(trial)
                        switched.add(p)
                        rec["switches"] += 1
                        # update unresolved from this union's after-state
                        for t in list(unresolved):
                            if t in tin and after.get(t, 0) == 0:
                                unresolved.discard(t)
                        improved = True
                        break
                    if improved:
                        break
                if improved:
                    break
            if improved:
                break
        if not improved:
            rec["state"] = "unsolved_exhausted"
            return False
        if not unresolved:
            # FINAL global verification (Addendum 2.2)
            big = set()
            for r in rows:
                for cc in (tuple(r["q1"]), tuple(r["q2"])):
                    for p in corner_pixels(cc):
                        big |= patch_of(p, H, W)
            for p in switched:
                big |= patch_of(p, H, W)
            baseline = tester.census(big, {})
            final = tester.census(big, overrides)
            rqb = tester.retained_quads(big, {})
            rqa = tester.retained_quads(big, overrides)
            tkeys = set(targets)
            if (all(final.get(t, 0) == 0 for t in tkeys)
                    and all(c <= baseline.get(k, 0)
                            for k, c in final.items() if k not in tkeys)
                    and rqb <= rqa):
                rec["state"] = "solved"
                return True
            rec["state"] = "unsolved_exhausted"
            return False
    rec["state"] = "censored_timeout_or_cap"
    return False


_G = None


def _worker(job):
    diag, key, rows = job
    bands, idx_of, cand, actionable_set, max_switches = _G
    with tempfile.TemporaryDirectory() as td:
        tester = PatchTester(bands, Path(td))
        W = tester.W
        pixset = set()
        for r in rows:
            for c in (r["q1"], r["q2"]):
                pixset |= set(corner_pixels(tuple(c)))
        touch = any((pv * W + pu) in actionable_set for pv, pu in pixset)
        rec = {"diag": diag, "event": key, "n_rows": len(rows),
               "contested_touch": touch, "state": "unsolved_exhausted",
               "switches": 0, "engine_calls": 0}
        if touch:
            try:
                greedy_solve(rows, tester, idx_of, cand, max_switches,
                             rec, diag)
            except EngineError as e:
                rec["state"] = "engine_error"
                rec["error"] = str(e)
        rec["engine_calls"] = tester.n_calls
        return rec


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None) -> int:
    global _G
    ap = argparse.ArgumentParser()
    ap.add_argument("--export", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--csv-dir", default="out/fiesta/unroll")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--max-switches", type=int, default=64)
    args = ap.parse_args(argv)

    export = Path(args.export)
    bands = tuple(np.asarray(tifffile.imread(export / f"{a}.tif"))
                  for a in ("x", "y", "z"))
    H, W = bands[0].shape

    csize = Path(args.candidates).stat().st_size
    if csize % CAND_DTYPE.itemsize != 0:
        raise SystemExit("candidate stream size not record-divisible")
    cand = np.fromfile(args.candidates, dtype=CAND_DTYPE)
    if not (np.isfinite(cand["x"]).all() and np.isfinite(cand["y"]).all()
            and np.isfinite(cand["z"]).all()):
        raise SystemExit("non-finite candidate coordinates")
    if cand["pi"].max() >= H * W or cand["face"].min() < 0:
        raise SystemExit("candidate bounds violation")
    order = np.argsort(cand["pi"], kind="stable")
    cand = cand[order]
    pis, starts, counts = np.unique(cand["pi"], return_index=True,
                                    return_counts=True)
    # first streamed XYZ equals emitted first-cover at covered px
    fx = bands[0].reshape(-1)[pis[:1000]]
    if not np.allclose(cand["x"][starts[:1000]], fx, atol=0, rtol=0):
        raise SystemExit("first-cover mismatch vs export")
    # actionable: >= 2 DISTINCT xyz tuples per pixel
    xyz = np.stack([cand["x"], cand["y"], cand["z"]], axis=1)
    first = np.repeat(xyz[starts], counts, axis=0)
    differs = (xyz != first).any(axis=1)
    diff_per_pixel = np.add.reduceat(differs.astype(np.int64), starts)
    actionable = pis[(counts >= 2) & (diff_per_pixel > 0)]
    actionable_set = set(int(p) for p in actionable)
    idx_of = {int(p): (int(s), int(c))
              for p, s, c in zip(pis, starts, counts)}
    print(f"candidates: {len(cand):,} records, {len(pis):,} px, "
          f"{int((counts >= 2).sum()):,} raw contested, "
          f"{len(actionable):,} actionable", flush=True)

    jobs = []
    for diag in (0, 1):
        by_event = load_rows_events(Path(args.csv_dir), diag, H, W)
        for key, rows in sorted(by_event.items()):
            jobs.append((diag, list(key),
                         [{"q1": tuple(r["q1"]), "q2": tuple(r["q2"]),
                           "t1": r["t1"], "t2": r["t2"]} for r in rows]))
    binding = {"procedure": PROCEDURE_VERSION,
               "candidates_sha256": sha256_file(args.candidates),
               "x_tif_sha256": sha256_file(export / "x.tif"),
               "y_tif_sha256": sha256_file(export / "y.tif"),
               "z_tif_sha256": sha256_file(export / "z.tif"),
               "csv_d0_sha256": sha256_file(Path(args.csv_dir) / "s12345_d0.csv"),
               "csv_d1_sha256": sha256_file(Path(args.csv_dir) / "s12345_d1.csv")}
    ckpt = Path(args.out).with_suffix(".ckpt.jsonl")
    done = {}
    if ckpt.exists():
        lines = ckpt.read_text().splitlines()
        if lines and json.loads(lines[0]) != binding:
            raise SystemExit("checkpoint binding mismatch; remove it")
        for line in lines[1:]:
            if line:
                r = json.loads(line)
                done[(r["diag"], tuple(r["event"]))] = r
    else:
        ckpt.write_text(json.dumps(binding) + "\n")
    pending = [j for j in jobs if (j[0], tuple(j[1])) not in done]
    print(f"{len(jobs)} events; checkpoint done={len(done)} "
          f"pending={len(pending)}; workers={args.workers}", flush=True)

    import multiprocessing as mp
    ctx = mp.get_context("fork")
    _G = (bands, idx_of, cand, actionable_set, args.max_switches)
    t0 = time.time()
    results = list(done.values())
    with ctx.Pool(args.workers) as pool, open(ckpt, "a") as ck:
        for i, res in enumerate(pool.imap_unordered(_worker, pending,
                                                    chunksize=1)):
            if res["state"] == "engine_error":
                raise SystemExit(f"ENGINE ERROR, aborting: {res}")
            results.append(res)
            ck.write(json.dumps(res) + "\n")
            ck.flush()
            if (i + 1) % 25 == 0 or i + 1 == len(pending):
                print(f"  {i+1}/{len(pending)} ({time.time()-t0:.0f}s)",
                      flush=True)

    doc = {"prereg": "notes/PREREG-COVER-APPLICABILITY.md (addenda 1-2)",
           "binding": binding,
           "candidate_stats": {
               "records": int(len(cand)), "pixels": int(len(pis)),
               "raw_contested": int((counts >= 2).sum()),
               "actionable_contested": int(len(actionable))},
           "diagonals": {}}
    g1_pass, g2_lower_pass, g2_upper_fail, indeterminate = [], [], [], []
    for diag in (0, 1):
        evs = [r for r in results if r["diag"] == diag]
        n = len(evs)
        g1 = sum(1 for r in evs if r["contested_touch"])
        solved = sum(1 for r in evs if r["state"] == "solved")
        censored = sum(1 for r in evs
                       if r["state"] == "censored_timeout_or_cap")
        lower = solved / n
        upper = (solved + censored) / n
        f1 = g1 / n
        doc["diagonals"][f"d{diag}"] = {
            "n_events": n, "g1_touch_actionable": g1, "g1_fraction": f1,
            "solved": solved, "censored": censored,
            "unsolved_exhausted": n - solved - censored
            - sum(1 for r in evs if not r["contested_touch"]
                  and r["state"] != "solved") if False else
            sum(1 for r in evs if r["state"] == "unsolved_exhausted"),
            "g2_lower": lower, "g2_upper": upper,
            "engine_calls": sum(r["engine_calls"] for r in evs),
            "events": evs}
        g1_pass.append(f1 >= 0.50)
        g2_lower_pass.append(lower >= 0.30)
        g2_upper_fail.append(upper < 0.30)
        indeterminate.append(lower < 0.30 <= upper)
        print(f"d{diag}: n={n} G1={f1:.3f} (>=0.50 {f1 >= 0.50}) | "
              f"solved={solved} censored={censored} "
              f"G2 lower={lower:.3f} upper={upper:.3f}", flush=True)
    if all(g1_pass) and all(g2_lower_pass):
        verdict = "QUALIFIES"
    elif (not all(g1_pass)) or any(g2_upper_fail):
        verdict = "FAILS"
    else:
        verdict = "INDETERMINATE"
    doc["verdict"] = verdict
    Path(args.out).write_text(json.dumps(doc, indent=1))
    print("VERDICT:", verdict, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
