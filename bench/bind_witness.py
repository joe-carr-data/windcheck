"""Bind every witness to a FRESH census, by identity and not by count.

The packaged witness says where each input self-intersects. Until now it
was tied to its trace only by per-diagonal contact COUNTS matching the
certificate. Two entirely different contact sets with the same counts
would both pass -- and "where the defects are" is the benchmark's
distinguishing asset, so that is the wrong place to be approximate.

This recensuses each input from scratch with the frozen parameters and
requires the witness to equal the fresh census as a MULTISET of canonical
contact identities:

    (diagonal, sorted[(v1, u1, tri1), (v2, u2, tri2)])

sorted because a contact is an unordered pair, a multiset because the
same pair of quads can legitimately contact more than once across
triangle combinations, and per-diagonal because d0 and d1 are different
triangulations of the same lattice.

A count comparison is kept alongside, and reported separately, so the
record shows which of the two actually bit.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from crossing_census import census_one  # noqa: E402

EXCLUDE, CELL, MAXEDGE, THREADS = 1, 40.0, 60.0, 3


def multiset_from_csv(path: Path, diagonal: int) -> collections.Counter:
    c: collections.Counter = collections.Counter()
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["verdict"] != "transverse":
                continue
            a = (int(row["v1"]), int(row["u1"]), int(row["tri1"]))
            b = (int(row["v2"]), int(row["u2"]), int(row["tri2"]))
            c[(diagonal, *sorted((a, b)))] += 1
    return c


def multiset_from_witness(npz: Path) -> collections.Counter:
    c: collections.Counter = collections.Counter()
    if not npz.is_file():
        return c
    z = np.load(npz)
    for d, v1, u1, t1, v2, u2, t2 in zip(
            z["diagonal"], z["v1"], z["u1"], z["tri1"],
            z["v2"], z["u2"], z["tri2"]):
        a = (int(v1), int(u1), int(t1))
        b = (int(v2), int(u2), int(t2))
        c[(int(d), *sorted((a, b)))] += 1
    return c


def negative_control(fresh: collections.Counter,
                     shipped: collections.Counter) -> dict:
    """A mutation that PRESERVES every count and changes one identity.

    This is the case the previous binding could not see: move one contact
    to a different pair of quads and the per-diagonal totals are
    unchanged, so a count comparison still passes. If the identity
    comparison does not fail here, it is not doing anything.
    """
    if not shipped:
        return {"ran": False, "why": "this trace has no contacts to mutate"}
    victim = next(iter(shipped))
    mutated = collections.Counter(shipped)
    mutated[victim] -= 1
    if mutated[victim] == 0:
        del mutated[victim]
    d, qa, qb = victim[0], victim[1], victim[2]
    moved = (d, tuple(sorted(((qa[0] + 7, qa[1] + 11, qa[2]), qb)))[0],
             tuple(sorted(((qa[0] + 7, qa[1] + 11, qa[2]), qb)))[1])
    mutated[moved] += 1
    counts_still_agree = sum(mutated.values()) == sum(fresh.values())
    identity_catches = bool((fresh - mutated) or (mutated - fresh))
    return {"ran": True,
            "counts_still_agree_after_mutation": counts_still_agree,
            "identity_comparison_catches_it": identity_catches,
            "caught": bool(counts_still_agree and identity_catches),
            "what": ("one contact moved to a different pair of quads; the "
                     "totals are unchanged, so a count-only binding would "
                     "still pass")}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", default="out/benchmark")
    ap.add_argument("--work", default="out/benchmark-evidence/witness-recensus")
    ap.add_argument("--indices", nargs="+",
                    default=["out/release/index.json",
                             "out/release/expansion_index.json"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    entries: dict = {}
    for p in a.indices:
        entries.update(json.loads(Path(p).read_text())["segments"])
    bench = json.loads((Path(a.benchmark) / "index.json").read_text())
    scored = [t for t in bench["traces"] if t.get("scored")]
    if a.limit:
        scored = scored[:a.limit]

    work = Path(a.work)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    rows, t0 = [], time.time()
    for i, t in enumerate(scored, 1):
        seg = t["segment"]
        e = entries[seg]
        td = Path(a.benchmark) / "traces" / seg
        r: dict = {"segment": seg}
        try:
            res = census_one(Path(e["input_mesh"]), seg[:40], EXCLUDE, CELL,
                             THREADS, MAXEDGE, work=work)
        except Exception as exc:                       # noqa: BLE001
            r.update(ok=False, error=f"recensus raised: {exc!r}")
            rows.append(r)
            print(f"[{i:3d}] {seg:52s} RECENSUS FAILED")
            continue
        if res is None:
            r.update(ok=False, error="the engine declined this input")
            rows.append(r)
            print(f"[{i:3d}] {seg:52s} DECLINED")
            continue

        fresh: collections.Counter = collections.Counter()
        for d in (0, 1):
            fresh += multiset_from_csv(Path(res[f"csv_d{d}"]), d)
        shipped = multiset_from_witness(td / "witness.npz")

        only_fresh = fresh - shipped
        only_shipped = shipped - fresh
        r.update(
            ok=not only_fresh and not only_shipped,
            n_fresh=int(sum(fresh.values())),
            n_shipped=int(sum(shipped.values())),
            counts_agree=sum(fresh.values()) == sum(shipped.values()),
            n_only_in_fresh=int(sum(only_fresh.values())),
            n_only_in_witness=int(sum(only_shipped.values())),
            example_only_in_fresh=[list(map(list, k[1:]))
                                   for k in list(only_fresh)[:2]],
            recensus_counts={f"d{d}": res[f"d{d}"]["transverse"]
                             for d in (0, 1)},
        )
        # Run the count-preserving mutation once, on the first trace that
        # actually has contacts to mutate.
        if shipped and not any(x.get("negative_control") for x in rows):
            r["negative_control"] = negative_control(fresh, shipped)
        rows.append(r)
        print(f"[{i:3d}] {seg:52s} "
              f"{'identical' if r['ok'] else 'DIFFERS'} "
              f"({r['n_shipped']} contacts)", flush=True)
        for d in (0, 1):
            Path(res[f"csv_d{d}"]).unlink(missing_ok=True)

    bad = [r for r in rows if not r.get("ok")]
    # The point of the exercise: did identity catch anything counts missed?
    count_only = [r for r in bad if r.get("counts_agree")]
    out = {
        "test": "witness bound to a fresh census by canonical identity",
        "identity": ("(diagonal, sorted[(v1,u1,tri1), (v2,u2,tri2)]) as a "
                     "multiset"),
        "census_parameters": {"exclude": EXCLUDE, "cell": CELL,
                              "maxedge": MAXEDGE, "touch_tol": 1e-3,
                              "diagonals": [0, 1]},
        "engine": "engines/selfcross",
        "n_traces": len(rows),
        "n_identical": len(rows) - len(bad),
        "n_differ": len(bad),
        "n_differ_while_counts_agree": len(count_only),
        "note": ("n_differ_while_counts_agree is the number of traces the "
                 "old count-only binding would have passed and this one "
                 "does not"),
        "seconds": round(time.time() - t0, 1),
        "negative_control": next((r["negative_control"] for r in rows
                                  if r.get("negative_control")), None),
        "failures": bad[:20],
        "traces": rows,
    }
    if a.out:
        Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    nc = out.get("negative_control") or {}
    print(f"\n{out['n_identical']}/{out['n_traces']} witnesses are "
          f"identical to a fresh census, by identity")
    if nc.get("ran"):
        print(f"negative control: a count-preserving identity change is "
              f"{'CAUGHT' if nc.get('caught') else 'MISSED'}"
              f" (counts still agree: {nc.get('counts_still_agree_after_mutation')})")
        if not nc.get("caught"):
            return 1
    if count_only:
        print(f"{len(count_only)} of the differences have MATCHING COUNTS -- "
              "the previous binding would have accepted them")
    for r in bad[:8]:
        print(f"  {r['segment']}: {r.get('error') or r}")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
