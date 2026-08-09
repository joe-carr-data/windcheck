"""Verify a packaged benchmark INDEPENDENTLY OF THE PACKER.

Not "from the shipped files alone": it also needs the referenced
operands, the release indices and the certificates. What it does not need
is `pack_benchmark`, and that is the point -- a checker that imports the
thing it checks agrees with it by construction.

It re-derives, from `masks.npz` and the referenced operands:

  1. every recorded sha256 in `index.json` matches the file on disk;
  2. `removed | retained` is exactly the input's valid set, and the two
     do not overlap;
  3. the reference derivative is reconstructible: input with `x = y = z =
     -1` stamped at `removed` reproduces it over the COMPLETE arrays,
     shape and dtype included -- not merely at the cells one of them
     calls valid;
  4. the witness contact count equals the certificate's `census_before`;
  5. every witness contact names quads that are non-adjacent under the
     frozen `exclude = 1` rule, and both quads exist in the input.

Deliberately does NOT re-use pack_benchmark: a checker that imports the
thing it checks agrees with it by construction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import tifffile

EXCLUDE = 1
MAXEDGE = 60.0


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def planes(d: Path):
    x, y, z = (np.asarray(tifffile.imread(d / f"{a}.tif"), np.float32)
               for a in "xyz")
    valid = ~((x == -1) & (y == -1) & (z == -1))
    valid &= np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    m = d / "mask.tif"
    if m.is_file():
        mm = np.asarray(tifffile.imread(m))
        if mm.ndim == 3:
            mm = mm[..., 0]
        if mm.shape == valid.shape:
            valid &= mm >= 255
    return x, y, z, valid


def check_trace(td: Path, row: dict, indices: dict) -> list[str]:
    bad: list[str] = []
    seg = row["segment"]

    for name, meta in (row.get("files") or {}).items():
        f = td / name
        if not f.is_file():
            bad.append(f"{name}: missing")
            continue
        if sha256_file(f) != meta["sha256"]:
            bad.append(f"{name}: sha256 differs from the index")

    prov = json.loads((td / "provenance.json").read_text())
    e = indices[seg]
    inp = Path(e["input_mesh"])
    ref = Path(e.get("output_mesh") or e["input_mesh"])
    if not inp.is_dir():
        return bad + [f"input not on disk: {inp}"]

    xi, yi, zi, vi = planes(inp)
    xr, yr, zr, vr = planes(ref)
    m = np.load(td / "masks.npz")
    removed = m["removed"].astype(bool)
    retained = m["retained"].astype(bool)
    changed = m["changed"].astype(bool)

    if removed.shape != vi.shape:
        return bad + [f"masks grid {removed.shape} != input {vi.shape}"]
    if (removed & retained).any():
        bad.append("removed and retained overlap")
    if not np.array_equal(removed | retained, vi):
        bad.append("removed | retained is not the input's valid set")
    if changed.any():
        bad.append(f"{int(changed.sum())} retained cells have moved "
                   "coordinates, which the frozen policy forbids")

    for a, (si, sr) in zip("xyz", ((xi, xr), (yi, yr), (zi, zr))):
        if si.shape != sr.shape or si.dtype != sr.dtype:
            bad.append(f"{a}: shape/dtype differ between input and "
                       f"reference ({si.shape}/{si.dtype} vs "
                       f"{sr.shape}/{sr.dtype})")
            continue
        rebuilt = np.where(removed, si.dtype.type(-1.0), si)
        if not np.array_equal(rebuilt, sr):
            n = int((rebuilt != sr).sum())
            bad.append(f"{a}: the reference is not reconstructible from "
                       f"input + removed mask over the complete array "
                       f"({n} cells differ)")

    w = json.loads((td / "witness.json").read_text())
    cb = w.get("bound_to_certificate_census_before") or {}
    cert = json.loads(Path(e["certificate"]).read_text())
    if cb != (cert.get("census_before") or {}):
        bad.append("witness.json's recorded census_before is not the "
                   "certificate's")
    n_total = w["n_contacts_total"]
    want = sum((cert.get("census_before") or {}).get(f"d{d}", {})
               .get("transverse", 0) for d in (0, 1))
    if n_total != want:
        bad.append(f"witness carries {n_total} contacts, the certificate's "
                   f"census_before says {want}")

    if n_total:
        wp = w.get("payload") or {}
        f = td / wp.get("file", "witness.npz")
        if not f.is_file():
            bad.append("witness payload missing")
        else:
            if sha256_file(f) != wp.get("sha256"):
                bad.append("witness payload sha256 differs from witness.json")
            a = np.load(f)
            sizes = {k: int(a[k].size) for k in a.files}
            if len(set(sizes.values())) != 1:
                bad.append(f"witness arrays have unequal lengths: {sizes}")
            if a["v1"].size != n_total:
                bad.append(f"witness payload has {a['v1'].size} rows, "
                           f"witness.json says {n_total}")
            for k in ("diagonal", "tri1", "tri2"):
                vals = set(int(v) for v in np.unique(a[k]))
                if not vals <= {0, 1}:
                    bad.append(f"witness {k} carries {sorted(vals)[:5]}; "
                               "require exactly {0,1}")
            per = w.get("n_contacts_per_diagonal") or {}
            for d in (0, 1):
                got_n = int((a["diagonal"] == d).sum())
                if per.get(f"d{d}") is not None and got_n != per[f"d{d}"]:
                    bad.append(f"witness payload has {got_n} d{d} contacts, "
                               f"witness.json says {per[f'd{d}']}")
            for i, smp in enumerate(w.get("sample_first_50") or []):
                if (int(a["diagonal"][i]) != smp["diagonal"]
                        or [int(a["v1"][i]), int(a["u1"][i])] != smp["quad_a"]
                        or [int(a["v2"][i]), int(a["u2"][i])] != smp["quad_b"]
                        or int(a["tri1"][i]) != smp["tri_a"]
                        or int(a["tri2"][i]) != smp["tri_b"]):
                    bad.append(f"sample_first_50[{i}] does not match the "
                               "payload")
                    break
            nv, nu = vi.shape
            for k, lim in (("v1", nv - 1), ("v2", nv - 1),
                           ("u1", nu - 1), ("u2", nu - 1)):
                v = a[k]
                if v.size and (v.min() < 0 or v.max() >= lim):
                    bad.append(f"witness {k} out of range [0,{lim})")
            # exclude=1: contacts between quads within one lattice step of
            # each other are adjacency, not defect, and must not appear.
            near = ((np.abs(a["v1"].astype(np.int64) - a["v2"]) <= EXCLUDE)
                    & (np.abs(a["u1"].astype(np.int64) - a["u2"]) <= EXCLUDE))
            if near.any():
                bad.append(f"{int(near.sum())} witness contacts are between "
                           f"quads within exclude={EXCLUDE}; those are "
                           "adjacency, not self-intersection")
            # Both quads must be quads the census would actually build:
            # four valid corners AND all six pairwise corner distances
            # within maxedge. Checking only the ORIGIN vertex let a
            # witness name a quad the census never instantiated.
            P = np.stack([xi, yi, zi], axis=-1).astype(np.float64)
            Q = (vi[:-1, :-1] & vi[:-1, 1:] & vi[1:, :-1] & vi[1:, 1:])
            p00, p10, p01, p11 = (P[:-1, :-1], P[1:, :-1],
                                  P[:-1, 1:], P[1:, 1:])
            e = np.zeros(Q.shape)
            for c1, c2 in ((p00, p01), (p01, p11), (p11, p10), (p00, p10),
                           (p00, p11), (p10, p01)):
                e = np.maximum(e, np.linalg.norm(c1 - c2, axis=-1))
            Q = Q & (e <= MAXEDGE)
            for s in ("1", "2"):
                ok = Q[a[f"v{s}"], a[f"u{s}"]]
                if not ok.all():
                    bad.append(f"{int((~ok).sum())} witness contacts name a "
                               f"quad {s} the census would not have built "
                               "(invalid corners or over maxedge)")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", default="out/benchmark")
    ap.add_argument("--indices", nargs="+",
                    default=["out/release/index.json",
                             "out/release/expansion_index.json"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    root = Path(a.benchmark)
    idx = json.loads((root / "index.json").read_text())
    entries: dict = {}
    for p in a.indices:
        entries.update(json.loads(Path(p).read_text())["segments"])

    results, n = [], 0
    for row in idx["traces"]:
        if not row.get("scored"):
            continue
        if a.limit and n >= a.limit:
            break
        n += 1
        td = root / "traces" / row["segment"]
        bad = check_trace(td, row, entries)
        results.append({"segment": row["segment"], "problems": bad})
        print(f"[{n:3d}] {row['segment']:52s} "
              f"{'ok' if not bad else 'PROBLEMS: ' + bad[0]}", flush=True)

    failed = [r for r in results if r["problems"]]
    out = {"benchmark": str(root), "n_traces": len(results),
           "n_ok": len(results) - len(failed), "n_failed": len(failed),
           "independent_of": "bench/pack_benchmark.py",
           "not_self_contained": ("also reads the referenced operands, the "
                                  "release indices and the certificates"),
           "checks": ["recorded sha256s",
                      "removed|retained == input valid",
                      "reference reconstructible over the COMPLETE arrays",
                      "witness count bound to the certificate",
                      "witness arrays equal length, diagonals/triangles "
                      "in {0,1}, per-diagonal counts, sample matches "
                      "payload",
                      "witness contacts non-adjacent, in range, and on "
                      "quads the census would build"],
           "failures": failed[:20], "results": results}
    if a.out:
        Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"\n{out['n_ok']}/{out['n_traces']} traces verify independently "
          f"of the packer")
    for f in failed[:10]:
        print(f"  {f['segment']}: {f['problems'][0]}")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
