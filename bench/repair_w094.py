"""The w094 demonstrator: one certified sub-voxel repair, end to end.

Round-17 locked spec: symmetric harmonic repair of w094's single crossing
event (union of its d0 and d1 regions), strict one-voxel fidelity tier,
per-triangle swept-AABB obstacle discovery, tifxyz emission, reload, global
C++ recensus under BOTH diagonals, and a certificate that separates "local
both-diagonal clearance" (Python predicate) from "globally validated"
(reloaded C++ census). Acceptance for w094 is ZERO global transverse
crossings under both diagonals -- every original crossing belongs to the
target event.

    uv run python bench/repair_w094.py
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, "src")
sys.path.insert(0, "bench")
from windcheck import tifxyz                                    # noqa: E402
from windcheck.check import load_pairs                          # noqa: E402
from windcheck.intrinsic import SurfaceGraph, oriented_events   # noqa: E402
from windcheck.repair import (certified_repair, event_clearance,  # noqa: E402
                              merge_events)
from crossing_census import census_one                          # noqa: E402

VOLUME = "20241024131839"
WORK = Path("out/crossing")
OUT = Path("out/repaired")
BUDGET_STRICT_VX = 1.0          # 7.91 um -- the strict tier
SUPPORT_VX = 40.0


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def clip_rec(rec, nv, nu):
    return rec[(rec["v1"] < nv - 1) & (rec["v2"] < nv - 1)
               & (rec["u1"] < nu - 1) & (rec["u2"] < nu - 1)]


def canonical_pairs(csv: Path) -> list:
    rec = load_pairs(csv)
    return sorted((int(r["v1"]), int(r["u1"]), int(r["v2"]), int(r["u2"]))
                  for r in rec)


def contact_keys(csv: Path) -> set:
    """Canonical keys of grazing/coplanar rows (round 17: a repair must not
    convert a crossing into a touch)."""
    out = set()
    for line in csv.read_text().splitlines()[1:]:
        p = line.split(",")
        if p[4] in ("grazing", "coplanar"):
            out.add((int(p[0]), int(p[1]), int(p[2]), int(p[3])))
    return out


def main() -> None:
    seg = next(Path("data/scroll5_tifxyz").glob("*w094*"))
    mesh = sorted(seg.glob(f"mesh/*{VOLUME}*.tifxyz"))[0]
    s = tifxyz.read(mesh)
    nv, nu = s.shape
    print(f"w094: {seg.name}  grid {nv}x{nu}")

    evs = {}
    for diag in (0, 1):
        rec = clip_rec(load_pairs(WORK / f"{seg.name[:40]}_d{diag}.csv"),
                       nv, nu)
        cand = [e for e in oriented_events(rec) if not e["ambiguous"]]
        assert len(cand) == 1, f"w094 must have exactly one d{diag} event"
        evs[diag] = cand[0]
        print(f"  d{diag}: {len(rec)} pairs, regions "
              f"{sorted(cand[0]['region_a'])} x "
              f"{sorted(cand[0]['region_b'])}")
    union = merge_events(evs[0], evs[1])

    g = SurfaceGraph(s.points, s.valid, 0)
    gv, gu = np.nonzero(g.idx >= 0)
    X = np.empty((g.n, 3))
    X[g.idx[gv, gu]] = np.asarray(s.points, np.float64)[gv, gu]
    r = event_clearance(g, X, union, 7.91, t_max_vx=BUDGET_STRICT_VX)
    assert r is not None, "no rigid candidate within the strict tier"
    print(f"  candidate: side {r['side']}  {r['t_vx']:.6f} vx "
          f"({r['t_mm']*1000:.2f} um)  d0 proposal")
    ev = union if r["side"] == "a" else \
        {"region_a": union["region_b"], "region_b": union["region_a"]}

    out = certified_repair(g, s.points, s.valid, ev, r["direction"],
                           r["t_vx"], mode="symmetric",
                           budget_vx=BUDGET_STRICT_VX,
                           support_vx=SUPPORT_VX)
    assert out is not None, "certification failed within the strict tier"
    P32, rep = out
    print(f"  local both-diagonal clearance: applied "
          f"{rep['applied_relative_vx']:.6f} vx, max point "
          f"{rep['max_disp_vx']:.6f} vx ({rep['max_disp_vx']*7.91:.2f} um)")
    for k in ("quads_newly_dropped", "quads_newly_retained",
              "triangle_inversions_d0", "triangle_inversions_d1"):
        assert rep[k] == 0, f"quality gate failed: {k}={rep[k]}"

    # ---- obstacle discovery: per-triangle swept AABBs (round 17 Q4)
    P0 = np.asarray(s.points, np.float64)
    qv, qu = np.nonzero(g.Q)
    C0 = np.stack([P0[qv, qu], P0[qv + 1, qu], P0[qv, qu + 1],
                   P0[qv + 1, qu + 1]])
    qlo, qhi = C0.min(axis=0), C0.max(axis=0)

    def discover(P_end):
        delta = np.linalg.norm(P_end - P0, axis=-1)
        mv, mu = np.nonzero(delta > 0)
        inc = {(v + dv, u + du) for v, u in zip(mv.tolist(), mu.tolist())
               for dv in (-1, 0) for du in (-1, 0)}
        inc = {(v, u) for v, u in inc
               if 0 <= v < nv - 1 and 0 <= u < nu - 1 and g.Q[v, u]}
        hits = set()
        for v, u in inc:
            c0 = np.array([P0[v, u], P0[v + 1, u], P0[v, u + 1],
                           P0[v + 1, u + 1]])
            c1 = np.array([P_end[v, u], P_end[v + 1, u], P_end[v, u + 1],
                           P_end[v + 1, u + 1]])
            blo = np.minimum(c0.min(0), c1.min(0))
            bhi = np.maximum(c0.max(0), c1.max(0))
            m = np.all(qhi >= blo, axis=1) & np.all(qlo <= bhi, axis=1)
            for a, b in zip(qv[m].tolist(), qu[m].tolist()):
                if max(abs(a - v), abs(b - u)) > 1 and (a, b) not in inc:
                    hits.add((a, b))
        hits -= set(map(tuple, ev["region_a"]))             | set(map(tuple, ev["region_b"]))
        return inc, hits

    inc, hits = discover(P32)
    print(f"  obstacle discovery (preliminary): {len(inc)} swept quads, "
          f"{len(hits)} third-party hits")

    # ---- engine-oracle loop: the local predicate is blind to the census's
    # scale-aware grazing tolerance (~0.025 vx here); grow the applied
    # displacement until the C++ census itself is clean (round 17: a
    # crossing must not become a touch)
    from windcheck.repair import apply_field, displacement_stats
    t_rel = rep["applied_relative_vx"]
    t_dirty = 0.0
    dst = OUT / f"{seg.name}_repaired.tifxyz"
    before = census_one(mesh, "w094_before", 1, 40.0, 0, 60.0, OUT / "work")
    contacts_before = {d: contact_keys(OUT / "work" / f"w094_before_d{d}.csv")
                       for d in (0, 1)}

    def emit(P_arr):
        if dst.exists():
            shutil.rmtree(dst)
        dst.mkdir(parents=True)
        Q32 = P_arr.astype(np.float32)
        keep = ~s.valid
        bands = []
        for i, ax in enumerate(("x", "y", "z")):
            a = Q32[..., i].copy()
            a[keep] = np.asarray(s.points, np.float32)[keep, i]
            tifffile.imwrite(dst / f"{ax}.tif", a)
            bands.append(a)
        for extra in ("mask.tif", "mask.png"):
            srcf = mesh / extra
            if srcf.exists():
                shutil.copy(srcf, dst / extra)
        # Displacement moves coordinates; the source bbox no longer describes
        # what is written here, and consumers filter on it.
        tifxyz.write_meta(mesh, dst, np.stack(bands, axis=-1), s.valid)

    after = None
    for attempt in range(12):
        if t_rel > BUDGET_STRICT_VX:
            print("  BUDGET EXHAUSTED before the engine oracle was clean")
            sys.exit(1)
        if attempt > 0:
            P2, _ = apply_field(g, s.points, s.valid, ev, r["direction"],
                                t_rel, "symmetric", SUPPORT_VX)
            P32 = P2.astype(np.float32).astype(np.float64)
        emit(P32)
        s2 = tifxyz.read(dst)
        assert s2.shape == s.shape and (s2.valid == s.valid).all()
        after = census_one(dst, "w094_after", 1, 40.0, 0, 60.0, OUT / "work")
        contacts_after = {d: contact_keys(OUT / "work" / f"w094_after_d{d}.csv")
                          for d in (0, 1)}
        trans = sum(after[f"d{d}"]["transverse"] for d in (0, 1))
        n_con = sum(len(contacts_after[d]) for d in (0, 1))
        print(f"  engine oracle @ {t_rel:.4f} vx: transverse {trans}, "
              f"contacts {n_con}")
        if trans == 0 and n_con == 0:      # round 18: ALL contacts empty
            break
        t_dirty, t_rel = t_rel, t_rel + 0.08
    else:
        print("  engine oracle never clean"); sys.exit(1)

    # engine-driven bracket refinement (round 18: do not spend fidelity on
    # the coarse step) to a stated resolution
    RES = 0.01
    lo, hi = t_dirty, t_rel
    def engine_clean(t):
        P2, _ = apply_field(g, s.points, s.valid, ev, r["direction"],
                            t, "symmetric", SUPPORT_VX)
        Pq = P2.astype(np.float32).astype(np.float64)
        emit(Pq)
        a = census_one(dst, "w094_after", 1, 40.0, 0, 60.0, OUT / "work")
        c = sum(len(contact_keys(OUT / "work" / f"w094_after_d{d}.csv"))
                for d in (0, 1))
        t_ = sum(a[f"d{d}"]["transverse"] for d in (0, 1))
        return Pq, a, (t_ == 0 and c == 0)
    while hi - lo > RES:
        mid = 0.5 * (lo + hi)
        _, _, okm = engine_clean(mid)
        if okm:
            hi = mid
        else:
            lo = mid
    t_rel = hi
    P32, after, okf = engine_clean(t_rel)
    assert okf
    contacts_after = {d: contact_keys(OUT / "work" / f"w094_after_d{d}.csv")
                      for d in (0, 1)}
    print(f"  engine-certified relative: {t_rel:.4f} vx "
          f"(bracket [{lo:.4f}, {hi:.4f}], resolution {RES})")
    rep2 = displacement_stats(s.points, P32, s.valid, 60.0)
    rep.update({f"final_{k}": v for k, v in rep2.items()})
    rep["engine_certified_relative_vx"] = t_rel
    rep["oracle_bracket_vx"] = [lo, hi]
    rep["oracle_resolution_vx"] = RES
    inc, hits = discover(P32)              # round 18: certify FINAL geometry
    print(f"  obstacle discovery (final): {len(inc)} swept quads, "
          f"{len(hits)} third-party hits")

    # deformation-quality metrics over affected quads (round 18: disclose)
    qual = {}
    for diag_ in (0, 1):
        areas0, areas1, rots, easp0, easp1 = [], [], [], [], []
        for v, u in inc:
            combos = (((v, u), (v, u + 1), (v + 1, u + 1)),
                      ((v, u), (v + 1, u + 1), (v + 1, u))) if diag_ == 0                 else (((v, u), (v, u + 1), (v + 1, u)),
                      ((v, u + 1), (v + 1, u + 1), (v + 1, u)))
            for c in combos:
                n0 = np.cross(P0[c[1]] - P0[c[0]], P0[c[2]] - P0[c[0]])
                n1 = np.cross(P32[c[1]] - P32[c[0]], P32[c[2]] - P32[c[0]])
                a0, a1 = np.linalg.norm(n0) / 2, np.linalg.norm(n1) / 2
                areas0.append(a0); areas1.append(a1)
                cosang = float(n0 @ n1) / max(np.linalg.norm(n0)
                                              * np.linalg.norm(n1), 1e-30)
                rots.append(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
        ratio = np.array(areas1) / np.maximum(areas0, 1e-30)
        qual[f"d{diag_}"] = {
            "area_ratio_min": float(ratio.min()),
            "area_ratio_max": float(ratio.max()),
            "max_normal_rotation_deg": float(np.max(rots)),
        }
    edges0, edges1 = [], []
    for v, u in inc:
        for a, b in (((v, u), (v, u + 1)), ((v, u), (v + 1, u)),
                     ((v, u), (v + 1, u + 1))):
            edges0.append(np.linalg.norm(P0[a] - P0[b]))
            edges1.append(np.linalg.norm(P32[a] - P32[b]))
    qual["max_edge_length_change_vx"] = float(
        np.max(np.abs(np.array(edges1) - np.array(edges0))))
    tw0 = [np.linalg.norm(P0[v, u] - P0[v + 1, u] - P0[v, u + 1]
                          + P0[v + 1, u + 1]) for v, u in inc]
    tw1 = [np.linalg.norm(P32[v, u] - P32[v + 1, u] - P32[v, u + 1]
                          + P32[v + 1, u + 1]) for v, u in inc]
    qual["twist_median_before_after_vx"] = [float(np.median(tw0)),
                                            float(np.median(tw1))]
    qual["twist_max_before_after_vx"] = [float(np.max(tw0)),
                                         float(np.max(tw1))]
    print(f"  quality: max normal rotation d0 "
          f"{qual['d0']['max_normal_rotation_deg']:.1f} deg, d1 "
          f"{qual['d1']['max_normal_rotation_deg']:.1f} deg")
    pairs_before = {d: canonical_pairs(OUT / "work" / f"w094_before_d{d}.csv")
                    for d in (0, 1)}
    pairs_after = {d: canonical_pairs(OUT / "work" / f"w094_after_d{d}.csv")
                   for d in (0, 1)}
    ok = True
    for d in (0, 1):
        b, a = before[f"d{d}"], after[f"d{d}"]
        print(f"  global d{d}: transverse {b['transverse']} -> "
              f"{a['transverse']}, coplanar {b['coplanar']} -> "
              f"{a['coplanar']}, grazing {b['grazing']} -> {a['grazing']}, "
              f"quads_dropped {b['quads_dropped']} -> {a['quads_dropped']}")
        ok &= a["transverse"] == 0
        ok &= a["quads_dropped"] == b["quads_dropped"]
    ok &= all(not (contacts_after[d] - contacts_before[d]) for d in (0, 1))

    porcelain = subprocess.run(["git", "status", "--porcelain"],
                               capture_output=True, text=True).stdout
    compiler = subprocess.run(["clang++", "--version"], capture_output=True,
                              text=True).stdout.splitlines()[0]
    cert = {
        "claim": (("The emitted w094 mesh is globally recensus-clean for "
                   "non-adjacent transverse, coplanar, and grazing contacts "
                   "under both triangulations and the stated validator "
                   "parameters, after moving "
                   f"{rep['final_quantized_moved_vertices']} vertices by at "
                   f"most {rep['final_max_disp_vx']:.3f} voxel. This is an "
                   "endpoint certificate; texture fidelity and "
                   "collision-free deformation through time are not "
                   "certified.") if ok else "REJECTED"),
        "quality_metrics": qual,
        "segment": seg.name, "volume": VOLUME, "voxel_um": 7.91,
        "mode": "symmetric", "budget_tier": "strict_1_voxel",
        "support_vx": SUPPORT_VX,
        "event": {"region_a": sorted(map(list, ev["region_a"])),
                  "region_b": sorted(map(list, ev["region_b"])),
                  "candidate_side": r["side"],
                  "direction": r["direction"]},
        "repair": rep,
        "obstacle_discovery": {"swept_quads": len(inc),
                               "third_party_hits": sorted(map(list, hits))},
        "census_before": {f"d{d}": before[f"d{d}"] for d in (0, 1)},
        "census_after": {f"d{d}": after[f"d{d}"] for d in (0, 1)},
        "pairs_before": {str(d): pairs_before[d] for d in (0, 1)},
        "pairs_after": {str(d): pairs_after[d] for d in (0, 1)},
        "hashes": {
            "input": {f"{ax}.tif": sha(mesh / f"{ax}.tif")
                      for ax in ("x", "y", "z")},
            "output": {f"{ax}.tif": sha(dst / f"{ax}.tif")
                       for ax in ("x", "y", "z")},
            "engine_binary": sha(Path("engines/selfcross")),
            "engine_source": sha(Path("engines/selfcross.cpp")),
            "census_csvs": {f"{w}_d{d}": sha(OUT / "work" / f"w094_{w}_d{d}.csv")
                            for w in ("before", "after") for d in (0, 1)},
        },
        "code_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True,
            text=True).stdout.strip(),
        "git_status_porcelain": porcelain,
        "executor_sha256": sha(Path("bench/repair_w094.py")),
        "uv_lock_sha256": sha(Path("uv.lock")),
        "compiler": compiler,
        "meta_json_sha256": sha(mesh / "meta.json")
            if (mesh / "meta.json").exists() else None,
        "params": {"cell": 40.0, "exclude": 1, "maxedge": 60.0,
                   "touch_tol": 1e-3},
        "note": ("local clearance is the Python both-diagonal predicate; "
                 "global validation is the reloaded C++ census. Endpoint "
                 "proof only: no collision-free-path claim."),
    }
    cpath = OUT / f"{seg.name}_certificate.json"
    cpath.write_text(json.dumps(cert, indent=1))
    print(f"\n{'ACCEPTED' if ok else 'REJECTED'}: certificate {cpath}")


if __name__ == "__main__":
    main()
