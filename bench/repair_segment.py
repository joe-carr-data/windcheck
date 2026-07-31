"""Slice-1 corpus executor: certified sub-voxel repair of any segment whose
sole crossing event passes the round-18 gate.

Derived from the w094 demonstrator (kept intact as the historical artifact);
parameterized over corpus, volume and resolution. One segment per
invocation:

    uv run python bench/repair_segment.py --segment <name-substring>

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

import argparse
import re

CORPORA = [
    ("Scroll 1", "data/scroll1_tifxyz", "20230205180739", "out/crossing_s1"),
    ("Scroll 5", "data/scroll5_tifxyz", "20241024131839", "out/crossing"),
    ("PHerc0139", "data/PHerc0139_tifxyz", "20250728140407", "out/crossing_0139"),
    ("PHerc0814", "data/PHerc0814_tifxyz", "20250804134230", "out/crossing_0814"),
    ("PHerc1667", "data/PHerc1667_tifxyz", "20231117161658", "out/crossing_1667"),
]
RES_UM = re.compile(r"-(\d+\.?\d*)um\.tifxyz$")
OUT = Path("out/repaired")
# Round-19 tier semantics: the strict tier binds the MAXIMUM EMITTED-VERTEX
# displacement (<=1 source voxel, measured after float32 emission), which is
# the CT-fidelity quantity. Relative branch motion is a search bound derived
# from the field: t_rel * max|w| <= 1, so symmetric mode (max|w|=0.5) may
# search to 2 vx relative. Both numbers go on the certificate.
STRICT_POINT_VX = 1.0
SEARCH_REL_VX = 2.0             # = STRICT_POINT_VX / 0.5 (symmetric field)
FIELD_MAX_COEFF = 0.5
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


def contact_multiset(csv: Path):
    """Multiset over TRIANGLE identities (schema v2, round 19): key =
    (v1,u1,tri1,v2,u2,tri2,verdict), canonical order guaranteed by the
    engine. Legacy 7-column CSVs are REFUSED for transactional acceptance
    — quad-level comparison hides one combo replacing another."""
    from collections import Counter
    lines = csv.read_text().splitlines()
    if lines and "tri1" not in lines[0]:
        raise ValueError(f"legacy census CSV (no triangle identities): {csv}")
    c = Counter()
    for line in lines[1:]:
        p = line.split(",")
        if p[4] in ("transverse", "grazing", "coplanar"):
            c[(int(p[0]), int(p[1]), int(p[7]),
               int(p[2]), int(p[3]), int(p[8]), p[4])] += 1
    return c


def is_target(key, ev) -> bool:
    q1, q2 = (key[0], key[1]), (key[3], key[4])
    A = set(map(tuple, ev["region_a"]))
    B = set(map(tuple, ev["region_b"]))
    return (q1 in A and q2 in B) or (q1 in B and q2 in A)


def accepted(before, after, ev) -> tuple[bool, str]:
    """Round-19 multiset-relative acceptance: every target key absent under
    every verdict; no key's multiplicity increases (which also bans any new
    key of any verdict); unrelated pre-existing contacts may remain."""
    for k, n in after.items():
        if is_target(k, ev):
            return False, f"target key survives: {k} x{n}"
        if n > before.get(k, 0):
            return False, f"new/increased contact: {k} {before.get(k,0)}->{n}"
    return True, "clean"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--segment", required=True,
                    help="unique substring of the segment directory name")
    args = ap.parse_args()
    exact, subs = [], []
    for corpus, root, volume, work in CORPORA:
        for d in Path(root).iterdir() if Path(root).exists() else []:
            if d.is_dir() and args.segment in d.name:
                (exact if d.name == args.segment else subs).append(
                    (corpus, d, volume, Path(work)))
    cands = exact or subs
    assert len(cands) == 1, \
        f"{args.segment}: {len(cands)} matches {[c[1].name for c in cands]}"
    hit = cands[0]
    corpus, seg, VOLUME, WORK = hit
    mesh = sorted(seg.glob(f"mesh/*{VOLUME}*.tifxyz"))[0]
    m = RES_UM.search(mesh.name)
    voxel_um = float(m.group(1)) if m else 7.91
    s = tifxyz.read(mesh)
    nv, nu = s.shape
    print(f"{corpus} {seg.name}  grid {nv}x{nu}  voxel {voxel_um} um")

    evs = {}
    for diag in (0, 1):
        rec = clip_rec(load_pairs(WORK / f"{seg.name[:40]}_d{diag}.csv"),
                       nv, nu)
        cand = [e for e in oriented_events(rec) if not e["ambiguous"]]
        assert len(cand) == 1, f"slice-1 gate: exactly one d{diag} event"
        e = cand[0]
        assert not e["self_touching"], "slice-1 gate: self-touching"
        evs[diag] = cand[0]
        print(f"  d{diag}: {len(rec)} pairs, regions "
              f"{sorted(cand[0]['region_a'])} x "
              f"{sorted(cand[0]['region_b'])}")
    union = merge_events(evs[0], evs[1])

    g = SurfaceGraph(s.points, s.valid, 0)
    gv, gu = np.nonzero(g.idx >= 0)
    X = np.empty((g.n, 3))
    X[g.idx[gv, gu]] = np.asarray(s.points, np.float64)[gv, gu]
    r = event_clearance(g, X, union, voxel_um, t_max_vx=SEARCH_REL_VX)
    if r is not None:
        print(f"  candidate: side {r['side']}  {r['t_vx']:.6f} vx "
              f"({r['t_mm']*1000:.2f} um)  d0 proposal")
        ev = union if r["side"] == "a" else \
            {"region_a": union["region_b"], "region_b": union["region_a"]}
    else:
        # round-20: the rigid-with-obstacles proposal finding nothing is not
        # the end -- the ranked search with the actual-field pre-check is
        print("  no rigid proposal; going straight to ranked search")
        ev = union
        r = {"side": "a", "direction": None, "t_vx": None}

    out = None if r["direction"] is None else \
        certified_repair(g, s.points, s.valid, ev, r["direction"],
                         r["t_vx"], mode="symmetric",
                         budget_vx=SEARCH_REL_VX,
                         support_vx=SUPPORT_VX)
    from windcheck.repair import search_repair
    ranked = None
    attempts = []                     # every candidate tried, with outcome
    if out is not None:
        P32, rep = out
        ranked = [{"direction": r["direction"], "lam": 0.5,
                   "t_rel": rep["applied_relative_vx"], "P32": P32,
                   "field_report": rep, "max_point_vx": rep["max_disp_vx"],
                   "provenance": "primary_certified"}]
    else:
        # round-20 fallback: ranked joint direction x split search with the
        # actual-field local pre-check; candidates go to the oracle in order
        print("  primary candidate failed; running ranked search")
        ranked = search_repair(g, s.points, s.valid, ev,
                               budget_point_vx=STRICT_POINT_VX,
                               support_vx=SUPPORT_VX)
        assert ranked, "no locally-clean candidate in the searched family"
        c = ranked[0]
        P32 = c["P32"]
        from windcheck.repair import displacement_stats
        rep = displacement_stats(s.points, P32, s.valid, 60.0)
        rep.update(c["field_report"])
        rep.update({"applied_relative_vx": c["t_rel"],
                    "candidate_lp_exit_vx": None,
                    "search_lambda": c["lam"],
                    "search_direction": c["direction"]})
        rep.setdefault("final_max_disp_vx", rep["max_disp_vx"])
        c["provenance"] = "search"
        r = {"side": "a", "direction": c["direction"],
             "t_vx": c["t_rel"], "t_mm": c["t_rel"] * voxel_um / 1000.0}
    print(f"  local both-diagonal clearance: applied "
          f"{rep['applied_relative_vx']:.6f} vx, max point "
          f"{rep['max_disp_vx']:.6f} vx ({rep['max_disp_vx']*voxel_um:.2f} um)")
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
    active = (ranked[0] if ranked else
              {"direction": r["direction"], "lam": 0.5,
               "provenance": "primary_certified"})
    t_rel = rep["applied_relative_vx"]
    t_dirty = 0.0
    dst = OUT / f"{seg.name}_repaired.tifxyz"
    tag = hashlib.sha256(seg.name.encode()).hexdigest()[:12]
    wdir = OUT / f"work_{tag}"                  # round 19: parallel-safe,
                                                # immune to name truncation
    before = census_one(mesh, f"{tag}_before", 1, 40.0, 0, 60.0, wdir)
    before_ms = {d: contact_multiset(wdir / f"{tag}_before_d{d}.csv")
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
        if t_rel > SEARCH_REL_VX:
            print("  BUDGET EXHAUSTED before the engine oracle was clean")
            sys.exit(1)
        if attempt > 0:
            P2, _ = apply_field(g, s.points, s.valid, ev,
                                active["direction"], t_rel,
                                active["lam"], SUPPORT_VX)
            P32 = P2.astype(np.float32).astype(np.float64)
        emit(P32)
        s2 = tifxyz.read(dst)
        assert s2.shape == s.shape and (s2.valid == s.valid).all()
        after = census_one(dst, f"{tag}_after", 1, 40.0, 0, 60.0, wdir)
        after_ms = {d: contact_multiset(wdir / f"{tag}_after_d{d}.csv")
                    for d in (0, 1)}
        verdicts = [accepted(before_ms[d], after_ms[d], ev) for d in (0, 1)]
        print(f"  engine oracle @ {t_rel:.4f} vx: "
              f"d0 {verdicts[0][1]} | d1 {verdicts[1][1]}")
        if all(v[0] for v in verdicts):
            break
        # round-20 direction fallback: on oracle failure, prefer switching
        # to the next locally-clean candidate over growing displacement
        # (growth along a direction that induces contacts cannot converge)
        if ranked is None or len(ranked) <= 1:
            print("  oracle failed; running ranked search for alternatives")
            more = search_repair(g, s.points, s.valid, ev,
                                 budget_point_vx=STRICT_POINT_VX,
                                 support_vx=SUPPORT_VX)
            ranked = (ranked or []) + [c for c in more
                                       if np.dot(c["direction"],
                                                 r["direction"]) < 0.999]
        if len(ranked) > 1:
            attempts.append({"direction": active["direction"],
                             "lam": active.get("lam", 0.5),
                             "t_rel": t_rel, "outcome": "oracle_failed"})
            ranked.pop(0)
            c = ranked[0]
            c.setdefault("provenance", "search")
            active = c
            P32 = c["P32"]
            r["direction"] = c["direction"]
            t_rel = c["t_rel"]
            t_dirty = 0.0
            print(f"  switching to candidate: lam {c['lam']} "
                  f"t {c['t_rel']:.3f} dir "
                  f"{[round(x,2) for x in c['direction']]}")
            continue
        t_dirty, t_rel = t_rel, t_rel + 0.08
    else:
        print("  engine oracle never clean"); sys.exit(1)

    # engine-driven bracket refinement (round 18: do not spend fidelity on
    # the coarse step) to a stated resolution
    RES = 0.01
    lo, hi = t_dirty, t_rel
    def engine_clean(t):
        P2, _ = apply_field(g, s.points, s.valid, ev, active["direction"],
                            t, active.get("lam", 0.5), SUPPORT_VX)
        Pq = P2.astype(np.float32).astype(np.float64)
        emit(Pq)
        a = census_one(dst, f"{tag}_after", 1, 40.0, 0, 60.0, wdir)
        oks = all(accepted(before_ms[d], contact_multiset(
            wdir / f"{tag}_after_d{d}.csv"), ev)[0] for d in (0, 1))
        return Pq, a, oks
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
    after_ms = {d: contact_multiset(wdir / f"{tag}_after_d{d}.csv")
                for d in (0, 1)}
    print(f"  engine-certified relative: {t_rel:.4f} vx "
          f"(bracket [{lo:.4f}, {hi:.4f}], resolution {RES})")
    rep2 = displacement_stats(s.points, P32, s.valid, 60.0)
    rep.update({f"final_{k}": v for k, v in rep2.items()})
    rep["engine_certified_relative_vx"] = t_rel
    rep["oracle_bracket_vx"] = [lo, hi]
    rep["oracle_resolution_vx"] = RES
    rep["winning_candidate"] = {"direction": active["direction"],
                                "lam": active.get("lam", 0.5),
                                "provenance": active.get("provenance", "?")}
    rep["attempted_candidates"] = attempts
    # gates re-enforced from the FINAL emitted mesh (round 21): a direction
    # switch must never inherit the primary's gate results
    for kk in ("final_quads_newly_dropped", "final_quads_newly_retained",
               "final_triangle_inversions_d0", "final_triangle_inversions_d1"):
        assert rep[kk] == 0, f"final gate failed: {kk}={rep[kk]}"
    assert rep["final_max_disp_vx"] <= STRICT_POINT_VX, "point tier exceeded"
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
    pairs_before = {d: canonical_pairs(wdir / f"{tag}_before_d{d}.csv")
                    for d in (0, 1)}
    pairs_after = {d: canonical_pairs(wdir / f"{tag}_after_d{d}.csv")
                   for d in (0, 1)}
    ok = True
    for d in (0, 1):
        b, a = before[f"d{d}"], after[f"d{d}"]
        print(f"  global d{d}: transverse {b['transverse']} -> "
              f"{a['transverse']}, coplanar {b['coplanar']} -> "
              f"{a['coplanar']}, grazing {b['grazing']} -> {a['grazing']}, "
              f"quads_dropped {b['quads_dropped']} -> {a['quads_dropped']}")
        ok &= a["quads_dropped"] == b["quads_dropped"]
    ok &= all(accepted(before_ms[d], after_ms[d], ev)[0] for d in (0, 1))
    # THE strict-tier gate: max emitted-vertex displacement (round 19)
    ok &= rep["final_max_disp_vx"] <= STRICT_POINT_VX

    porcelain = subprocess.run(["git", "status", "--porcelain"],
                               capture_output=True, text=True).stdout
    compiler = subprocess.run(["clang++", "--version"], capture_output=True,
                              text=True).stdout.splitlines()[0]
    cert = {
        "claim": ((f"The emitted {seg.name} mesh is globally recensus-clean for "
                   "non-adjacent transverse, coplanar, and grazing contacts "
                   "under both triangulations and the stated validator "
                   "parameters, after moving "
                   f"{rep['final_quantized_moved_vertices']} vertices by at "
                   f"most {rep['final_max_disp_vx']:.3f} voxel. This is an "
                   "endpoint certificate; texture fidelity and "
                   "collision-free deformation through time are not "
                   "certified.") if ok else "REJECTED"),
        "quality_metrics": qual,
        "corpus": corpus, "segment": seg.name, "volume": VOLUME,
        "voxel_um": voxel_um,
        "mode": "symmetric",
        "fidelity_tier": "strict_1vx_max_emitted_vertex_displacement",
        "tier_binds": "max point displacement after float32 emission",
        "field_max_coeff": FIELD_MAX_COEFF,
        "search_relative_cap_vx": SEARCH_REL_VX,
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
            "census_csvs": {f"{w}_d{d}": sha(wdir / f"{tag}_{w}_d{d}.csv")
                            for w in ("before", "after") for d in (0, 1)},
        },
        "code_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True,
            text=True).stdout.strip(),
        "git_status_porcelain": porcelain,
        "executor_sha256": sha(Path("bench/repair_segment.py")),
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
