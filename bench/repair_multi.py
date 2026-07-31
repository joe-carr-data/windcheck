"""Slice-2 multi-event executor: transactional certified repair of every
crossing event in a segment, with mechanical dispositions.

Per PLAN-AUGUST W2 + round-22 amendments:
- d0/d1 events matched with match_events (matched / unmatched / ambiguous,
  scores recorded); merge_events only orients an already-matched pair.
- Min-exit PREFILTER on the CURRENT geometry: rigid search is skipped when
  any diagonal constituent's L_safe exceeds the admissible relative bound.
  Diagonal bounds are never averaged. Sweep/prefilter results are
  scheduling hints; final labels come only from the final accepted mesh.
- Transactional per-event loop, isolated-first then lowest-bound-first:
  accept -> the emitted mesh becomes the base and the WHOLE segment is
  recensused and rematched from it (never stale rows); fail -> rollback,
  run event_min_exit (cheap), label.
- Acceptance = round-19 triangle-identity multiset rule per diagonal:
  target keys absent under every verdict, no multiplicity increases,
  unrelated pre-existing contacts may remain or decrease (collateral
  resolution kept). Certified-infeasible events stay in recensus.
- Hard gates (pre-registered, any failure rejects the transaction):
  identical validity mask; zero newly dropped/retained quads; unchanged
  component count under both diagonals; zero triangle inversions; per-event
  and CUMULATIVE-vs-original max emitted-vertex displacement <= strict
  tier. Warning metrics (normal rotation, area ratio, edge/twist change)
  are recorded as continuous values, never invented pass thresholds.
- Instrumentation: wall-time split (graph builds / min-exit / search /
  oracle censuses / emit), candidates tried, fallback engagements.

    uv run python bench/repair_multi.py --segment <name-substring>
    uv run python bench/repair_multi.py --segment <name> --plan-only
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, "src")
sys.path.insert(0, "bench")
from windcheck import tifxyz                                    # noqa: E402
from windcheck.check import load_pairs                          # noqa: E402
from windcheck.intrinsic import SurfaceGraph, oriented_events   # noqa: E402
from windcheck.repair import (apply_field, certified_repair,    # noqa: E402
                              displacement_stats, event_clearance,
                              field_weights, kernel_profile_reset,
                              kernel_profile_snapshot, match_events,
                              merge_events, quantize32,
                              remaining_budget_interval, search_repair)
from crossing_census import census_one                          # noqa: E402
from repair_segment import (CORPORA, RES_UM, SEARCH_REL_VX,     # noqa: E402
                            STRICT_POINT_VX, SUPPORT_VX, accepted,
                            clip_rec, contact_multiset, sha)
from min_exit_sweep import ADMISSIBLE_REL_VX, event_bounds      # noqa: E402

OUT = Path("out/repaired/multi")


class SegmentNotCensusable(RuntimeError):
    """Mesh below the census validity threshold: a disposition, not a bug."""
CENSUS = {"cell": 40.0, "exclude": 1, "maxedge": 60.0}
GEODESIC = Path("engines/geodesic")


class Timers:
    def __init__(self):
        self.t = Counter()
        self.n = Counter()

    def add(self, key: str, dt: float):
        self.t[key] += dt
        self.n[key] += 1

    def timed(self, key: str, fn, *a, **kw):
        t0 = time.time()
        out = fn(*a, **kw)
        self.add(key, time.time() - t0)
        return out

    def report(self) -> dict:
        return {k: {"seconds": round(self.t[k], 2), "calls": self.n[k]}
                for k in sorted(self.t)}


def signature(ev: dict, P: np.ndarray):
    """Canonical event identity: region quads (symmetric under swap) PLUS a
    hash of their corner coordinates. A repair elsewhere that changes this
    event's regions OR their geometry changes the signature, so the event
    is re-evaluated instead of inheriting a stale label (round-22
    amendment 2: recompute after every accepted repair)."""
    a = tuple(sorted(map(tuple, ev["region_a"])))
    b = tuple(sorted(map(tuple, ev["region_b"])))
    regs = tuple(sorted((a, b)))
    h = hashlib.sha256()
    nv, nu = P.shape[:2]
    for reg in regs:
        for v, u in reg:
            for c in ((v, u), (v + 1, u), (v, u + 1), (v + 1, u + 1)):
                if 0 <= c[0] < nv and 0 <= c[1] < nu:
                    h.update(np.asarray(P[c], np.float32).tobytes())
    return (regs, h.hexdigest()[:16])


def build_graph(P, V, diag, timers: Timers):
    g = timers.timed("graph_build", SurfaceGraph,
                     np.asarray(P, np.float64), V, diag)
    gv, gu = np.nonzero(g.idx >= 0)
    X = np.empty((g.n, 3))
    X[g.idx[gv, gu]] = np.asarray(P, np.float64)[gv, gu]
    return g, X


def census_state(mesh: Path, tag: str, wdir: Path, timers: Timers) -> dict:
    """Full census of one mesh: counts, triangle-identity multisets, raw
    records, oriented events per diagonal."""
    counts = timers.timed("oracle_census", census_one, mesh, tag, 1,
                          CENSUS["cell"], 0, CENSUS["maxedge"], wdir)
    if counts is None:
        # census_one declines meshes below its validity threshold (a
        # degenerate upstream duplicate like *_copy has too few valid
        # cells). That is a DISPOSITION, not a crash.
        raise SegmentNotCensusable(
            f"{mesh}: below the census validity threshold")
    st = {"counts": counts, "ms": {}, "rec": {}, "evs": {}, "tag": tag}
    s = tifxyz.read(mesh)
    nv, nu = s.shape
    for d in (0, 1):
        csv = wdir / f"{tag}_d{d}.csv"
        st["ms"][d] = contact_multiset(csv)
        rec = clip_rec(load_pairs(csv), nv, nu)
        st["rec"][d] = rec
        st["evs"][d] = oriented_events(rec) if len(rec) else []
    st["surf"] = s
    return st


def plan_units(st: dict) -> list[dict]:
    """Match d0/d1 events into work units; compute per-constituent fresh
    min-exit bounds; order isolated-first then lowest-bound-first."""
    evs0 = [e for e in st["evs"][0] if not e["ambiguous"]]
    evs1 = [e for e in st["evs"][1] if not e["ambiguous"]]
    amb = [(d, e) for d in (0, 1) for e in st["evs"][d] if e["ambiguous"]]
    m = match_events(evs0, evs1)
    units = []
    for i, j, score, _sw in m["matched"]:
        try:
            ev = merge_events(evs0[i], evs1[j])
        except ValueError:
            units.append({"kind": "orientation_ambiguous", "ev": None,
                          "constituents": [(0, i, evs0[i]), (1, j, evs1[j])],
                          "match_score": score})
            continue
        units.append({"kind": "matched", "ev": ev, "match_score": score,
                      "constituents": [(0, i, evs0[i]), (1, j, evs1[j])]})
    # the reciprocal matcher owns the full partition (round-23 blocker 2):
    # unmatched lists already exclude every ambiguous-group member
    for i in m["unmatched_d0"]:
        units.append({"kind": "d0_only",
                      "ev": {"region_a": set(evs0[i]["region_a"]),
                             "region_b": set(evs0[i]["region_b"])},
                      "constituents": [(0, i, evs0[i])]})
    for j in m["unmatched_d1"]:
        units.append({"kind": "d1_only",
                      "ev": {"region_a": set(evs1[j]["region_a"]),
                             "region_b": set(evs1[j]["region_b"])},
                      "constituents": [(1, j, evs1[j])]})
    for g0, g1 in m["ambiguous"]:
        units.append({"kind": "match_ambiguous", "ev": None,
                      "constituents": [(0, i, evs0[i]) for i in g0]
                      + [(1, j, evs1[j]) for j in g1]})
    for d, e in amb:
        units.append({"kind": "parity_ambiguous", "ev": None,
                      "constituents": [(d, None, e)]})
    P = st["surf"].points
    for u in units:
        u["signature"] = signature(
            u["ev"] if u["ev"] is not None else u["constituents"][0][2], P)
        u["self_touching"] = any(c[2].get("self_touching")
                                 for c in u["constituents"])
    # isolated-first: a unit whose 1-ring-inflated quad set meets no other
    # unit's quads goes ahead of entangled ones (round 18)
    footprints = []
    for u in units:
        quads = set()
        for _, _, e in u["constituents"]:
            quads |= set(map(tuple, e["region_a"]))
            quads |= set(map(tuple, e["region_b"]))
        infl = {(v + dv, u_ + du) for v, u_ in quads
                for dv in (-1, 0, 1) for du in (-1, 0, 1)}
        footprints.append((quads, infl))
    for k, u in enumerate(units):
        u["isolated"] = not any(footprints[k][1] & footprints[j][0]
                                for j in range(len(units)) if j != k)
    return units


def unit_bounds(u: dict, graphs: dict, st: dict, timers: Timers) -> None:
    """Fresh certificate-grade min-exit bound per diagonal constituent on
    the CURRENT geometry. Never averaged across diagonals."""
    u["bounds"] = {}
    for d, _, e in u["constituents"]:
        g, X = graphs[d]
        t0 = time.time()
        L, L_safe, wit, npairs = event_bounds(g, X, st["rec"][d], e)
        timers.add("min_exit", time.time() - t0)
        u["bounds"][d] = {"L_vx": round(L, 6), "L_safe_vx": round(L_safe, 6),
                          "witness": wit, "n_tri_pairs": npairs}
    u["max_L_safe"] = max(b["L_safe_vx"] for b in u["bounds"].values())
    u["prefilter_infeasible"] = u["max_L_safe"] > ADMISSIBLE_REL_VX


def emit_mesh(P_arr, surf, src_mesh: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    Q32 = P_arr.astype(np.float32)
    keep = ~surf.valid
    for i, ax in enumerate(("x", "y", "z")):
        a = Q32[..., i].copy()
        a[keep] = np.asarray(surf.points, np.float32)[keep, i]
        tifffile.imwrite(dst / f"{ax}.tif", a)
    for extra in ("meta.json", "mask.tif", "mask.png"):
        srcf = src_mesh / extra
        if srcf.exists():
            shutil.copy(srcf, dst / extra)


def quality_metrics(P0, P32, inc) -> dict:
    """Continuous warning metrics over affected quads (round 18: disclose;
    W2: recorded, no invented pass thresholds)."""
    qual = {}
    if not inc:
        return {"affected_quads": 0}
    for diag_ in (0, 1):
        areas0, areas1, rots = [], [], []
        for v, u in inc:
            combos = (((v, u), (v, u + 1), (v + 1, u + 1)),
                      ((v, u), (v + 1, u + 1), (v + 1, u))) if diag_ == 0 \
                else (((v, u), (v, u + 1), (v + 1, u)),
                      ((v, u + 1), (v + 1, u + 1), (v + 1, u)))
            for c in combos:
                n0 = np.cross(P0[c[1]] - P0[c[0]], P0[c[2]] - P0[c[0]])
                n1 = np.cross(P32[c[1]] - P32[c[0]], P32[c[2]] - P32[c[0]])
                areas0.append(np.linalg.norm(n0) / 2)
                areas1.append(np.linalg.norm(n1) / 2)
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
    qual["twist_max_before_after_vx"] = [float(np.max(tw0)),
                                         float(np.max(tw1))]
    qual["affected_quads"] = len(inc)
    return qual


def swept_incident_quads(P0, P32, Q, nv, nu) -> set:
    delta = np.linalg.norm(P32 - P0, axis=-1)
    mv, mu = np.nonzero(delta > 0)
    inc = {(v + dv, u + du) for v, u in zip(mv.tolist(), mu.tolist())
           for dv in (-1, 0) for du in (-1, 0)}
    return {(v, u) for v, u in inc
            if 0 <= v < nv - 1 and 0 <= u < nu - 1 and Q[v, u]}


def transact(base: dict, unit: dict, P_orig, wdir: Path, txn_id: int,
             timers: Timers, log) -> dict:
    """One transactional repair attempt against the current base mesh.

    Returns {"outcome": "accepted", "mesh": path, ...} or
    {"outcome": "failed", "reason": ...}. Base mesh is untouched on
    failure; the emitted mesh lives in its own txn directory.

    Wraps the actual attempt in a per-transaction kernel-profile scope
    (reviewer spec A4): counters reset at entry, snapshot stored under
    "kernel_profile" in every result, so the certificate keeps it."""
    kernel_profile_reset()
    res = _transact(base, unit, P_orig, wdir, txn_id, timers, log)
    res["kernel_profile"] = kernel_profile_snapshot()
    if not res.get("budget_skipped_candidates"):
        res.pop("budget_skipped_candidates", None)
    return res


def _transact(base: dict, unit: dict, P_orig, wdir: Path, txn_id: int,
              timers: Timers, log) -> dict:
    surf = base["st"]["surf"]
    P, V = surf.points, surf.valid
    nv, nu = surf.shape
    ev = unit["ev"]
    # propose on the diagonal where the crossing exists: a d1-only event
    # has no intersecting pairs under the d0 triangulation, so proposing
    # on d0 would return no candidates at all (certification is always
    # both-diagonal regardless of the proposal graph)
    diag_pick = 1 if unit["kind"] == "d1_only" else 0
    g, _X = base["graphs"][diag_pick]
    before_ms = base["st"]["ms"]
    result = {"txn": txn_id, "signature": str(unit["signature"]),
              "kind": unit["kind"], "attempted_candidates": []}
    # (direction, lam) combos excluded by an empty remaining-budget
    # interval before any probe (round-24 Part C: skipped AND labelled);
    # stripped by transact() when empty
    skipped: list = []
    result["budget_skipped_candidates"] = skipped

    # ---- rigid proposal (ranking only) + primary certified candidate
    t0 = time.time()
    try:
        r = event_clearance(g, _X, ev, base["voxel_um"],
                            t_max_vx=SEARCH_REL_VX)
    except ValueError as e:
        timers.add("search", time.time() - t0)
        return {**result, "outcome": "failed",
                "reason": f"proposal error: {e}"}
    ranked, rep, P32 = None, None, None
    try:
        if r is not None:
            evd = ev if r["side"] == "a" else \
                {"region_a": ev["region_b"], "region_b": ev["region_a"]}
            out = certified_repair(g, P, V, evd, r["direction"], r["t_vx"],
                                   mode="symmetric",
                                   budget_vx=SEARCH_REL_VX,
                                   support_vx=SUPPORT_VX,
                                   P_orig=P_orig,
                                   budget_point_vx=STRICT_POINT_VX)
            if out is not None:
                P32, rep = out
                ranked = [{"direction": r["direction"], "lam": 0.5,
                           "t_rel": rep["applied_relative_vx"], "P32": P32,
                           "field_report": rep,
                           "provenance": "primary_certified"}]
                ev = evd
        if ranked is None:
            log("    no primary candidate; ranked search")
            ranked = search_repair(g, P, V, ev,
                                   budget_point_vx=STRICT_POINT_VX,
                                   support_vx=SUPPORT_VX,
                                   P_orig=P_orig,
                                   rel_cap_vx=SEARCH_REL_VX,
                                   skipped_log=skipped)
            if not ranked:
                timers.add("search", time.time() - t0)
                return {**result, "outcome": "failed",
                        "reason": "no locally-clean candidate in the "
                                  "searched direction x split family"}
            c = ranked[0]
            c["provenance"] = "search"
            P32 = c["P32"]
            rep = displacement_stats(P, P32, V, CENSUS["maxedge"])
            rep.update(c["field_report"])
            rep.update({"applied_relative_vx": c["t_rel"],
                        "candidate_lp_exit_vx": None})
    except ValueError as e:
        timers.add("search", time.time() - t0)
        return {**result, "outcome": "failed", "reason": str(e)}
    timers.add("search", time.time() - t0)
    result["n_candidates"] = len(ranked)

    active = ranked[0]
    t_rel = rep["applied_relative_vx"]
    t_dirty = 0.0
    dst = wdir / f"txn_{txn_id}.tifxyz"
    tag = f"{base['tag']}_t{txn_id}"
    # deformation graph and harmonic fields are invariant across all
    # probes of this transaction (same base mesh, same event); ws is the
    # apply_field patch-not-copy buffer (round 24), valid for this base P
    gf = SurfaceGraph(np.asarray(P, np.float64), V, diagonal=-1,
                      maxedge=g.maxedge)
    fcache: dict = {}
    ws: dict = {}

    def probe(P_arr):
        timers.timed("emit", emit_mesh, P_arr, surf, base["mesh"], dst)
        timers.timed("oracle_census", census_one, dst, tag, 1,
                     CENSUS["cell"], 0, CENSUS["maxedge"], wdir)
        ams = {d: contact_multiset(wdir / f"{tag}_d{d}.csv") for d in (0, 1)}
        return ams, all(accepted(before_ms[d], ams[d], ev)[0]
                        for d in (0, 1))

    # ---- remaining-budget intervals (round-24 Part C) -------------------
    # exact per-vertex admissible t interval vs the ORIGINAL mesh, per
    # (direction, lambda) -- the weights change with the split. Applied to
    # candidate construction (inside search/certified), to GROWTH below,
    # and (by the t_rel <= cap invariant) to every bracket probe.
    _iv_cache: dict = {}

    def active_interval(cand):
        lam = cand.get("lam", 0.5)
        key = (tuple(float(x) for x in cand["direction"]),
               lam if isinstance(lam, str) else float(lam))
        if key not in _iv_cache:
            w = field_weights(gf, ev, lam, SUPPORT_VX, fcache)
            _iv_cache[key] = remaining_budget_interval(
                w, gf, P, P_orig, cand["direction"], STRICT_POINT_VX, lam)
        return _iv_cache[key]

    # ---- oracle-failure state machine (round-24 Part A ruling) ----------
    # EXPLICIT phases replace the accidental break-driven flow:
    #   HEAD -- consume unique legacy-head direction switches in legacy
    #           order: the initial ranked list, then (once it exhausts to
    #           one) the one-shot doubled-budget re-search's non-parallel
    #           head candidates. Selection and order match the legacy scan
    #           (search_repair's head/tail guarantee).
    #   GROW -- the legacy grow-t path for the ACTIVE candidate, bounded
    #           by the relative cap and the remaining-budget interval.
    #   TAIL -- ONLY when HEAD+GROW reach their budget/attempt failure,
    #           the re-search's "extra": True candidates are tried, each
    #           attempt marked rescue_extension=true. Previously ACCEPTED
    #           goldens are identical by construction: TAIL engages only
    #           after the legacy path would have FAILED.
    # A seen-set over (direction, lambda, t, base-mesh-hash) means no
    # identical probe ever runs twice ((iv) of the ruling): probes are
    # deterministic, so a duplicate could only re-fail; it is recorded as
    # duplicate_probe_skipped and consumes no attempt. Exhaustion falls
    # through to the existing "engine oracle never clean" failure.
    MAX_ATTEMPTS = 12        # legacy probe budget (HEAD + GROW)
    MAX_TAIL_ATTEMPTS = 12   # additional rescue-extension probe budget
    base_id = base.get("x_sha", "")

    def probe_key(cand_dir, lam, t):
        return (tuple(float(x) for x in cand_dir),
                lam if isinstance(lam, str) else float(lam),
                float(t), base_id)

    def cand_record(cand, t, outcome, rescue_flag=False):
        rec = {"direction": list(map(float, cand["direction"])),
               "lam": cand.get("lam", 0.5), "t_rel": t, "outcome": outcome}
        if rescue_flag:
            rec["rescue_extension"] = True
        return rec

    after_ms, ok = probe(P32)
    seen = {probe_key(active["direction"], active.get("lam", 0.5), t_rel)}
    fallbacks = 0
    attempts = 0             # HEAD + GROW probes after the initial one
    tail_attempts = 0
    fb_pool = None           # the one-shot re-search pool (round 24)
    tail_pool: list = []     # its "extra": True candidates, rescue-only
    legacy_failure = None    # the failure the legacy path reached
    rescue = False           # active candidate is a rescue extension
    phase = "HEAD"
    while not ok:
        if phase == "HEAD":
            if attempts >= MAX_ATTEMPTS:
                legacy_failure = ("engine oracle never clean within "
                                  f"{MAX_ATTEMPTS} attempts")
                phase = "TAIL"
                continue
            if len(ranked) > 1:
                nxt = ranked[1]
                nxt.setdefault("provenance", "search")
                k = probe_key(nxt["direction"], nxt.get("lam", 0.5),
                              nxt["t_rel"])
                if k in seen:            # never retry an identical probe
                    ranked.pop(1)
                    result["attempted_candidates"].append(
                        cand_record(nxt, nxt["t_rel"],
                                    "duplicate_probe_skipped"))
                    continue
                result["attempted_candidates"].append(
                    cand_record(active, t_rel, "oracle_failed"))
                ranked.pop(0)
                active = ranked[0]
                P32 = active["P32"]
                t_rel = active["t_rel"]
                t_dirty = 0.0
                fallbacks += 1
                seen.add(k)
                log(f"    fallback -> lam {active['lam']} t {t_rel:.3f}")
                after_ms, ok = probe(P32)
                attempts += 1
            elif fb_pool is None:
                t0 = time.time()
                fb_pool = search_repair(g, P, V, ev,
                                        budget_point_vx=STRICT_POINT_VX,
                                        support_vx=SUPPORT_VX,
                                        extra_candidates=6,
                                        P_orig=P_orig,
                                        rel_cap_vx=SEARCH_REL_VX,
                                        skipped_log=skipped)
                timers.add("search", time.time() - t0)
                d0v = np.asarray(active["direction"])
                ranked += [c for c in fb_pool if not c.get("extra")
                           and np.dot(c["direction"], d0v) < 0.999]
                tail_pool = [c for c in fb_pool if c.get("extra")]
            else:
                phase = "GROW"           # head switches exhausted
            continue
        if phase == "GROW":
            if t_rel > SEARCH_REL_VX:
                legacy_failure = ("relative search budget exhausted before "
                                  "the engine oracle was clean")
                phase = "TAIL"
                continue
            if attempts >= MAX_ATTEMPTS:
                legacy_failure = ("engine oracle never clean within "
                                  f"{MAX_ATTEMPTS} attempts")
                phase = "TAIL"
                continue
            iv = active_interval(active)
            t_next = t_rel + 0.08
            if iv is None or t_next > iv[1]:
                legacy_failure = ("remaining per-vertex budget interval "
                                  "exhausted before the engine oracle was "
                                  "clean")
                phase = "TAIL"
                continue
            t_dirty, t_rel = t_rel, t_next
            P2, _ = apply_field(g, P, V, ev, active["direction"], t_rel,
                                active.get("lam", 0.5), SUPPORT_VX,
                                gf=gf, field_cache=fcache, workspace=ws)
            P32 = quantize32(P2)
            seen.add(probe_key(active["direction"], active.get("lam", 0.5),
                               t_rel))
            after_ms, ok = probe(P32)
            attempts += 1
            continue
        # phase == "TAIL": rescue extensions, only after legacy failure
        nxt = None
        while tail_pool:
            c = tail_pool.pop(0)
            k = probe_key(c["direction"], c.get("lam", 0.5), c["t_rel"])
            if k in seen:                # never retry an identical probe
                result["attempted_candidates"].append(
                    cand_record(c, c["t_rel"], "duplicate_probe_skipped",
                                rescue_flag=True))
                continue
            nxt = (c, k)
            break
        if nxt is None or tail_attempts >= MAX_TAIL_ATTEMPTS:
            break                        # exhausted -> failure below
        result["attempted_candidates"].append(
            cand_record(active, t_rel, "oracle_failed", rescue_flag=rescue))
        active, k = nxt
        active.setdefault("provenance", "search")
        rescue = True
        P32 = active["P32"]
        t_rel = active["t_rel"]
        t_dirty = 0.0
        fallbacks += 1
        tail_attempts += 1
        seen.add(k)
        log(f"    rescue tail -> lam {active['lam']} t {t_rel:.3f} "
            f"(after: {legacy_failure})")
        after_ms, ok = probe(P32)
    if not ok:
        result["attempted_candidates"].append(
            cand_record(active, t_rel, "oracle_failed", rescue_flag=rescue))
        reason = (legacy_failure if tail_attempts == 0 else
                  f"engine oracle never clean within {MAX_ATTEMPTS} attempts")
        out = {**result, "outcome": "failed", "reason": reason}
        if tail_attempts:
            out["rescue_extension_attempts"] = tail_attempts
            out["legacy_failure"] = legacy_failure
        return out
    result["fallbacks"] = fallbacks
    if tail_attempts:
        result["rescue_extension_attempts"] = tail_attempts
        result["legacy_failure_before_rescue"] = legacy_failure

    # ---- engine-driven bracket refinement (round 18)
    # Part C cap on bracket probes: every mid lies in [t_dirty, t_rel] and
    # t_rel was admitted by the winner's remaining-budget interval (search/
    # certified construction, GROW bound), so no bracket probe can exceed
    # the cap; assert the invariant rather than trust it.
    _iv = active_interval(active)
    assert _iv is not None and t_rel <= _iv[1] * (1 + 1e-9) + 1e-12, \
        "bracket start exceeds the remaining-budget cap"
    RESOLUTION = 0.01
    lo, hi = t_dirty, t_rel

    def clean_at(t):
        P2, _ = apply_field(g, P, V, ev, active["direction"], t,
                            active.get("lam", 0.5), SUPPORT_VX,
                            gf=gf, field_cache=fcache, workspace=ws)
        Pq = quantize32(P2)
        ams, okq = probe(Pq)
        return Pq, ams, okq

    while hi - lo > RESOLUTION:
        mid = 0.5 * (lo + hi)
        _, _, okm = clean_at(mid)
        if okm:
            hi = mid
        else:
            lo = mid
    t_rel = hi
    P32, after_ms, okf = clean_at(t_rel)
    if not okf:
        return {**result, "outcome": "failed",
                "reason": "bracket endpoint not clean on re-probe"}

    # ---- hard gates from the FINAL emitted mesh (pre-registered)
    s2 = tifxyz.read(dst)
    if s2.shape != surf.shape or not (s2.valid == V).all():
        return {**result, "outcome": "failed",
                "reason": "validity mask changed"}
    rep2 = displacement_stats(P, P32, V, CENSUS["maxedge"])
    gates = {"quads_newly_dropped": rep2["quads_newly_dropped"],
             "quads_newly_retained": rep2["quads_newly_retained"],
             "triangle_inversions_d0": rep2["triangle_inversions_d0"],
             "triangle_inversions_d1": rep2["triangle_inversions_d1"]}
    for d in (0, 1):
        g2 = SurfaceGraph(np.asarray(s2.points, np.float64), V, d)
        gates[f"component_delta_d{d}"] = g2.ncomp - base["graphs"][d][0].ncomp
    cum = np.linalg.norm(np.asarray(s2.points, np.float64)
                         - P_orig, axis=-1)
    gates["max_point_disp_vx"] = rep2["max_disp_vx"]
    gates["cumulative_max_disp_vx"] = float(cum.max())
    hard_fail = ([k for k, v in gates.items()
                  if k.startswith(("quads_", "triangle_", "component_"))
                  and v != 0]
                 + (["max_point_disp_vx"]
                    if rep2["max_disp_vx"] > STRICT_POINT_VX else [])
                 + (["cumulative_max_disp_vx"]
                    if gates["cumulative_max_disp_vx"] > STRICT_POINT_VX
                    else []))
    if hard_fail:
        return {**result, "outcome": "failed", "gates": gates,
                "reason": f"hard gate failed: {hard_fail}"}

    from windcheck.intrinsic import retained_quads
    P0 = np.asarray(P, np.float64)
    inc = swept_incident_quads(P0, P32,
                               retained_quads(P0, V, CENSUS["maxedge"]),
                               nv, nu)
    # Round-23 blocker 3: the repair block must describe the WINNING
    # candidate, never the failed primary (the round-21 bug repeated).
    # Re-apply the winner's field once at the certified t to get the
    # authoritative report; rep2 already holds the final-mesh stats.
    _, win_report = apply_field(g, P, V, ev, active["direction"], t_rel,
                                active.get("lam", 0.5), SUPPORT_VX,
                                gf=gf, field_cache=fcache, workspace=ws)
    rep = dict(win_report)
    rep.update({"applied_relative_vx": t_rel,
                "candidate_lp_exit_vx":
                    active.get("t_rel")
                    if active.get("provenance") == "primary_certified"
                    else None})
    rep.update({f"final_{k}": v for k, v in rep2.items()})
    result.update({
        "outcome": "accepted", "mesh": str(dst),
        "engine_certified_relative_vx": t_rel,
        "oracle_bracket_vx": [lo, hi],
        "oracle_resolution_vx": RESOLUTION,
        "winning_candidate": {
            "direction": list(map(float, active["direction"])),
            "lam": active.get("lam", 0.5),
            "provenance": active.get("provenance", "?"),
            **({"rescue_extension": True} if rescue else {})},
        "repair": {k: v for k, v in rep.items() if k != "P32"},
        "gates": gates,
        "quality_warning_metrics": quality_metrics(P0, P32, inc),
        "after_ms": after_ms,
        "census_csvs": {f"d{d}": str(wdir / f"{tag}_d{d}.csv")
                        for d in (0, 1)},
    })
    return result


def resolve_segment(sub: str):
    exact, subs = [], []
    for corpus, root, volume, work in CORPORA:
        for d in Path(root).iterdir() if Path(root).exists() else []:
            if d.is_dir() and sub in d.name:
                (exact if d.name == sub else subs).append(
                    (corpus, d, volume))
    cands = exact or subs
    assert len(cands) == 1, \
        f"{sub}: {len(cands)} matches {[c[1].name for c in cands]}"
    return cands[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--segment", required=True)
    ap.add_argument("--plan-only", action="store_true",
                    help="census + match + prefilter, no repairs")
    ap.add_argument("--max-rounds", type=int, default=12)
    args = ap.parse_args()
    corpus, seg, VOLUME = resolve_segment(args.segment)
    mesh = sorted(seg.glob(f"mesh/*{VOLUME}*.tifxyz"))[0]
    m = RES_UM.search(mesh.name)
    voxel_um = float(m.group(1)) if m else 7.91
    tag = hashlib.sha256(seg.name.encode()).hexdigest()[:12]
    # probe/census scratch may be redirected to fast storage (tmpfs) on
    # highly concurrent hosts; certificates and archives stay under OUT
    import os
    work_root = Path(os.environ.get("WINDCHECK_WORK_ROOT", str(OUT)))
    wdir = work_root / f"work_{tag}"
    wdir.mkdir(parents=True, exist_ok=True)
    timers = Timers()
    t_start = time.time()

    def log(msg):
        print(msg, flush=True)

    log(f"{corpus} {seg.name}  voxel {voxel_um} um  workdir {wdir}")
    try:
        orig = census_state(mesh, f"{tag}_orig", wdir, timers)
    except SegmentNotCensusable as e:
        cert = {"claim": ("Segment below the census validity threshold; "
                          "no audit or repair is defined on it."),
                "segment_class": "not_censusable", "corpus": corpus,
                "segment": seg.name, "volume": VOLUME, "reason": str(e),
                "hashes": {"input": {f"{ax}.tif": sha(mesh / f"{ax}.tif")
                                     for ax in ("x", "y", "z")}},
                "code_commit": subprocess.run(
                    ["git", "rev-parse", "HEAD"], capture_output=True,
                    text=True).stdout.strip()}
        cpath = OUT / f"{seg.name}_multi_certificate.json"
        cpath.write_text(json.dumps(cert, indent=1))
        log(f"\n[NOT_CENSUSABLE] certificate {cpath}")
        return
    P_orig = np.asarray(orig["surf"].points, np.float64)
    n_orig = {d: len([e for e in orig["evs"][d] if not e["ambiguous"]])
              for d in (0, 1)}
    log(f"  original events: d0 {n_orig[0]}  d1 {n_orig[1]}")

    labelled: dict = {}          # signature -> disposition record
    label_history: list = []     # superseded labels, keyed by base mesh
    txn_log: list = []
    txn_walls: list = []         # wall seconds per transaction attempt
    cur_mesh, cur_st = mesh, orig
    txn_id = 0

    for rnd in range(args.max_rounds):
        base_hash = sha(cur_mesh / "x.tif")
        graphs = {d: build_graph(cur_st["surf"].points,
                                 cur_st["surf"].valid, d, timers)
                  for d in (0, 1)}
        base = {"st": cur_st, "graphs": graphs, "mesh": cur_mesh,
                "voxel_um": voxel_um, "tag": tag,
                "x_sha": base_hash}   # seen-set probe identity (round 25)
        units = plan_units(cur_st)
        if not units:
            log(f"  round {rnd}: no events remain")
            break
        for u in units:
            if u["ev"] is not None:
                unit_bounds(u, graphs, cur_st, timers)
        attemptable = [u for u in units
                       if u["ev"] is not None
                       and u["signature"] not in labelled
                       and not u["prefilter_infeasible"]]
        attemptable.sort(key=lambda u: (not u["isolated"], u["max_L_safe"]))
        for u in units:
            if u["signature"] in labelled:
                continue
            if u["ev"] is None:
                labelled[u["signature"]] = {
                    "disposition": "inconclusive", "kind": u["kind"],
                    "reason": "ambiguous match/orientation/parity -- "
                              "reported, never guessed"}
                log(f"  labelled {u['kind']}: {u['signature'][0][:2]}...")
            elif u["prefilter_infeasible"]:
                labelled[u["signature"]] = {
                    "disposition": "certified_infeasible",
                    "kind": u["kind"], "bounds": u["bounds"],
                    "model_scope": "rules out rigid relative core "
                                   "translation <= admissible only; not "
                                   "non-rigid deformation, remeshing, or "
                                   "cutting",
                    "admissible_rel_vx": ADMISSIBLE_REL_VX,
                    "provisional": True}
                log(f"  prefilter: certified-infeasible "
                    f"(max L_safe {u['max_L_safe']:.3f} vx) "
                    f"{sorted(u['ev']['region_a'])[:2]}...")
        if args.plan_only:
            log("  plan-only: units " + json.dumps([
                {"kind": u["kind"], "isolated": u["isolated"],
                 "self_touching": u["self_touching"],
                 "bounds": {f"d{d}": b["L_safe_vx"]
                            for d, b in u.get("bounds", {}).items()},
                 "prefilter_infeasible": u.get("prefilter_infeasible")}
                for u in units], indent=1))
            log("  plan-only: instrumentation "
                + json.dumps(timers.report()))
            return
        accepted_this_round = False
        for u in attemptable:
            log(f"  attempt {u['kind']} max_L_safe {u['max_L_safe']:.3f} "
                f"isolated {u['isolated']} "
                f"A{sorted(u['ev']['region_a'])[:2]} "
                f"B{sorted(u['ev']['region_b'])[:2]}")
            t0_txn = time.time()
            res = transact(base, u, P_orig, wdir, txn_id, timers, log)
            txn_walls.append(time.time() - t0_txn)
            txn_id += 1
            if res["outcome"] == "accepted":
                log(f"    ACCEPTED at {res['engine_certified_relative_vx']:.4f} "
                    f"vx relative, max point "
                    f"{res['gates']['max_point_disp_vx']:.4f} vx")
                txn_log.append({k: v for k, v in res.items()
                                if k != "after_ms"})
                cur_mesh = Path(res["mesh"])
                cur_st = census_state(cur_mesh, f"{tag}_r{rnd+1}", wdir,
                                      timers)
                # Round-23 blocker 1: an accepted repair can change
                # obstacle geometry, harmonic support, or remaining
                # budget far from an event's own corners. EVERY
                # provisional label is invalidated and re-evaluated
                # against the new base; old attempts kept as history.
                for kk, vv in labelled.items():
                    label_history.append({"base_mesh_x_sha256": base_hash,
                                          "signature": str(kk), **vv})
                labelled.clear()
                accepted_this_round = True
                break                      # full recensus + rematch (round 18)
            else:
                # cheap lower bound on every failure (disposition machine)
                bounds = u.get("bounds", {})
                worst = max((b["L_safe_vx"] for b in bounds.values()),
                            default=0.0)
                disp = ("certified_infeasible"
                        if worst > ADMISSIBLE_REL_VX else "inconclusive")
                labelled[u["signature"]] = {
                    "disposition": disp, "kind": u["kind"],
                    "reason": res["reason"], "bounds": bounds,
                    "attempt": {k: v for k, v in res.items()
                                if k not in ("after_ms",)},
                    "label": "no displacement repair found under the "
                             "stated finite search and budget"}
                log(f"    FAILED: {res['reason']}")
        if not accepted_this_round:
            break

    # ---- final state: labels come only from the final accepted mesh.
    # Fresh census so c_atlas.bin and the CSVs both describe THIS mesh
    # (probe censuses of rejected candidates may have overwritten them).
    # THE DELIVERABLE: the final accepted mesh must live in the durable
    # output tree, never only in the (possibly tmpfs) workdir -- the
    # certificate alone cannot reconstruct coordinates
    mesh_out = None
    if cur_mesh != mesh:
        mesh_out = OUT / "meshes" / f"{seg.name}_repaired.tifxyz"
        if mesh_out.exists():
            shutil.rmtree(mesh_out)
        mesh_out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(cur_mesh, mesh_out)
        cur_mesh = mesh_out
    final = census_state(cur_mesh, f"{tag}_final", wdir, timers)
    graphs = {d: build_graph(final["surf"].points, final["surf"].valid, d,
                             timers) for d in (0, 1)}
    residual_units = plan_units(final)
    for u in residual_units:
        if u["ev"] is not None:
            unit_bounds(u, graphs, final, timers)
            u["final_disposition"] = (
                "certified_infeasible"
                if u["max_L_safe"] > ADMISSIBLE_REL_VX else "inconclusive")
        else:
            u["final_disposition"] = "inconclusive_ambiguous"
        prior = labelled.get(u["signature"])
        u["prior_label"] = prior["disposition"] if prior else None
        u["prior_reason"] = prior.get("reason") if prior else None

    # residual spectrum at intrinsic scale via the promoted engine
    # (census_state just ran on the final mesh, so c_atlas.bin matches)
    spectrum = {}
    if GEODESIC.exists():
        for d in (0, 1):
            outj = wdir / f"final_spectrum_d{d}.json"
            csv = wdir / f"{final['tag']}_d{d}.csv"
            if csv.exists():
                t0 = time.time()
                subprocess.run([str(GEODESIC), str(wdir / "c_atlas.bin"),
                                str(csv), str(outj), "0", str(d),
                                str(CENSUS["maxedge"]), str(voxel_um)],
                               check=True, capture_output=True)
                timers.add("geodesic_rematch", time.time() - t0)
                spectrum[f"d{d}"] = json.loads(outj.read_text())

    none_new = all(
        n <= orig["ms"][d].get(k, 0)
        for d in (0, 1) for k, n in final["ms"][d].items())
    n_final = {d: len(final["evs"][d]) for d in (0, 1)}
    fully_clean = all(len(final["ms"][d]) == 0 for d in (0, 1))
    cum = np.linalg.norm(np.asarray(final["surf"].points, np.float64)
                         - P_orig, axis=-1)
    # round-23: orthogonal status fields; "partial" is reserved for
    # residual TRANSVERSE events, never for leftover grazing alone
    n_trans = {d: sum(n for k, n in final["ms"][d].items()
                      if k[-1] == "transverse") for d in (0, 1)}
    n_graze = {d: sum(n for k, n in final["ms"][d].items()
                      if k[-1] == "grazing") for d in (0, 1)}
    transverse_clean = sum(n_trans.values()) == 0
    seg_class = ("clean" if fully_clean else
                 "transverse_clean" if transverse_clean else
                 "partial" if txn_log and none_new else
                 "unresolved")
    archive = OUT / "archive" / seg.name
    archive.mkdir(parents=True, exist_ok=True)
    csv_records = {}
    for w, t_ in (("orig", orig["tag"]), ("final", final["tag"])):
        for d in (0, 1):
            src = wdir / f"{t_}_d{d}.csv"
            if src.exists():
                dstc = archive / f"{w}_d{d}.csv"
                shutil.copy(src, dstc)
                csv_records[f"{w}_d{d}"] = {"path": str(dstc),
                                            "sha256": sha(dstc)}
    porcelain = subprocess.run(["git", "status", "--porcelain"],
                               capture_output=True, text=True).stdout
    cert = {
        "claim": (
            f"Segment {seg.name}: {len(txn_log)} crossing event(s) removed "
            f"transactionally; residual events d0 {n_final[0]} / d1 "
            f"{n_final[1]}, each dispositioned from the final mesh. "
            + ("The final mesh is globally recensus-clean for non-adjacent "
               "transverse, coplanar and grazing contacts under both "
               "triangulations." if fully_clean else
               "No contact multiplicity increased at any point "
               "(residuals enumerated below).")
            + " Endpoint certificates only; texture fidelity and "
              "collision-free deformation through time are not certified."),
        "segment_class": seg_class,
        "transverse_status": "clean" if transverse_clean else "residual",
        "contact_status": ("clean" if fully_clean else
                           "grazing_residual"
                           if sum(n_graze.values()) and transverse_clean
                           else "contact_residual"),
        "residual_transverse": n_trans,
        "residual_grazing": n_graze,
        "census_csv_archive": csv_records,
        "label_history": label_history,
        "corpus": corpus, "segment": seg.name, "volume": VOLUME,
        "voxel_um": voxel_um,
        "fidelity_tier": "strict_1vx_max_emitted_vertex_displacement",
        "tier_binds": "max point displacement after float32 emission, "
                      "per transaction AND cumulative vs the original mesh",
        "search_relative_cap_vx": SEARCH_REL_VX,
        "support_vx": SUPPORT_VX,
        "admissible_rel_vx": ADMISSIBLE_REL_VX,
        "original_events": n_orig,
        "final_events": n_final,
        "none_new": none_new,
        "cumulative_max_disp_vx": float(cum.max()),
        "transactions": txn_log,
        "event_dispositions": [
            {"signature": str(u["signature"]),
             "kind": u["kind"],
             "final_disposition": u["final_disposition"],
             "bounds": {f"d{d}": b for d, b in u.get("bounds", {}).items()},
             "prior_label": u["prior_label"],
             "prior_reason": u["prior_reason"]}
            for u in residual_units],
        "labelled_during_run": {str(k): v for k, v in labelled.items()},
        "census_before": {f"d{d}": orig["counts"][f"d{d}"] for d in (0, 1)},
        "census_final": {f"d{d}": final["counts"][f"d{d}"] for d in (0, 1)},
        "residual_spectrum": spectrum,
        "instrumentation": {**timers.report(),
                            "wall_seconds": round(time.time() - t_start, 1),
                            "transactions_attempted": txn_id,
                            "transactions_accepted": len(txn_log),
                            "transaction_wall_s": (
                                {"p50": round(float(np.percentile(
                                     txn_walls, 50)), 3),
                                 "p95": round(float(np.percentile(
                                     txn_walls, 95)), 3),
                                 "max": round(float(max(txn_walls)), 3),
                                 "n": len(txn_walls)}
                                if txn_walls else None)},
        "hashes": {
            "input": {f"{ax}.tif": sha(mesh / f"{ax}.tif")
                      for ax in ("x", "y", "z")},
            "output_mesh_path": str(mesh_out) if mesh_out else None,
            "output": ({f"{ax}.tif": sha(cur_mesh / f"{ax}.tif")
                        for ax in ("x", "y", "z")}
                       if cur_mesh != mesh else None),
            "engine_binary": sha(Path("engines/selfcross")),
            "engine_source": sha(Path("engines/selfcross.cpp")),
            "geodesic_binary": sha(GEODESIC) if GEODESIC.exists() else None,
        },
        "code_commit": subprocess.run(["git", "rev-parse", "HEAD"],
                                      capture_output=True,
                                      text=True).stdout.strip(),
        "git_status_porcelain": porcelain,
        "executor_sha256": sha(Path("bench/repair_multi.py")),
        "uv_lock_sha256": sha(Path("uv.lock")),
        "params": {**CENSUS, "touch_tol": 1e-3},
        "note": ("Acceptance authority is the reloaded selfcross census "
                 "(triangle-identity multisets, round 19). The geodesic "
                 "engine is used only for the residual spectrum. Min-exit "
                 "bounds rule out rigid relative core translation within "
                 "the admissible budget only."),
    }
    cpath = OUT / f"{seg.name}_multi_certificate.json"
    cpath.write_text(json.dumps(cert, indent=1))
    log(f"\n[{seg_class.upper()}] {len(txn_log)} repaired, "
        f"residual d0 {n_final[0]} / d1 {n_final[1]}; certificate {cpath}")


if __name__ == "__main__":
    main()
