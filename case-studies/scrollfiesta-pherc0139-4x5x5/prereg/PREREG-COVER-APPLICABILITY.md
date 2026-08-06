# Pre-registration: alternate-cover applicability (topology-aware export)

Written 2026-08-05, BEFORE the candidate instrumentation ran, per house
rule. Round-41 rank-1 gate: the topology-aware exporter (choose among
Fiesta's OWN candidate covers instead of first-cover-wins) is built
ONLY if alternate covers exist where the crossings are. This experiment
measures existence and nothing else. No selector is built here.

## Instrumentation

Env-gated candidate tap in the strip rasterizer (`tifxyz_export.c`),
same discipline as the proven facemap tap: `SF_TXZ_CANDIDATES=<path>`
streams one binary record per cover attempt at every in-triangle pixel
(pixel index u64; face index i32; candidate world x, y, z as f32).
NEW file only; with the variable unset the code path is a single
getenv. Proof obligation: with candidates ON, every standard output of
the steps-12345 export is byte-identical to the already-recorded
canonical run (the pipeline's determinism is established).

## Definitions (frozen)

- Contested pixel: a valid output pixel with >= 2 recorded candidates.
- Event: a final-export region-pair event (the published clustering,
  verbatim), with its constituent transverse rows.
- Local patch of a pixel: the up-to-4 grid cells having that pixel as
  a corner, plus the 8-neighbouring cells of each (3x3 cell patch per
  changed pixel).
- A candidate switch replaces one pixel's (x,y,z) with another recorded
  candidate's coordinates for that same pixel, unchanged from Fiesta.
  Validity is never changed; no pixel is invalidated or invented.
- Local acceptance of a switch: within the union of local patches of
  all switched pixels, under BOTH diagonals, (a) the event's transverse
  rows strictly decrease, (b) no new transverse row appears, (c) no
  quad newly fails the frozen maxedge gate (60.0).

## Existence procedure (frozen; conservative by construction)

Per event, a deterministic greedy search: iterate the event's
transverse rows in engine emission order; for each of the row's up to 8
corner pixels, in (v, u) order, try every alternate candidate for that
pixel in recorded stream order; accept the first locally-accepted
switch and restart the row scan; stop when the event has zero
transverse rows in its local union (SOLVED) or a full pass yields no
accepted switch (UNSOLVED). Greedy failure does NOT prove
nonexistence, so the measured pass rates UNDERCOUNT true existence;
the thresholds below are therefore conservative.

## Gates (frozen, from the accepted round-41 plan)

Evaluated on the canonical final export, each diagonal's event set
separately, both diagonals required:

- G1: >= 50% of crossing events touch >= 1 contested pixel (a corner
  pixel of a constituent row with >= 2 candidates).
- G2: >= 30% of events are SOLVED by the existence procedure.

BOTH pass -> the selector is designed and pre-registered separately.
EITHER fails -> the exporter line is killed; the fallback (ground-truth
merge-edge benchmark) proceeds; the negative is published.

## Outputs

Candidate-stream statistics (pixels, candidates per pixel
distribution, contested count vs the exporter's own multi/conflict
stats as a consistency check); per-event: contested-touch flag,
SOLVED/UNSOLVED, switches used; the two gate fractions per diagonal;
the verdict. All counts reported, including zeros.

---

## ADDENDUM 1 — 2026-08-06, before any per-event outcome was seen

Recorded while the first run was in flight (150/654 events complete;
only aggregate progress counts observed, no per-event outcome, no gate
fraction). Projected completion 30-50 h, which consumes the gate's own
3-day budget. Three amendments, all scheduling/implementation policy in
the round-28 tradition (the frozen 600 s rule: scheduling, never a
hardness claim):

1. PER-EVENT WALL-CLOCK BUDGET: 15 minutes. An event exceeding it
   reports UNSOLVED_TIMEOUT, counted as UNSOLVED for G2 and reported
   separately. Conservative by construction: timeouts can only lower
   the measured pass rate, never raise it.
2. CLARIFICATION of the row scan: rows of the event that are already
   non-transverse in the current switched state are skipped when
   choosing switch sites (a switch there cannot remove a remaining
   row's transversality via the acceptance rule, which requires strict
   event-row decrease). Trial order is otherwise unchanged.
3. CHECKPOINTING: per-event results are appended to a JSONL as they
   complete; a restart skips completed events. No measured value
   changes.

The first run's partial results are discarded unread; the restarted run
recomputes every event under the amended procedure.

## ADDENDUM 2 — 2026-08-06, after external review, before any outcome read

The Addendum-1 run was stopped and its checkpoint DISCARDED UNREAD on
the external reviewer's finding that the implementation did not enforce the
frozen acceptance rule. Corrections, all committed before any per-event
outcome existed:

1. ACCEPTANCE, corrected and strengthened (reviewer's conservative
   form): a switch is accepted iff, over the trial's local union under
   BOTH diagonals, (a) the event's unresolved target rows strictly
   decrease; (b) NO transverse, grazing or coplanar contact key is new
   or increased in multiplicity — keys are the engine's schema-v2
   triangle-level identities (diag, v1,u1,tri1, v2,u2,tri2, verdict),
   compared as multisets; resurrection of previously resolved event
   rows is NOT forgiven; (c) the retained-quad set (the engine's own
   maxedge criterion) satisfies before ⊆ after.
2. LOCAL UNION, clarified to the frozen text's letter: per-trial
   acceptance is evaluated on the union of 3x3 patches of the pixels
   switched so far plus the trial pixel (not the whole event
   neighbourhood, which made large events intractable and was never in
   the frozen text). SOLVED additionally requires a FINAL engine census
   over the union of all event-row patches plus all switched patches,
   showing zero event rows and conditions (b)-(c) against that same
   union's baseline.
3. SKIP RULE (replaces Addendum 1.2): preserve original row and pixel
   order; skip a candidate pixel only when none of its at-most-four
   incident quads is an endpoint of a currently unresolved target row
   (exact incidence argument; verdict-safe).
4. OUTCOMES: SOLVED / UNSOLVED_EXHAUSTED / CENSORED_TIMEOUT_OR_CAP
   (the 64-switch cap is censoring too). G2 bounds per diagonal:
   lower = solved/all; upper = (solved+censored)/all. lower >= 0.30
   passes; upper < 0.30 fails; otherwise indeterminate -> a
   pre-registered rescue on censored events only. Engine errors ABORT.
5. G1, tightened: an ACTIONABLE contested pixel has >= 2 DISTINCT
   float32 XYZ tuples (identical re-covers are not alternatives); both
   raw and actionable counts reported; G1 uses actionable.
6. HYGIENE: per-subprocess timeouts; checkpoint bound to sha256 of the
   export bands, candidate stream, census CSVs and procedure version;
   candidate-stream schema and bounds validated; assert export shape.
7. SCOPE STATEMENTS for any positive: per-event, per-diagonal
   existence only — no claim of a single mutually compatible global
   assignment; acceptance is local and the eventual selector still
   needs global obstacle discovery plus a full reloaded census.
8. The 15-minute budget stands (no monster exemptions).

## ADDENDUM 3 — 2026-08-06, after review of the first Addendum-2 run

The reviewer found accept_ok exempted target keys from the no-increase
rule, permitting INTERMEDIATE resurrection of resolved target rows,
contradicting Addendum 2. The 12 SOLVED witnesses stand (the final
global verification is strict), but UNSOLVED outcomes are not provably
conservative. Therefore:

1. accept_ok applies the no-increase rule to EVERY key (target keys
   included); the strict-decrease requirement on unresolved targets is
   additional. Checkpoint binding extended to y.tif, z.tif and both
   census CSVs. The Addendum-2 run's outcomes are superseded; a fresh
   primary rerun decides G2.
2. ONE pre-registered grazing sensitivity, run after the primary,
   clearly labelled, cannot qualify the selector: identical procedure
   except acceptance allows NEW GRAZING keys (transverse and coplanar
   still forbidden). Purpose: determine whether the engine's plane-sign
   grazing semantics (which can class AABB-overlapping but disjoint
   pairs as grazing; fixed upstream by distance confirmation) caused
   the negative. It answers that question only.
3. Publication wording (frozen): "Under the preregistered monotone
   single-pixel selector and contact-non-regression policy, only X% of
   events were resolved. This rejects that selector; it does not rule
   out coordinated regional cover optimization." The broader claims
   ("cover choice is not fixable", "rasterization causes the defect")
   are NOT made.
