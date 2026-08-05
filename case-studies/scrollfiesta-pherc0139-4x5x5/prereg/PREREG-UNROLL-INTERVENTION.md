# Pre-registration: post-export transaction on Fiesta's final tifxyz

Written 2026-08-05, AFTER the unroll-export gate returned POSITIVE
(PREREG-UNROLL-GATE.md; 17,832/17,470 cross-source transverse events
d0/d1 in the final export) and BEFORE any intervention ran. This is
level 1 of the two-level intervention plan accepted in the external review: the
ALREADY-FROZEN windcheck transform (round-28 policy, unchanged: same
validator, thresholds, selection rules, invalidation, budgets) applied
to the final export as a post-export transaction. Level 2 (GroupGraph
edge intervention) is NOT covered here and would need its own prereg.

## Operation

Input: the final staged export (`unroll_s12345/tifxyz_run`), pinned by
a one-entry corpus_bases/v1 manifest (hashes of x/y/z/mask; voxel_um
9.362, from the source volume name `20250728140407-9.362um-...`,
recorded as name-derived). Operator: `bench/excise_shadow.py
--certificate` under `FROZEN_POLICY` — not one parameter may change.
The 600 s process limit applies; if the solve exceeds it, that IS the
result (the policy's scheduling, reported as such).

## Endpoints (fixed now)

- PRIMARY: post-reload census of the transformed export = ZERO
  transverse under BOTH canonical triangulations (the transform's own
  acceptance bar; `residual_transverse` is a valid, reportable outcome).
- Preservation: retained fraction of canonical area >= 0.99 for the
  STRONG claim (the external review bar). Below 0.99 the result is still
  reported in full; the claim downgrades to "a validator-clean
  derivative exists at X% retention".
- Reported regardless: fragmentation (components before/after, core
  gate), cut boundary length, solver status/method mix, wall time.
- Fiesta's own audit gate is UPSTREAM of the export and untouched by a
  post-export transform; it is recorded as unchanged by construction,
  not claimed as re-verified.

## Outcomes

- STRONG: transverse 0/0 after reload AND retained >= 0.99. Claim:
  "the final exported surface of the champion pipeline's canonical run
  contains ~18k assembly-attributable transverse crossings per
  triangulation; a frozen, certified transform removes ALL of them for
  under 1% of represented area, as a transaction the pipeline could
  adopt."
- PARTIAL: transverse 0/0 but retained < 0.99, or residual_transverse.
  Report numbers; no "improvement" framing beyond what they support.
- FAIL: the operator errors or times out; reported as a limitation of
  the frozen policy at this grid size.

## What does NOT count

- Tuning any transform parameter for this mesh.
- Removing the same-source events from the claim denominator: the
  transform acts on ALL transverse events; cross/same-source split
  stays a diagnostic of the gate, not a filter of the fix.
- Any statement that the removed geometry was physically wrong papyrus.
