# Pre-registration: Fiesta's own controls on the unroll-export finding

Written 2026-08-05, BEFORE any control run, per house rule. the external review
rank-1 outstanding control on the POSITIVE unroll-export result
(PREREG-UNROLL-GATE.md: 25,939/25,072 transverse rows, 327/327
region-pair events in the final default export): a reviewer will ask
whether Fiesta's own ownership machinery, or its collision-safe atlas
terminal, already avoids the crossing load.

## Facts pinned from source (before running)

- The DEFAULT canonical path runs SeamOwn diagnostic-only; `--own-apply`
  is an explicit LEGACY opt-in ("Legacy opt-in only" in scroll_unroll.c)
  that compacts the piece set to kept faces and applies rehomed UVs, so
  stages 4-5 and the export see ownership-filtered geometry.
- `--export-atlas` splits the registered strip into per-(wrap x z-slab)
  tifxyz segments that tile the strip disjointly by construction, each a
  standalone VC3D segment on a shared lattice (export_atlas.h). This is
  the multi-piece "collision-safe" terminal.
- The facemap tap exists only in the strip tifxyz exporter; atlas
  segments carry no source attribution (their endpoints here are
  counts, not attribution).

## Runs (frozen)

A. `--steps 12345 --own-apply --export-tifxyz` + facemap tap, canonical
   flags otherwise, fresh out dir. Census the export with the
   corpus-frozen parameters (both diagonals, exclude=1, cell=40,
   maxedge=60, touch_tol=1e-3); classify rows by the frozen label rule;
   cluster events with bench/crossing_analyse.py's region-pair
   definition VERBATIM (8-connected quad regions; an event is an
   unordered pair of distinct regions).
B. `--steps 12345 --export-atlas` (atlas defaults, diagnostic step 3 as
   in the canonical path), fresh out dir. Census EVERY exported atlas
   segment the same way; segments below the census validity floor are
   recorded as such, never dropped silently.

## Endpoints

- A primary: transverse rows and region-pair events (d0/d1) in the
  ownership-applied export, against the default baseline
  (25,939/25,072 rows; 327/327 events).
- B primary: the number of atlas segments with >= 1 transverse event,
  and total rows/events across segments.
- Secondary, descriptive: label-disjoint splits (A only), coverage and
  export stats, per-segment tables (B).

## Interpretation branches (fixed now)

- A comparable to baseline: ownership application does not remove the
  crossing load; the default-path finding stands and the legacy path is
  not a fix.
- A clean or dramatically lower: Fiesta's own machinery already
  resolves most of the load when applied; per the R38 kill condition
  the finding is presented as an unverified-export case study, with the
  default-vs-legacy distinction stated plainly.
- B all-clean: the collision-safe atlas terminal avoids the defect; the
  strip-export finding must be scoped to the single-strip export, never
  generalized to "Fiesta's outputs".
- B dirty: crossings survive the per-wrap split; the finding broadens.

## Kill criteria

One working day. A run that cannot complete is recorded as such. No
parameter tuning; no additional stage combinations beyond A and B; no
event-rule adjustment after seeing counts.
