# Pre-registration: the unroll-export gate on Scroll Fiesta's final tifxyz

> **Post-run terminology correction:** In the frozen text below,
> "event" was sometimes used for one reported triangle-pair contact.
> Public results reserve "contact" or "row" for that unit and
> "region-pair event" for the subsequent clustering. Accordingly,
> 17,832/17,470 are label-disjoint raw contacts, not clustered events.
> Label-disjoint provenance does not by itself establish a causal
> assembly error, a faulty registration edge, or physical
> incorrectness. The preregistered positive/null outcome is unchanged.


Written 2026-08-05, BEFORE any export was produced or censused, per house
rule. Follows the Fiesta placed-frame experiment (PREREG-FIESTA-GATE.md),
which concluded in a pre-registered NULL with a structural explanation:
cross-cube 3D interpenetration is impossible by construction in the
placed representation (meshes box-clipped; registration is
parametrization-only). an external review round ranked this follow-up first;
approved 2026-08-05.

## Hypothesis

Scroll Fiesta's FINAL exported surface -- the tifxyz written by
`scroll_unroll --export-tifxyz`, which rasterizes the solved UV layout
into one regular grid by first-cover-wins after seam welding, ownership,
CT snapping and relaxation -- contains transverse self-intersections,
and some of them are CROSS-SOURCE: the two participating grid cells were
covered by faces from different (cube, winding-group) sources, making
them assembly artifacts rather than inherited per-piece geometry.

## Why this is the right surface

The registration's decisions (GGEdges) first become ONE testable
geometry in this export: individually disjoint cube meshes are placed
onto a shared (u, v) canvas, and pixels from different cubes/wraps
become grid-adjacent samples of a single claimed surface. This is the
artifact the pipeline actually hands to VC3D.

## Instrumentation (before any measurement)

Env-gated tap in `tifxyz_export.c` (mirroring the SF_GG_DUMP pattern):
`SF_TXZ_FACEMAP=<prefix>` writes, as NEW sidecar files only,

- `<prefix>_facemap.i32`: int32[W*H], the piece-set face index whose
  first cover wrote each pixel; -1 = empty;
- `<prefix>_faces.tsv`: per kept face f: cube id (`ids[face_cube[f]]`)
  and the face's winding-group id under the FROZEN attribution rule from
  the placed experiment (gid iff all three vertex gids equal and >= 0,
  else -2), using the piece set's own per-vertex `gid` array.

Proof obligation before any census is read: with the env var unset the
instrumented binary's outputs are byte-identical to the uninstrumented
build's; with it set, every STANDARD output (x/y/z/mask/meta/winding/
provenance TIFs, stats, logs excluded) is byte-identical to the
env-unset run.

## Protocol

1. Input: the reproduced PHerc0139 4x5x5 fixture, day-1 placed dir
   (`output/run_placed`), untouched. Raw cubes: `PHerc0139-4x5x5/
   cubes_RAW`. Their canonical flags (README): `--raw <cubes_RAW>
   --id run`, all other parameters at their defaults (du = dv = raster
   default; no tuning).
2. Five cumulative-stage runs, each to its own out dir, each with
   `--export-tifxyz` + the facemap tap:
   `--steps 0` (placed UVs as-is), `12`, `123`, `1234`, `12345`
   (1 bake, 2 UV-weld + banded relax, 3 SeamOwn ownership, 4 dark-vertex
   snap, 5 final relax). Stage geometry accumulates; `--steps 0` is the
   no-stage baseline export.
3. Census EVERY export as a tifxyz grid with the corpus-frozen
   parameters: both diagonals, exclude=1, cell=40.0, maxedge=60.0,
   touch_tol=1e-3 (`engines/selfcross`, the sole authority).
4. Attribution of a census event (frozen now): a grid cell (v, u) has
   source set = { (cube, gid) of facemap at its four corner pixels,
   omitting empty corners }. An event (two cells) is
   - CROSS-SOURCE iff the two source sets are non-empty and DISJOINT;
   - SAME-SOURCE iff the sets intersect;
   - UNATTRIBUTED iff either set is empty (reported, never dropped).
5. Facemap consistency check: facemap == -1 exactly where the export's
   own provenance.tif == 0 (empty). Any disagreement stops the analysis
   as an instrumentation bug.
6. Positive control (already demonstrated at the placed stage, repeated
   here at grid level): a copy of the final export with a +0.5 vx
   x-shift planted into a band of columns must produce transverse
   events under the same census parameters.

## Endpoints

- PRIMARY: the number of CROSS-SOURCE transverse events in the FINAL
  (`--steps 12345`) export.
- Secondary (descriptive, no claims hang on them): total transverse per
  stage export; the first stage at which each final-export event's cell
  pair is transverse; same-source and unattributed counts; the export's
  own multi/conflict pixel statistics; valid-pixel coverage per stage.

## Outcomes (fixed now)

- POSITIVE: >= 1 cross-source transverse event in the final export.
  Claim allowed: "the final exported surface of the champion pipeline's
  canonical run contains transverse self-intersections created by
  assembly (cross-source)". The INTERVENTION (post-export windcheck
  transform as a transaction; GroupGraph edge intervention only where an
  event maps unambiguously to an admitted edge) gets its OWN
  pre-registration, written only after this outcome, before it runs.
- SAME-SOURCE-ONLY: transverse events exist but none cross-source.
  Claim allowed: "the export inherits/creates per-piece topology
  defects; assembly-attributable crossings were not observed". The
  post-export transform intervention remains available; edge
  intervention is NOT justified.
- NULL: the final export is transverse-clean. Published briefly as a
  null; the placed-frame structural finding plus this null closes the
  Fiesta line for August. No second fixture, no parameter-space
  analysis, no threshold tuning.

## Kill criteria (fixed now)

- 3 working days from 2026-08-05 for the applicability answer (through
  step 6). If the pipeline cannot produce the five exports, or the
  byte-neutrality proof fails and cannot be fixed inside the budget,
  stop and record.
- No selection: every export censused, every event classified, all
  counts reported including zeros and unattributed.

## What does NOT count

- Gating anything Fiesta's canonical run did not produce.
- Any claim about ink, text, or physical sheet correctness.
- Extrapolation beyond this fixture without qualification.
- Re-defining attribution after seeing events (the rule above is
  frozen; if it proves inadequate the inadequacy is REPORTED).
