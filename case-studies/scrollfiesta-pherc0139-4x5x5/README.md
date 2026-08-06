# Case study: transverse self-intersection in ScrollFiesta's tifxyz strip export

ScrollFiesta ([Hob3rMallow/scrollfiesta_public](https://github.com/Hob3rMallow/scrollfiesta_public))
won Best of Month in July 2026 for automatic surface construction on
Herculaneum scrolls. This case study reproduces its
documented pipeline end-to-end on its own canonical PHerc0139 4×5×5
fixture, measures one topology property of its optional
`scroll_unroll --export-tifxyz` strip output that its metrics do not
currently track, and demonstrates a certified, opt-in cleanup
transaction verified by volume-cartographer's own merged validator.

The primary gate, intervention, ownership/atlas controls, resolution
control and held-out replication were pre-registered before they ran
([`prereg/`](prereg)), and every number regenerates from the released
artifacts ([reproduce](#reproducing)). This work identifies a missing
post-export invariant rather than contradicting ScrollFiesta's existing
metrics: ScrollFiesta already marks the whole-grid artifact unverified
and records contested coverage; windcheck adds a deterministic check for
transverse self-intersection in the resulting tifxyz geometry.

## The finding

> On ScrollFiesta's documented, explicitly unverified
> `scroll_unroll --export-tifxyz` output for its canonical PHerc0139
> 4×5×5 fixture, windcheck found 25,939 and 25,072 transverse
> triangle-pair contacts under the two grid triangulations, clustering
> into 327 region-pair events under each triangulation. Of the raw
> contacts, 17,832 and 17,470 joined quads whose four-corner
> first-cover `(cube,gid)` label sets were disjoint. This establishes
> intersecting topology in the shared-grid export, not an error in any
> particular registration edge or a physically incorrect sheet
> assignment.

And the transaction:

> Applying the preregistered frozen excision policy produced a
> reload-verified derivative with zero transverse contacts under both
> triangulations while retaining 99.443% of canonical represented area.
> The operation increased connected components from 117 to 140 and left
> all retained coordinates bit-identical. It demonstrates a post-export
> topology-cleaning transaction for this artifact; it does not
> establish physical correctness or repair ScrollFiesta's winding
> solver.

The faithful-reproduction baseline: 99 of 100 cubes meshed (one cube
crashed deterministically under Linux/GCC but succeeds in
ScrollFiesta's Windows reference run), and ScrollFiesta's own audit
gate passes on our run (turn-off 2.92% vs their documented 2.98%), so
the transverse contacts coexist with a passing native quality gate — the
two measurements address different properties.

## Stage-by-stage counts

Cumulative pipeline stages, each exported and censused identically
(rows = transverse triangle-pair contacts; events = unordered pairs
of 8-connected participating-quad regions; d0/d1 = the two grid
triangulations):

| stages | rows d0 | rows d1 | events d0 | events d1 |
|---|---:|---:|---:|---:|
| `0` (placed UVs as-is) | 11,232 | 10,646 | 116 | 109 |
| `12` (bake + weld/relax) | 21,651 | 20,947 | 261 | 265 |
| `123` (+ ownership, diagnostic) | 21,651 | 20,947 | 261 | 265 |
| `1234` (+ CT snap) | 23,741 | 22,681 | 286 | 284 |
| `12345` (+ final relax) | 25,939 | 25,072 | 327 | 327 |

On this fixed fixture, export resolution and configuration, the
measured transverse-contact and clustered-event counts increased after
each geometry-changing cumulative stage. Step 3 produced byte-identical
geometry because SeamOwn is diagnostic-only under the default
configuration. Some transverse contacts already exist at stage 0; the events are not
independently adjudicated physical errors.

Among label-disjoint raw contacts in the final export, the median
intersection-segment length is 0.225 voxel (p90 0.522, maximum 16.9),
the median crossing angle is 68.4°, and the participating sources span
62 of 99 cubes. The exporter uses first-cover-wins rasterization and
reports 21.5% of valid pixels as contested, providing a candidate
mechanism for investigation. These measurements do not establish
causation, identify the correct cover, or implicate a particular
registration edge.

## Where the property lives — and where it cannot

Cross-cube 3D contacts cannot appear in the placed representation: all
99 placed cube meshes lie strictly inside their own 128³ carve boxes,
and registration changes only the parametrization (`vt`), not their
world coordinates. A pre-registered gate over all 927 admitted
registration edges therefore found zero 3D contacts, with zero edge
pairs sharing any 3D AABB overlap
([`results/gate_edges.json`](results/gate_edges.json)). The shared-grid
export is the first tested artifact in which the solved UV layout
becomes one connected tifxyz surface.

## Controls

- **Ownership applied** (`--apply-ownership`, the legacy opt-in;
  verified applied, piece set compacted 189,841 → 185,852 faces): rows
  36,932 / 35,859, events 513 / 534 — higher than the default. One
  topology invariant was measured, not overall surface quality.
- **Per-wrap atlas export** (`--export-atlas`, the per-(wrap × z-slab)
  multi-segment representation; ScrollFiesta's documented terminal
  *stage* is `atlas_ribbon_fit`, which was not tested here): of 47
  exported segments, 11 fall below the census validity floor and 18 of
  the remaining 36 contain within-segment transverse crossings
  (9,171 / 8,886 rows, 187 / 179 events; the worst atlas segment carries 86 d0 events and 80 d1 events).
- **Resolution** (`du = dv = 2`, quarter the pixels): 11,309 / 10,937
  rows, 182 / 174 events. Counts are raster-resolution dependent; the
  phenomenon is robust.
- **Planted positive control**: a +0.5 voxel shift produces detections
  on top of baseline, so the instrument detects at this scale.

## Held-out replication

Pre-registered before any carve on a z-adjacent, voxel-disjoint 4×5×5
window of the same sample (one pre-named fallback, not needed):

- their audit gate passes there too (turn-off 3.48%);
- final export: 11,585 / 10,421 rows, 273 / 278 events;
- transaction: 1,424 of 5,212,978 valid pixels excised, 99.685% of
  represented area retained, and the output reloads clean through the
  official validator.

The finding replicates on a disjoint window of the same sample. Other
samples remain untested.

## The transaction, end to end

The adapter ([`bench/fiesta_adapter.py`](../../bench/fiesta_adapter.py))
consumes the masked export directly, runs the unchanged frozen excision
policy, preserves ScrollFiesta's grid-aligned sidecars
(`provenance.tif` and the face-provenance map) bit-for-bit
at retained pixels and invalidates them at excised pixels, hashes every
emitted file, and refuses to hand over unless retained pixels are
byte-identical everywhere and the shipped output censuses clean.

Independent verification through the **official pipeline's own code**:
`vc_tifxyz_selfcross`, merged into volume-cartographer
([ScrollPrize/villa#1303](https://github.com/ScrollPrize/villa/pull/1303)),
reproduces the transverse counts on the unmodified export exactly
(25,939 / 25,072, clean flag `false`) and reports the transformed
output as transverse-clean (0 / 0 across 36.4 billion tested pairs, clean flag `true`).
Reports: [`results/`](results).

## Scope, stated plainly

- All claims concern topology of the exported tifxyz grid — never the
  physical correctness of any sheet, and never a defect in any
  particular registration decision.
- The tested strip export is documented but optional, and explicitly
  unverified in ScrollFiesta's own writeup; `atlas_ribbon_fit`, their
  documented terminal stage, was not tested.
- Two fixtures of one sample were measured; nothing here extrapolates
  beyond them.
- Cleaned means validator-clean under the stated validator, at a
  measured retention cost, with fragmentation reported
  ([`results/transaction_canonical.json`](results/transaction_canonical.json)).

## Reproducing

See [`reproduce/REPRODUCE.md`](reproduce/REPRODUCE.md): the documented command sequence per
fixture regenerates the staged exports, censuses, event clustering,
provenance splits, transaction, and the official reload. Heavy
artifacts (original and transformed exports, face-provenance maps, raw
pair rows) are release assets; small records live in
[`results/`](results). The instrumentation used (an env-gated,
byte-neutral per-pixel first-cover face map in `tifxyz_export.c`, plus
a read-only registration-graph dump in `scroll_whole.c`) ships as a
patch against their MIT-licensed source in [`patches/`](patches),
with two small portability fixes included.

## Follow-up: can re-selecting the discarded covers fix it? (negative)

A pre-registered applicability experiment
([`prereg/PREREG-COVER-APPLICABILITY.md`](prereg/PREREG-COVER-APPLICABILITY.md),
with its addenda) asked whether choosing among the pipeline's own
discarded candidate covers could remove the crossings. Of 327 events
per triangulation, 314 (96.0%) on d0 and 311 (95.1%) on d1 touched a
corner pixel for which the pipeline generated at least two distinct
float32 XYZ candidates but retained only one. Under the preregistered
monotone single-pixel search, with no increase in any contact class
and no loss of retained quads, 5 d0 events (1.5%) and 7 d1 events
(2.1%) were resolved; one event per triangulation was censored, giving
upper bounds of 1.8% and 2.4% against the frozen 30% gate. This
rejects that selector; it does not rule out coordinated regional cover
optimization. A labelled sensitivity permitting new grazing contacts
also failed, with an upper bound of 2.1% under both triangulations.
Events were solved independently; d0 and d1 solutions are not
necessarily mutually compatible; G1 establishes spatial coincidence
between alternatives and event corners, not rasterization causation.

## The production remedy: `windcheck transaction`

The transaction demonstrated above is now one command:

```sh
uv run windcheck transaction candidate.tifxyz --out final.tifxyz \
    [--adapter scrollfiesta] [--official-validator vc_tifxyz_selfcross] \
    [--report report.json]
```

It stages the input, transforms only if needed under the unchanged
frozen policy, verifies (retained pixels byte-identical, sidecars
preserved/invalidated, reload census clean under both triangulations,
optionally the official validator), writes a hash-bound report, and
promotes atomically. Unknown sidecars are refused, never guessed.
Stable exit codes: 0 already-clean committed, 10 transformed-and-clean
committed, 3 refused, 2 invalid input, 1 internal.
