# Case study: transverse self-intersection in ScrollFiesta's tifxyz strip export

ScrollFiesta ([Hob3rMallow/scrollfiesta_public](https://github.com/Hob3rMallow/scrollfiesta_public),
Best of Month, July 2026) is the strongest published surface-construction
pipeline for Herculaneum scrolls. This case study reproduces its
documented pipeline end-to-end on its own canonical PHerc0139 4×5×5
fixture, measures one topology property of its optional
`scroll_unroll --export-tifxyz` strip output that its metrics do not
currently track, and demonstrates a certified, opt-in cleanup
transaction verified by volume-cartographer's own merged validator.

Every measurement below was pre-registered before it ran
([`prereg/`](prereg)), and every number regenerates from the released
artifacts ([reproduce](#reproducing)). This work identifies a missing
post-export invariant rather than contradicting ScrollFiesta's existing
metrics: ScrollFiesta already marks the whole-grid artifact unverified
and records contested coverage; windcheck adds an exact check on the
resulting tifxyz connectivity.

## The finding

> On ScrollFiesta's documented, explicitly unverified
> `scroll_unroll --export-tifxyz` output for its canonical PHerc0139
> 4×5×5 fixture, windcheck found 25,939 and 25,072 transverse
> triangle-pair contacts under the two grid triangulations, clustering
> into 327 region-pair events under either triangulation. Of the raw
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

The faithful-reproduction baseline: ScrollFiesta's own audit gate
passes on our run (turn-off 2.92% vs their documented 2.98%), so the
crossings coexist with a passing native quality gate — the two measure
different failures.

## Stage-by-stage counts

Cumulative pipeline stages, each exported and censused identically
(rows = transverse triangle-pair contacts; events = 8-connected
region-pair clusters; d0/d1 = the two grid triangulations):

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
configuration. Some defects exist at stage 0; the events are not
independently adjudicated physical errors.

Event character in the final export: median intersection-segment length
0.225 voxel (p90 0.522, maximum 16.9), median crossing angle 68.4°,
spread over 62 of 99 cubes. The likely mechanism is first-cover-wins
rasterization of contested coverage: the exporter itself reports 21.5%
of valid pixels contested (a later cover more than two voxels from the
first-cover coordinate).

## Where the property lives — and where it cannot

The placed per-cube stage is structurally incapable of showing this:
all 99 placed cube meshes lie strictly inside their own 128³ carve
boxes, and registration changes only the parametrization (`vt`), never
the placed world coordinates. A pre-registered gate over all 927
admitted registration edges accordingly found zero 3D contacts — with
zero of 927 edge pairs sharing any 3D bounding-box overlap
([`results/gate_edges.json`](results/gate_edges.json)). The property
first becomes measurable where the solved UV layout is rasterized into
one shared grid.

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
  (9,171 / 8,886 rows, 187 / 179 events; the worst single wrap carries
  86 events).
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
- transaction: 1,424 of 5,212,978 pixels excised — 99.685% retained —
  and the output reloads clean through the official validator.

The finding replicates on a disjoint window of the same sample. Other
samples remain untested.

## The transaction, end to end

The adapter ([`bench/fiesta_adapter.py`](../../bench/fiesta_adapter.py))
consumes the masked export directly, runs the unchanged frozen excision
policy, preserves ScrollFiesta's grid-aligned sidecars
(`provenance.tif`, `winding.tif`, the face-provenance map) bit-for-bit
at retained pixels and invalidates them at excised pixels, hashes every
emitted file, and refuses to hand over unless retained pixels are
byte-identical everywhere and the shipped output censuses clean.

Independent verification through the **official pipeline's own code**:
`vc_tifxyz_selfcross`, merged into volume-cartographer
([ScrollPrize/villa#1303](https://github.com/ScrollPrize/villa/pull/1303)),
reproduces the counts on the unmodified export exactly
(25,939 / 25,072, clean flag `false`) and certifies the transformed
output (0 / 0 across 36.4 billion tested pairs, clean flag `true`).
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

See [`reproduce/REPRODUCE.md`](reproduce/REPRODUCE.md): one command per
fixture regenerates the staged exports, censuses, event clustering,
provenance splits, transaction, and the official reload. Heavy
artifacts (original and transformed exports, face-provenance maps, raw
pair rows) are release assets; small records live in
[`results/`](results). The instrumentation used (an env-gated,
byte-neutral per-pixel first-cover face map in `tifxyz_export.c`, plus
a read-only registration-graph dump in `scroll_whole.c`) ships as a
patch against their MIT-licensed source in [`patches/`](patches),
with two small portability fixes included.
