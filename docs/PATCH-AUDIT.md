# Auditing 84,316 published patches

An independent self-intersection census of the `verified_patches` dataset
for PHercParis4, published on `dl.ash2txt.org` on 30 July 2026. Every
number on this page comes from `results/patches/summary.json`, which is
generated rather than written; the per-patch records are attached to the
release.

This is not our data. It was audited because its author published a
specific, checkable claim about it, and because a contact predicate needs
no labels, no threshold and no permission to test one.

## Result

| | |
|---|---|
| Patches censused | **84,316** (the entire published set) |
| Valid cells | 96,715,651 |
| Transverse-clean | **84,311** |
| Self-intersecting | **5** (1 in 16,863) |
| Coplanar contacts anywhere | 0 |
| Failures | 0 |
| Wall time | 40 minutes, one laptop, no GPU |

Both canonical quad triangulations, non-adjacent contacts only, the same
census that produced every other result in this repository.

## What "5" means, and the baseline that makes it mean anything

A count with nothing to compare it against is not a measurement. "We
found almost nothing" is equally consistent with a clean dataset and with
a test that cannot see anything.

Two baselines were computed. The first was rejected: a crossing rate per
valid cell over the published corpus predicts ~916 affected quads in a
sample this size, but crossings cluster heavily, so a per-cell rate
overstates what independent small surfaces would show.

The defensible comparison is **tiling**. Cut the published traces into
tiles the size of a patch, and count a tile as self-intersecting only
when **both** quads of a crossing land inside it, which is exactly the
condition under which a patch could exhibit that crossing at all.

| population | tiles | with a crossing | rate |
|---|---:|---:|---:|
| Scroll 1 only, matching the patch source | 124,856 | 780 | **0.6247%** |
| all five scrolls | 205,188 | 3,475 | 1.6936% |

At the Scroll 1 rate, **527** of the 84,316 patches would have carried a
self-intersection. Five did: about **100× cleaner than published traces**
of the same scroll, cut to the same size.

### The test is not vacuous

A small patch is nearly planar and has limited room to fold through
itself, so a null result could be an artifact of size. That was checked
before any claim was made: across the pinned corpus, **13.93%** of real
transverse crossings have both branches within a 24×24 window. A patch of
this size can exhibit one.

Patch sizes, in valid cells: min 164, median 571, p90 2,780, max 8,184.

## The five, and the single mechanism behind them

| patch | contacts | peak penetration | clean after 1-cell erosion |
|---|---:|---:|---:|
| `band-seed3224943-20260729-182133-568` | 19 | 16.88 vx | 83.5% retained |
| `band-seed7630648-20260729-181254-173` | 8 | 7.58 vx | 77.9% retained |
| `band-seed2570227-20260729-161313-704` | 7 | 8.36 vx | 83.8% retained |
| `band-seed12900803-20260729-195624-897` | 5 | 8.21 vx | 82.2% retained |
| `band-seed11817810-20260730-003948-986` | 2 | 3.49 vx | 84.5% retained |

All five appear under **both** triangulations, so none is an artifact of
one diagonal choice.

**They share one cause.** All 18 participating quads sit within 2.2 cells
of the boundary of the valid region, 13 of them at exactly 1.0, against
patch median depths of 3.6 to 5.0 cells. Grid separations are `dv ≤ 2`
with `du` from 2 to 120: these are local folds along the edge, not two
distant branches of a trace meeting.

The dataset's author had already identified this failure mode when
publishing, describing the last boundary vertex as curling inwards or
outwards and suggesting a boundary erosion. This audit adds three things
that were not previously known: that the curl sometimes folds the surface
**through itself** rather than merely distorting it; that it does so at a
rate of 1 in 16,863; and that **one** cell of erosion is sufficient.

**The remedy was executed, not predicted.** Eroding the validity mask by
a single cell takes all five to zero transverse contacts under both
diagonals, retaining 77.9% to 84.5% of valid cells. No coordinate is
moved: cells are dropped, so a patch that comes back clean is clean by
removal rather than by adjustment. Two cells also works and costs roughly
twice as much surface.

## Reproducing it

```sh
# the full sweep: ~40 min, ~6 GB transferred, resumable
uv run python bench/patch_audit.py --sample 0 --jobs 32 \
    --out out/patches/audit_full.jsonl

# diagnose whatever it flags: boundary depth and the erosion test
uv run python bench/patch_diagnose.py --audit out/patches/audit_full.jsonl

# the tile baseline, from this repository's own corpus
uv run python bench/patch_baseline_tiles.py

# does cleanliness compose: boxes, overlapping pairs, then a sample
uv run python bench/patch_pairs.py bboxes --jobs 32     # ~9 min
uv run python bench/patch_pairs.py pairs                # writes ~6 GB
uv run python bench/patch_pairs.py census --sample 400
```

The pair list is large enough that the census stage reservoir-samples it
in a streaming pass rather than loading it, so memory stays proportional
to the sample rather than the file.

Sampling is deterministic given `--seed` and the cached sorted index, so
a partial run is reproducible without fetching the whole set. Every
record carries the SHA-256 of all four files of its patch, so any single
finding can be rechecked or contested individually rather than having to
trust the aggregate.

`--jobs 32` is the measured plateau; more workers gain nothing and the
server is a shared community resource.

## Does cleanliness compose?

Patches are published to be assembled, and a per-patch verdict says
nothing about that. Two surfaces each free of self-intersection can still
pass through each other, and any merge of such a pair either
self-intersects or must discard part of one.

The census reads one surface at a time, so a pair is tested by stitching
both into a single grid separated by eight invalid rows, far wider than
the adjacency exclusion of one. Because each patch is independently
clean, every contact reported is necessarily *between* them. This is not
a merge: there is no reparametrisation and no seam. It answers the
load-bearing question, which is whether the two occupy the same space
transversally.

The method was validated in both directions before use. A copy of a patch
translated 50,000 voxels away reports **0** contacts; a copy rotated 90°
about its own centroid reports **112**. It neither misses nor invents.

Of 84,316 bounding boxes, **27,778,181 pairs overlap**. Of 459 sampled
and censused, **10 interpenetrate — 2.2%** (95% CI 0.8–3.5%).

### The angle is the finding

| | |
|---|---|
| inter-patch contacts measured | 52,081 |
| median crossing angle | **3.2°** |
| below 10° | **93.6%** |
| above 30° | 0.1% |
| median penetration | ~4.4 vx |
| grazing tolerance | 0.025 vx |

Penetration three orders of magnitude above the tolerance means these are
genuine contacts, not numerical noise. But shallow angle *with* real
penetration is not a fold: it is two nearly parallel surfaces weaving
through each other. The failure mode is **two patches disagreeing about
where the same sheet is, by a few voxels**, not a sheet folding through
itself. For contrast, the intra-patch crossings in the previous section
run 7.8° to 34.7°.

So a naive union of overlapping patches self-intersects roughly 2% of the
time, and what has to be reconciled is positional disagreement rather
than topology.

Detail, including the ten pairs: `results/patches/pairs_summary.json`.

## Limits

- **Self-intersection is a subset of tracing error.** A surface that
  leaves the correct sheet and never returns need not self-intersect.
  This audit says nothing about that case.
- **Bounding-box overlap is a weak proxy** for the pairs a merger would
  actually join, so the 2.2% is over box-overlapping pairs, not over
  merge candidates. The rate among genuine merge neighbours is
  unmeasured and could be higher or lower.
- **The pair test attributes nothing.** Two patches covering one region
  are supposed to be near-coincident; a disagreement between them says
  they are inconsistent, not which one is wrong.
- **The two populations differ by construction.** Patches are
  algorithmically selected; corpus tiles are arbitrary cuts of long
  traces. That difference *is* the effect being measured, but it means
  the comparison shows "cleaner than published traces", not "clean in
  any absolute sense".
- **Patches from one generator are not independent draws**, so the
  suppression factor is descriptive rather than a significance test.
- The census reports where a surface meets itself. It does not say which
  branch is wrong, and it establishes no cause.

## Artifacts

- `results/patches/summary.json` — every figure on this page, generated
- `results/patches/diagnosis.json` — the five, with quads, boundary
  depths, penetrations and the full erosion sweep
- `windcheck-patch-audit-records.jsonl.gz` (release asset) — one JSON
  object per patch for all 84,316, with per-file SHA-256
