# Every published trace

The open-data bucket publishes 46 samples. Twelve of them carry surface
traces as `tifxyz`; two more publish segment directories with no `tifxyz`
in them at all. **All 278 traces in those twelve samples have now been
censused** for non-adjacent transverse self-intersection under both
canonical triangulations.

Every figure here comes from `results/corpus/summary.json`, which is
generated rather than written.

## What was added

The audit released as `v0.1.0-beta` covered five samples and 179
censusable traces. This adds seven samples and 99 traces, with zero
fetch or census failures.

| sample | traces | self-intersecting | median valid cells |
|---|---:|---:|---:|
| PHerc0500P2 | 39 | 12 | 47,970 |
| PHerc0009B | 18 | 6 | 49,650 |
| PHerc1447 | 15 | 3 | 25,662 |
| PHercMANBp | 11 | 3 | 16,316 |
| PHerc0343P | 8 | 4 | 113,888 |
| PHerc0800 | 6 | 1 | 15,552 |
| PHerc0332 | 2 | 1 | 205,900 |

### Why these were not covered before

They publish under a second layout. The five original samples put the
mesh at `<segment>/mesh/<name>-on-<volume>-<res>.tifxyz/`; these seven
put it at `<segment>/mesh/intermediate/tifxyz_original/`. A filter
looking for a `.tifxyz` suffix reports zero traces for them, which is
exactly what happened until it was checked.

`tifxyz_original` is traced geometry in volume space and is censused.
`tifxyz_normalized`, which sits beside it in some samples, is output of
the flattening step rather than traced geometry, and is excluded — the
same rule this repository has applied since the first census.

## The rate depends mostly on size

The new samples self-intersect at 30% against 86% in the five original
ones. Read naively that says the new samples are three times cleaner. It
does not.

Their median trace is **5.8× smaller** — 49,650 valid cells against
286,784. Matched by size the difference largely disappears, and in the
band above 150,000 cells the new samples are if anything worse:

| valid cells | original samples | new samples |
|---|---:|---:|
| 5–20k | 0% (n=1) | 3% (n=35) |
| 20–50k | 50% (n=4) | 25% (n=12) |
| 50–150k | 72% (n=40) | 46% (n=41) |
| 150–400k | 91% (n=67) | **100%** (n=7) |

### A one-parameter model fits the whole corpus

If every valid cell independently carried a crossing with probability
`q`, then a trace of `N` cells would self-intersect with probability
`1 − (1−q)^N`. Fitting that single parameter over all 278 traces gives
**q = 7.2 × 10⁻⁶ per valid cell**, and the observed rates track it:

| valid cells | traces | observed | independent-cell null | difference |
|---|---:|---:|---:|---:|
| <5k | 4 | 0% | 1% | −1 |
| 5–20k | 36 | 3% | 8% | −5 |
| 20–50k | 16 | 31% | 23% | +9 |
| 50–150k | 81 | 59% | 50% | +9 |
| 150–400k | 74 | 92% | 81% | +11 |
| >400k | 67 | 93% | 99% | −6 |

Within about ten points across five orders of magnitude of surface size.

**So the size effect is mostly arithmetic, and that is the honest
reading.** A larger surface has more opportunities to fold through
itself; a small trace is clean substantially because it is small. This
is not a claim that any tracer is better or worse than another, and the
86% figure published for the original five samples is not a property of
those samples — it is what happens when you trace large surfaces.

It also gives a usable prior: a trace of 200,000 valid cells has roughly
a 4-in-5 chance of self-intersecting somewhere, before anyone looks at
it.

## Reproducing it

```sh
uv run python bench/corpus_expand.py plan   # discover samples and layouts
uv run python bench/corpus_expand.py run    # fetch and census
```

The plan step records, per sample, the layout, the chosen volume and the
volume ids seen. A census is only meaningful against one volume, because
the grid indices it emits are indices into the volume it read.

## Limits

- **Self-intersection is a subset of tracing error.** A surface that
  leaves the correct sheet and never returns need not self-intersect.
  Nothing here detects that case, and one attempt to close that gap is
  recorded as a null in [`PATCH-AUDIT.md`](PATCH-AUDIT.md).
- **These 99 traces are censused, not transformed.** The released
  transverse-clean outputs still cover the original 179.
- **Two samples publish no `tifxyz`.** That is a statement about what is
  published, not a gap in this audit.
- The census reports where a surface meets itself. It does not say which
  branch is wrong, and it establishes no cause.
