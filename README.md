# windcheck

A label-free consistency check for traced papyrus surfaces in the
[Vesuvius Challenge](https://scrollprize.org) open data.

It answers one question about a segmentation, using nothing but the mesh:

> **After going once around the scroll, does this trace land on the next sheet —
> or back on the one it just traced?**

No ground truth, no annotations, no ML, no GPU, no volume download. The whole
analysis runs on ~4 MB of `tifxyz` per segment, in seconds, on a laptop.

Targets [Open Problems](https://scrollprize.org/2026_open_problems) §2 and §3 —
*"sheet switches, where a traced mesh jumps from the intended sheet to a
neighboring wrap"* and mesh topology repair.

---

## Quick start

```sh
uv sync
clang++ -O3 -std=c++17 -pthread -o engines/atlas_query engines/atlas_query.cpp
uv run pytest -q                       # 5 tests, no data needed

# fetch one scroll's surfaces (~670 MB) and analyse every trace in it
uv run python -m windcheck.fetch --sample PHerc0172
uv run windcheck selfgap data/scroll5_tifxyz --json report.json
```

Output, one row per segment:

```
                             trace      u  <2.0vx  flag%     blob  blob%  top5   du p10   valid
20251115002745-auto_grown_20251115   3065    6.60   14.5   13,807   7.02   81%      500      OK
20251115002747-auto_grown_20251115   2817    5.71    8.3    2,555   1.52   46%      569      OK
20251115002740-auto_grown_20251115   2460    0.36    0.9       78   0.06   25%      209      OK
20250926114310-w080_20250926114310    571    2.30    4.8       16   2.84  100%       60  REJECT
```

`blob%` is the number to compare across traces. `REJECT` means the validity
filter refused to trust the row (see below) — it is not a defect finding.

---

## How it works

A multi-wrap trace is stored as a 2D grid of 3D points: grid index `(v,u)` maps
to a volume coordinate. Walking along `u` goes around the scroll.

For every point, measure the distance to the nearest **other part of the same
trace**, ignoring anything within `exclude_u` grid columns of the point's own
`u`. What survives that exclusion is the trace's **neighbouring wrap**.

```
   correct                              defective
   ─────────────────────────            ─────────────────────────
   wrap N+1  ───────────────            wrap N+1  ──────╮
                ↕ ~17 vx                                ╰──────  ← ~0 vx:
   wrap N    ───────────────            wrap N    ───────────────   same sheet
                                                                    traced twice
```

This never needs a scroll axis or a radius, which matters: a Herculaneum scroll
is crushed, not round, so anything built on "distance from the centre" fails
(see [notes/DECISIONS.md](notes/DECISIONS.md) — a radius fit has residual sd
79.8 vx against a 12–17 vx sheet spacing, on ground-truth sheets).

### Calibration

Measured on 126 pairs of the 44 labelled single-winding segments of Scroll 5:

| true separation | measured | ratio |
|---|---|---|
| 1 sheet | **17.0 vx** | 1.00 |
| 2 sheets | 34.1 vx | 2.00 |
| 3 sheets | 51.2 vx | 3.01 |

Independently re-derived during review as 16.75 / 33.76 / 49.87 (ratios 2.015 /
2.977). Per-*point* discrimination between a 1- and 2-sheet gap is only **77.8%**,
so nothing is claimed per point — verdicts are per region.

### The 2.0-voxel threshold is not ours

`volume-cartographer/core/src/GrowSurface.cpp:121` sets `same_surface_th = 2.0`,
and `:2919-2931` rejects any growth candidate landing within that distance of
already-traced surface (*"candidate rejected: nearest same-surface location had
no valid state"*). `<2.0vx` counts points in published output sitting inside
that same distance.

---

## Controls

**Null control.** A single-wrap segment has no previous wrap, so a correct
detector must find nothing. All 44 labelled windings of Scroll 5 report
**0.00% flagged**. `tests/test_geometry.py` pins this on a synthetic flat sheet.

**Validity filter.** Exclusion by `u`-index assumes `u` advances with angle. If a
trace's flagged partners sit at `|Δu|` barely above the cutoff, they may be
*same-wrap* neighbours and the row is an artifact. Any trace whose flagged
`|Δu|` 10th percentile falls within 2× the cutoff is rejected in code, not by eye.

It works: on Scroll 5 it rejects four labelled sheets that the raw detector fires
on, and on PHerc0814 it rejected the two *highest-scoring* traces (11.12% and
8.18%) — which would otherwise have been headline findings.

**Cross-check.** Flagged partners sit at `|Δu|` = 271–618 columns, matching
revolution periods derived independently by angular unwrapping (236–616) across
all nine traces.

---

## What this does NOT claim

Stated plainly, because three stronger versions of this were tested and died:

- **Not proven to be errors.** Three CT experiments failed to show flagged
  regions differ from clean ones (best: Cohen's *d* = 0.245, p = 0.0496 — too
  small to certify anything). The geometry is solid; the interpretation is not.
- **Not proven to violate the grower's criterion.** The published Scroll 5
  meshes trace back through private `s3://philodemos/...` artifacts, flattening
  and reconstruction. Nothing establishes they came from the audited
  `GrowSurface.cpp` revision. What *is* true: current VC code runs no final
  self-contact check after `optimize_surface_mapping`.
- **Not a large out-of-sample result.** The nine Scroll 5 `auto_grown` "traces"
  are partitions `_0`…`_8` of a *single* artifact
  (`auto_grown_20251115002740308`). Effective n is 1 run plus 8 PHerc0814
  segments — not 19 independent traces.
- **Violation *rate* is not distinctive.** It overlaps between scrolls
  (S5 0.36–6.60%, PHerc0814 0.32–6.77%). Only *concentration* separates them,
  and normalised that is single-digit, not the 58× an early draft claimed by
  mixing thresholds.

The surviving claim is deliberately small: **scattered self-coincidence is
normal; concentrated self-coincidence is not, and it ranks traces for review.**

---

## Limitations

- `exclude_u` is fixed. It should adapt to each trace's estimated turn period.
- Exclusion is by grid column, not mesh geodesic distance.
- Two scrolls is not enough for a pipeline-wide claim.
- Single-wrap patches (u-extent < 4·`exclude_u`) cannot be analysed at all.

## Layout

```
engines/atlas_query.cpp   point-to-surface distance kernel, ~27k queries/s (12 threads)
src/windcheck/
  tifxyz.py               reader, written from the format spec
  atlas.py                surface assembly + binary bridge to the engine
  selfgap.py              the analysis and the validity filter
  ct_check.py             the CT confirmation attempt (which failed)
  catalog.py, fetch.py    open-data catalog, download, SHA-256 manifest
notes/DECISIONS.md        append-only log of every decision and every kill
```

`tifxyz.py` deliberately does not use the upstream `vesuvius` package: this tool
audits that pipeline, so it must be able to disagree with it.

## Licence

MIT.
