# check-pairs: do two surfaces pass through each other?

A per-surface verdict says nothing about assembly. Two surfaces that are
each free of self-intersection can still interpenetrate, and any merge
that keeps both in full contains a self-intersection. Pipelines that build
a graph over patches and join along its edges need to know which proposed
edges are geometrically impossible BEFORE assembling.

`windcheck check-pairs` classifies a list of candidate edges:

```sh
windcheck check-pairs --root patches/ --edges edges.csv --out report/
```

`edges.csv` is either a two-column CSV with an `a,b` header or plain
whitespace/comma-separated pairs of surface directory names under
`--root`. Every edge gets one of three verdicts:

- `transverse_conflict` — the two surfaces interpenetrate; any merge that
  keeps both in full contains a self-intersection.
- `no_transverse_conflict` — the two surfaces do not pass through each
  other. Deliberately not "compatible": the test is silent about
  parametrisation, seams, scale agreement and every failure mode that is
  not interpenetration.
- `not_testable` — the pair cannot be decided (extents that cannot meet, a
  surface contributing no triangle, duplicate geometry under two names, or
  too little shared region to call the absence of contact decided).
  Reported, never silently dropped: a gate that discards what it cannot
  decide reports a fiction.

Each surface is also censused alone, and its own contacts are reported in
`self_a`/`self_b` rather than blamed on the pair — a patch that
self-intersects is a separate finding. A found contact outranks every
applicability heuristic: no sparseness rule can suppress a detected
crossing.

Report-only, like `check`: the command has no flag that writes geometry.

## How it runs

The C++ census loads every surface of a batch from one atlas, scopes
adjacency exclusion to within each surface, and tags each contact with the
two surface ids it joins, so one engine pass per triangulation answers a
whole batch of edges at once. Verdicts do not depend on how the edge list
is batched — the broad phase assigns each triangle pair to exactly one
owner cell, so a pair is tested exactly once wherever the grid origin
falls — and this is tested, not assumed: the suite asserts identical
results from one edge per batch up to every surface in one atlas, and
holds the batched path field-for-field against an independent
stitched-grid implementation.

Measured on the 459-pair example below: 459 edges in ~20 s end to end on a
laptop, including reading every surface off disk.

## The 459-pair example

`results/pairs/` holds a complete, reproducible run:

- `edges.csv` — 459 pairs of individually transverse-clean patches from
  the published `PHercParis4/verified_patches` dataset
  (`https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/verified_patches/`),
  drawn uniformly from the ~27.8M bounding-box-overlapping pairs among
  42,347 downloaded patches (reservoir sample; selection recorded in
  `bench/patch_pairs.py`).
- `pairs.json` — the tool's full output: 10 `transverse_conflict`,
  415 `no_transverse_conflict`, 34 `not_testable`.

To reproduce: download the patch directories named in `edges.csv` from the
URL above (each is `meta.json` + `x.tif`/`y.tif`/`z.tif`), then run the
command with `--root` pointing at them. The run is deterministic.

The 10 conflicts are the composition finding documented in
[PATCH-AUDIT.md](PATCH-AUDIT.md): individually clean patches do not
necessarily compose, and the contacts look like two patches disagreeing
about where the same sheet is by a few voxels (shallow angles at real
penetration), not like folds. The 10-of-459 rate is a property of
box-overlapping random pairs, NOT of any real merge candidate list; do not
extrapolate it.

## Relation to the upstream tool

The same census (single-surface form) is part of volume-cartographer as
the standalone app `vc_tifxyz_selfcross`
([ScrollPrize/villa#1303](https://github.com/ScrollPrize/villa/pull/1303),
merged 4 August 2026), so a surface can be checked right after export with
no external tooling.
One difference is disclosed there and here: upstream loads a surface with
`z <= 0` cells invalidated before the mask, and windcheck's published
counts treat them as valid. Measured across the pinned corpus the gap is
1.49% of contact rows and changes no verdict.

## Limits

The comparison is between coordinate tuples. Nothing verifies that two
surfaces share a volume and voxel scale; an edge list produced by one
workflow satisfies this by construction, one assembled by hand may not.
Contact counts are triangle-pair rows summed over both triangulations,
not distinct crossing events — they rank severity, they do not count
places.
