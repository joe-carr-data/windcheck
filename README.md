# Windcheck

**Certified topology auditing and transverse-clean `tifxyz` outputs for a
five-scroll corpus.**

Of 185 pinned trace artifacts, 179 were censusable. All 179 now have
reload-verified `tifxyz` outputs with zero non-adjacent transverse contacts
under both canonical triangulations: 154 were transformed and 25 required no
change. Six triangle-empty or invalid inputs have explicit terminal records.

Two figures qualify that outcome, and they are stated separately because they
measure different things:

- **Represented surface retained.** Unique-geometry-weighted retention of
  original represented surface area is **99.505%** across the 185 priced
  segments. The lowest per-segment retained fraction is **94.678%**
  (`20251001060526-auto_grown_20251001060526760`).
- **Fragmentation.** Five outputs do not meet the connectivity gate: every
  component of the input 99.9%-area core was required to keep at least 90% of
  its area in one descendant, and in five `auto_grown` segments it does not.
  The lowest core value is 0.034. This is material, and it is not a bookkeeping
  artifact.

Every emitted artifact was independently re-censused in a fresh working
directory: 154/154 transformed meshes re-hashed clean, 179/179 re-censused at
0/0 non-adjacent transverse contacts under both canonical triangulations, with
zero disagreements against the recorded census.

Targets [Open Problems](https://scrollprize.org/2026_open_problems) §2 and §3 —
sheet switches and mesh topology repair — on
[Vesuvius Challenge](https://scrollprize.org) open data.

---

## Install

```sh
uv sync
clang++ -O3 -std=c++17 -pthread -o engines/selfcross engines/selfcross.cpp
uv run pytest -q
```

The triangle-intersection kernel is the only compiled component: a
deterministic floating-point triangle-intersection predicate with scale-aware
tolerances, cross-validated against independent implementations — not exact
arithmetic. No GPU, no volume download and no ML model is involved: the
analysis reads only the `tifxyz` surface itself.

Fetch a segment to work on:

```sh
uv run python -m windcheck.fetch --help
```

---

## Quick start

The tool has two commands. **`check` never modifies anything**, and it is the
default posture. Transformation only happens when you ask for it by name.

### 1. Audit a segment (report-only)

```sh
uv run windcheck check data/scroll5_tifxyz/20251205115859-w094_20251205115859448_flatboi \
    --volume 20241024131839 --out out/check/w094
```

This runs the both-diagonal non-adjacent transverse-contact census, prints a
verdict, writes a machine-readable certificate, and emits a viewer-loadable
`PointCollection` of the crossing sites so you can jump straight to them in
VC3D. On this segment: **not clean**, 4 transverse contacts under diagonal 0 and
7 under diagonal 1, in 0.2 s.

### 2. Transform a segment (explicit)

Two worked examples, one per operator.

```sh
# w094 — resolved by bounded displacement
uv run windcheck transform data/scroll5_tifxyz/20251205115859-w094_20251205115859448_flatboi \
    --volume 20241024131839 --out out/transform/w094

# 20231005123336 — resolved by certified excision
uv run windcheck transform data/scroll1_tifxyz/20231005123336 \
    --out out/transform/20231005123336
```

Both finish in under three seconds on a laptop. Across the whole corpus the
median segment takes 36.8 s, the 90th percentile 140 s and the slowest 211 s,
against a 600 s policy limit — a full-corpus pass is an afternoon, not a cluster
job.

`transform` exits non-zero unless the *reloaded* output censuses clean, which
makes it usable as a post-export regression test or a CI gate.

### 3. Test candidate merge edges (report-only)

```sh
uv run windcheck check-pairs --root patches/ --edges edges.csv --out out/pairs
```

A per-surface verdict says nothing about assembly: two individually clean
surfaces can still pass through each other, and any merge that keeps both in
full then contains a self-intersection. `check-pairs` classifies each candidate
edge as `transverse_conflict`, `no_transverse_conflict` (deliberately not
"compatible") or `not_testable` — undecidable edges are reported, never
dropped. One batched census answers a whole edge list: the 459-pair example in
[`results/pairs/`](results/pairs) runs in ~20 s on a laptop.
[`docs/CHECK-PAIRS.md`](docs/CHECK-PAIRS.md) has the verdict semantics, the
reproducible example and the limits.

---

## Why both operators are necessary

The obvious repair is to nudge the mesh: move vertices until the surface stops
touching itself. That works, and it is tried first, because it removes no
material. But it is not universal.

Across the corpus, **1,477 crossing events carry certificates of rigid
infeasibility with explicit witnesses** — no bounded displacement under the
stated search and budget separates those contacts. For that population the only
certified route to a clean surface is excision. Excision finds a feasible
low-area mask. Minimum area is claimed only where the certificate records
proven optimality.

So the declared pipeline is **bounded displacement first, certified excision for
the residual**. 103 of the 185 segments enter excision on a
displacement-repaired base, each hash-verified against its repair certificate;
82 enter on the original published mesh.

The excision operator searches for low removed area subject to a
feasibility constraint, and records per component whether the area it found
is proven optimal. It does not optimise connectivity — which is exactly why
the five fragmented outputs above are reported rather than hidden.

---

## Certification

Nothing here is called clean because the optimiser thought it was finished. A
segment is clean only when the emitted `tifxyz` has been written to disk, read
back from disk, and re-censused from scratch at 0/0 under both canonical
triangulations.

Each certificate records:

- input and output mesh identity as a **content manifest**, not a hash of the
  coordinate planes alone: one explicit row per file a reader consumes
  semantically — `x.tif`, `y.tif`, `z.tif`, `mask.tif`/`mask.png`,
  `meta.json` — carrying path, byte size and SHA-256, plus a single directory
  digest over a frozen serialisation of those rows. A declared file that is
  absent is written in as an absent row rather than skipped, so "no mask" and
  "a mask I did not look at" hash differently. The mask is what decides which
  triangles exist, so it is inside the digest;
- grid shape, valid-vertex and retained-quad counts, voxel scale;
- the census engine's own binary and source hashes, and every census parameter;
- per-diagonal contact counts before and after, crossing events, maximum
  penetration;
- the frozen policy version and hash, and the per-component solve method
  (`exact_optimal`, `lp_improved`, `greedy_feasible`) — with **no minimum-area
  claim made unless every component is proven optimal**;
- whether the cut base was the original published mesh or a
  displacement-repaired one, decided by **semantic manifest equality** against
  the original mesh and never by comparing paths, since paths change between a
  downloaded archive and a fresh working directory while the bytes do not;
- removed-area accounting under both denominators, and the identity residual;
- explicit caveats.

Geometric certification and optimisation optimality are deliberately kept
apart. A certificate can say "this surface is certified transverse-clean
under the stated validator" and "this cut is not proven minimal" at the same
time, and often does.

### Verifying provenance without the private history

Nothing published cites a repository commit. A commit sha is not evidence to
someone holding a release: they have no repository to resolve it against. So
certificates and `out/release/index.json` record provenance that is
recomputable from the published files alone —

- the **code version** and the **frozen policy version and hash**;
- a **source-tree digest**: a SHA-256 over the canonical file manifest of
  every published file that can change a result (`pyproject.toml`, `uv.lock`,
  `src/windcheck/*.py`, `engines/*.cpp`, `engines/*.h`, `bench/*.py`),
  serialised by the same frozen rule the mesh manifests use;
- the engine binary and engine source SHA-256, and the lockfile SHA-256.

Recompute the digest and compare it against the one in the index:

```sh
uv run python -m windcheck.provenance
uv run python -m windcheck.provenance --verify <source_tree_digest from index.json>
uv run python -m windcheck.provenance --manifest       # every file that feeds it
uv run python -m windcheck.provenance --serialisation  # the exact bytes hashed
```

The serialisation is frozen and stated in `src/windcheck/manifest.py`: rows
sorted by path, each row `path\0size\0sha256` (or `path\0absent\0absent`),
joined by newlines with a trailing newline, hashed as UTF-8. A reader can
reimplement it in a dozen lines and check every digest in the release without
this repository's history.

---

## The corpus

185 pinned trace artifacts across five scrolls:

| corpus | artifacts |
|---|---|
| Scroll 1 | 55 |
| Scroll 5 | 53 |
| PHerc0139 | 38 |
| PHerc1667 | 20 |
| PHerc0814 | 19 |

Per-segment results — input hash to output hash, disposition, certificate path,
area retention under both denominators, fragmentation metrics and runtime — are
in [`docs/CORPUS.md`](docs/CORPUS.md), with the machine-readable index at
`out/release/index.json`. The six non-censusable inputs are listed there
explicitly, by name and by reason.

The **expansion corpus** extends the same treatment to the 99 further
traces the full census found across seven more samples: every trace
dispositioned under the identical frozen policy, every transformed output
independently re-censused clean under both triangulations, and the four
census refusals decline-confirmed. Records in
[`docs/CORPUS-EXPANSION.md`](docs/CORPUS-EXPANSION.md), machine-readable
index at `out/release/expansion_index.json`, meshes and certificates in
the `v0.2.0` release archives. Combined over both inventories — 274
censusable of 284 indexed artifacts — the area-weighted retained fraction
is **0.994566**, recomputed from per-segment areas.

---

## Reproducing

## Every published trace

Twelve samples in the open-data bucket publish surface traces as `tifxyz`.
**All 278 of them have been censused** — the release covers the first 179,
and [`docs/FULL-CORPUS.md`](docs/FULL-CORPUS.md) adds the remaining 99
across seven samples that use a second published layout. As of the
`v0.2.0` release those 99 are no longer census-only: each has a terminal
disposition, and every censusable one has a verified transverse-clean
output — see [`docs/CORPUS-EXPANSION.md`](docs/CORPUS-EXPANSION.md),
including the exact reconciliation of the 278-census and 284-index
denominators.

That page also reports what the corpus says about size: self-intersection
probability follows a one-parameter independent-cell model at
**q = 7.2 × 10⁻⁶ per valid cell**, within about ten points across five
orders of magnitude. A large surface is likely to fold through itself
mostly because it is large.

## The topology transaction

```sh
uv run windcheck transaction candidate.tifxyz --out final.tifxyz
```

One command stages a supported `tifxyz`, transforms it only if it is
not already transverse-clean (frozen policy, unchanged), verifies the
result including a from-disk reload census under both triangulations,
and commits atomically -- the input is never touched, and nothing is
promoted unless every transformation and verification gate succeeds.
The authoritative transaction report is committed inside the output
(`windcheck_transaction/`); `--report` writes a post-commit copy.
`--adapter scrollfiesta` preserves that pipeline's sidecars; unknown
sidecars are refused, never guessed. `--official-validator` additionally
requires volume-cartographer's merged `vc_tifxyz_selfcross` to certify
the staged result (a bare name resolves through PATH once; the resolved
binary is invoked and hashed). Stable exit codes make it a CI gate.
Run it from a source checkout: the command drives analysis scripts under
`bench/` that are not part of the installed package.

## Case studies

[`case-studies/scrollfiesta-pherc0139-4x5x5`](case-studies/scrollfiesta-pherc0139-4x5x5)
— a pre-registered, replicated measurement of transverse
self-intersection in ScrollFiesta's tifxyz strip export (July 2026 Best
of Month), with a certified opt-in cleanup transaction verified by
volume-cartographer's own merged validator: 327 clustered events on the
canonical fixture, 273–278 on a held-out disjoint window, removed at
99.44% / 99.69% surface retention.

## Auditing someone else's data

The same census was run over the whole `verified_patches` dataset for
PHercParis4, published elsewhere on 30 July 2026: **84,316 patches,
96.7 million valid cells, 40 minutes, no GPU**. 84,311 are transverse-clean
and five self-intersect, against 527 expected at the rate published traces
of the same scroll self-intersect when cut to the same size.

All five share one cause, and one cell of boundary erosion clears every one
of them. [`docs/PATCH-AUDIT.md`](docs/PATCH-AUDIT.md) has the method, the
baseline that makes the number mean something, and the limits. The dataset's
publisher reviewed the finding and removed the five flagged patches.

## Upstream

The single-surface census is now part of volume-cartographer as the
standalone app `vc_tifxyz_selfcross`
([ScrollPrize/villa#1303](https://github.com/ScrollPrize/villa/pull/1303),
merged 4 August 2026), so a surface can be checked right after export with
no external tooling —
report-only, both triangulations, a byte-reproducible JSON report, crossing
sites as a VC3D-loadable point collection, and an exit-code gate for
scripts. It reproduces this repository's census exactly on the published
trace it was verified against (contact identities, penetrations and
angles, not just counts). windcheck is also listed among the
[community projects](https://scrollprize.org/community_projects) on
scrollprize.org.

## The audit record

Before any mesh was transformed, this tool audited the published corpus and
reported what it found. That record is preserved as filed, in
[`docs/submission.md`](docs/submission.md), with the precomputed per-segment
certificates and viewer overlays in [`results/`](results). It measures the
**published** meshes; everything above measures the meshes this tool emits.

[`docs/REPRODUCE.md`](docs/REPRODUCE.md) gives a one-command spot reproduction of
a single segment, the corpus verification command, expected outputs and
runtimes, and how to re-census an emitted artifact independently of this tool.

---

## Limitations

Stated plainly, because weaker versions of several of these were tested and
died.

- **No physical branch identity.** A transverse contact is a property of the
  mesh. Nothing here establishes which physical sheet a surface belongs to, or
  that a flagged region corresponds to a segmentation error.
- **No claim that text improves.** Removing self-intersections changes the
  geometry a downstream renderer sees. Whether that helps ink detection or
  legibility is untested and is not claimed.
- **Nothing is claimed about coplanar or grazing contacts.** The census
  certifies non-adjacent *transverse* contacts under two canonical
  triangulations. Coplanar overlap and grazing contact lie outside what is
  certified.
- **Material fragmentation in five `auto_grown` outputs.** See the figure at the
  top. The search targets area, not connectivity, so this is a real
  property of those outputs and not a reporting choice.
- **Six inputs were never audited.** They are triangle-empty or otherwise
  invalid, and carry explicit terminal records rather than a verdict.
- **"Represented surface retained" is a geometric measure**, not a statement
  about how much papyrus or text survives.
- The corpus is 185 pinned artifacts. It is not the whole open dataset, and
  nothing is claimed outside it.

---

## Layout

```
engines/selfcross.cpp     deterministic float triangle-contact predicate
src/windcheck/
  cli.py                  the two user commands
  manifest.py             mesh content manifests, digests, base-kind rule
  provenance.py           release provenance a public reader can recompute
  pipeline.py             check / transform, end to end
  tifxyz.py               reader and writer, written from the format spec
  excise.py               certified excision under the frozen policy
  repair.py               bounded displacement
  intrinsic.py            quad retention and area accounting
  certificate.py          certificates and viewer PointCollections
  catalog.py, fetch.py    open-data catalog, download, SHA-256 manifest
bench/                    benchmark and corpus drivers (users never need these)
docs/CORPUS.md            per-segment corpus results
docs/REPRODUCE.md         reproduction instructions
docs/submission.md        the July 2026 audit record, as filed
docs/PATCH-AUDIT.md       independent audit of 84,316 published patches
docs/CHECK-PAIRS.md       pairwise conflict test for candidate merge edges
docs/FULL-CORPUS.md       every published trace, and the size model
docs/HISTORY.md           what is production, supporting, retired
results/                  precomputed audit results, all 179 segments
results/patches/          the patch audit summary and its five findings
results/corpus/           the full-corpus census records
```

`tifxyz.py` deliberately does not use the upstream `vesuvius` package: this tool
audits that pipeline, so it must be able to disagree with it.

## Licence

MIT.
