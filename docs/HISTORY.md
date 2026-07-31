# What is in this tree, and what it is for

This project reached its current form by discarding most of what it tried. The
discarded work is still in the tree — every file listed under
`src/windcheck/*.py`, `engines/*.cpp` and `bench/*.py` feeds the published
source-tree digest (see `README.md` and `windcheck.provenance`), so removing or
relocating a file would invalidate the provenance recorded in
`out/release/index.json`. Nothing here is deleted for tidiness.

That leaves a fair question for a reader: which of these files produced the
numbers in the release, and which are abandoned experiments left in place for
provenance? This document answers it explicitly.

---

## 1. The current production path

These are the files that produce every claim in the release. If you are
auditing the result, this is the whole surface you need to read.

**Entry points**

| File | Role |
| --- | --- |
| `src/windcheck/cli.py` | The `windcheck` command: `check`, `transform`, `status`. |
| `src/windcheck/pipeline.py` | The two user operations. `check` is report-only; `transform` applies the frozen policy and emits one aggregate `.tifxyz` plus its certificate. |
| `src/windcheck/check.py` | The single-segment `check` implementation and its printed verdict. |

**The cleanliness authority**

| File | Role |
| --- | --- |
| `engines/selfcross.cpp` | The census. A deterministic float triangle-contact predicate over a tifxyz quad mesh, run under **both** quad triangulations (diagonal 0 and diagonal 1) with Chebyshev adjacency exclusion. A mesh is clean when both diagonals report zero non-adjacent transverse contacts. Every "clean" in this release means exactly this and nothing else. |
| `src/windcheck/tifxyz.py` | Reading and writing the tifxyz grid format; the only I/O path for meshes. |
| `src/windcheck/objmesh.py` | OBJ ingest for meshes that do not arrive as tifxyz. |

**The transform**

| File | Role |
| --- | --- |
| `src/windcheck/excise.py` | The frozen scheduling policy (`round28-greedy-first-v1`): bounded displacement first where a displacement-repaired base exists, then one certified excision of every residual transverse contact. Removal is quad-level and hybrid-invalidating; retained coordinates are bit-identical to the input. |
| `src/windcheck/repair.py` | Turns a crossing event into a candidate rigid-clearance displacement. Its output is only ever a *candidate* — acceptance is the census reporting zero crossings under both triangulations. |
| `src/windcheck/clearance.py` | The exact rigid directional clearance primitive (the LP behind `repair.py`). |
| `src/windcheck/intrinsic.py` | Intrinsic separation along the represented surface between the two preimages of a contact, and the retained-quad accounting used for area retention. |
| `engines/geodesic.cpp` | A C++ port of `intrinsic.py` for corpus-scale runs. The Python module remains the reference implementation. |

### Building the engines

The three C++ engines are built from source; no binaries are committed.

```
clang++ -O3 -std=c++17 -pthread -o engines/selfcross   engines/selfcross.cpp
clang++ -O3 -std=c++17 -pthread -o engines/atlas_query engines/atlas_query.cpp
clang++ -O3 -std=c++17 -pthread -o engines/geodesic    engines/geodesic.cpp
```

`engines/selfcross` is required for anything that reports a verdict.
`engines/atlas_query` is required by `windcheck check` (§2).
`engines/geodesic` is only needed for corpus-scale intrinsic-separation runs;
`src/windcheck/intrinsic.py` covers the same ground in Python.

**Reporting and provenance**

| File | Role |
| --- | --- |
| `src/windcheck/classify.py` | The single source of truth for how a result is described. Every entry point that names a verdict goes through it, so one segment cannot receive two different verdicts. |
| `src/windcheck/certificate.py` | Per-segment certificates, plus VC3D-loadable point collections for the contact sites. |
| `src/windcheck/manifest.py` | The frozen mesh-identity manifest and its serialisation rule. |
| `src/windcheck/provenance.py` | The published source-tree digest. Recomputable from the release tree alone; cites no commit sha and no repository history. |
| `src/windcheck/catalog.py`, `fetch.py` | Locating and retrieving published segmentations. |

**Corpus drivers (bench)**

The corpus numbers were produced by these, in this order:

| File | Role |
| --- | --- |
| `bench/corpus_bases.py` | Pins each segment to exactly one base geometry and one volume; emits `out/corpus_bases.json`. |
| `bench/repair_corpus.py`, `bench/excise_corpus.py` | Apply the frozen policy across the corpus. |
| `bench/verify_corpus.py` | The independent re-census: every emitted mesh is reloaded from disk and re-censused from scratch under both diagonals. Nothing is certified on the strength of the run that produced it. |
| `bench/make_certificates.py` | Per-segment certificates for the released corpus. |
| `bench/headline_decision.py` | The area-retention accounting and its two declared denominators. |
| `bench/build_release_index.py` | Builds `out/release/index.json`, the released manifest. |
| `bench/period_cross_check.py` | The independent period cross-check (see §2). |

---

## 2. Supporting utilities still called by current code

These are not part of the contact census and they establish none of the release
claims, but current code imports them, so they are live and must be built.

**`src/windcheck/atlas.py` and `engines/atlas_query.cpp`**

`atlas.py` is the winding reference atlas: the contiguous single-winding
reference segments say where each sheet physically sits, so any 3D point can be
assigned a winding by nearest reference surface, with a confidence from how much
closer that surface is than the nearest surface of a *different* winding.
`engines/atlas_query.cpp` is its kernel — the nearest-surface lookup is a
tens-of-millions-of-triangles workload and does not belong in Python.

It is still on the live path in two places:

- `windcheck check` imports `atlas` (`src/windcheck/check.py`), and its **period
  cross-check invokes the `engines/atlas_query` binary** (`check.py`, the
  `atlas_query` engine path). Build it alongside the census engine:

  ```
  clang++ -O3 -std=c++17 -pthread -o engines/atlas_query engines/atlas_query.cpp
  ```

- `src/windcheck/pipeline.py` imports `atlas` on the same basis.

What the atlas contributes to a result is the `period_status` field
(`agreed` / `disagreed` / `unavailable`) that `classify.py` reports *alongside*,
and deliberately separate from, the crossing verdict. **The release claim — zero
non-adjacent transverse contacts under both canonical triangulations — is
established by `engines/selfcross` alone and does not depend on the atlas.** The
period cross-check is corroboration and is reported as such; where the atlas is
unavailable for a segment, the crossing verdict is unaffected.

---

## 3. Retired exploratory work

Everything in this section is superseded. It remains in the tree only because
the source-tree digest covers it. **None of it supports any claim in this
release.**

### 3.1 Proximity detection — superseded by exact contact census

The project's first year of measurement asked a proximity question: *are these
two parts of the trace within N voxels of each other?* The header of
`engines/selfcross.cpp` records why that line was abandoned: the question is
irreducibly ambiguous. Wraps in a crushed scroll genuinely lie very close
together, so a small distance can never by itself mean a trace is wrong, and
every result was therefore a threshold argument rather than a fact about the
data.

A transverse crossing is different in kind. An embedded surface cannot pass
through itself however tightly it is packed — there is no threshold to argue
about and no interpretation to get wrong. The whole proximity programme was
retired in favour of that binary predicate, and the release measures contacts,
never distances.

Retired with it: the threshold sweeps, separation-distribution studies and
band/strata taxonomies built on proximity. An earlier three-band classification
of segments was withdrawn outright — `classify.py` exists precisely because six
mutually inconsistent thresholds had accumulated across the codebase and had to
be replaced by one verdict.

### 3.2 `src/windcheck/selfgap.py` — retired

Self-gap analysis: for every point on a multi-wrap trace, the distance to the
nearest other part of the *same* trace, excluding a window around the point's
own parameter, compared against the tracer's own growth-rejection threshold.

Retired for two reasons. First, it is a proximity measurement and inherits every
objection in §3.1. Second, the comparison it rested on is not sound as stated:
the tracer's guard runs at growth time on a coarse grid via a stochastic search,
while the published tifxyz is post-optimisation and is rewritten with no
re-check, so "inside the pipeline's own rejection criterion" does not mean "the
pipeline would have rejected it".

Still imported by the retired bench scripts in §3.4 and by
`tests/test_geometry.py`, which pins its behaviour so the file does not rot. No
release number comes from it.

### 3.3 `src/windcheck/ct_check.py` — retired

The attempt to settle whether a self-gap violation corresponds to a *real*
segmentation error, by looking for a missed-sheet signature in the CT volume:
walk the surface normal, find the bright bands (each a papyrus sheet), and ask
whether a flagged region orphans a band that a clean region does not.

The test did not separate the populations. Flagged regions (n=99) orphaned bands
at a rate of 0.60 with 0.93 of samples orphaning at least one band; clean regions
(n=90) gave 0.54 and 0.92 against a near-identical band count (2.99 vs 2.92).
That is the outcome that says the self-gap signal is geometric only and does not
evidence a segmentation error. That result is what retired §3.2,
and having retired its own input, `ct_check.py` had nothing left to test. **No
module imports it.** It is dead code preserved for the digest.

The CT figure and any proximity-era CT findings are likewise not release
evidence and are not cited by `docs/CORPUS.md` or `docs/REPRODUCE.md`.

### 3.4 Retired bench scripts

These reproduce figures and studies from the proximity and self-gap era. They
run, they are covered by the digest, and they establish nothing in the release:

`bench/normal_profile.py`, `bench/normal_analyse.py`, `bench/texture_a1.py`,
`bench/texture_a2.py`, `bench/interp_ceiling.py`, `bench/interp_support.py`,
`bench/branch_id.py`, `bench/orphan_skip.py`, `bench/min_exit_sweep.py`,
`bench/lambda_filtration.py`, `bench/revolution_diag.py`,
`bench/revolution_summary.py`, `bench/twist_threshold.py`,
`bench/winding_control.py`, `bench/label_validation.py`,
`bench/clustering_sensitivity.py`, `bench/difficulty_sweep.py`,
`bench/inject_benchmark.py`.

The planted-defect work (`bench/inject_benchmark.py`, `bench/difficulty_sweep.py`
and `src/windcheck/inject.py`) is a special case: planted meshes were the
pre-registered development harness for the cutter, and `excise.py` still refers
to that contract. Its precision/recall figures belong to the retired
proximity-era detector and are **not** release claims. The release reports a
census outcome on real published segmentations, not a detection score.

---

## 4. A note on internal references

Some module docstrings cite design notes (`notes/…`) that record the review
rounds behind a decision. Those notes are working documents and are not
published. The citations are left in place unaltered because the files they sit
in are covered by the source-tree digest, and editing a docstring would change
the digest and invalidate the provenance recorded in `out/release/index.json`.
Nothing in the published claims depends on reading them: every rule they pin —
the frozen policy version and hash, the manifest serialisation, the census
parameters — is stated in `README.md`, `docs/CORPUS.md`, `docs/REPRODUCE.md` and
in `out/release/index.json` itself.
