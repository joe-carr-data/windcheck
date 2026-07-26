# Detecting sheet switches: published scroll meshes pass through themselves

**windcheck — a deterministic self-intersection validator, and an audit of 179
published segments across five Herculaneum scrolls.**

July 2026 progress submission. Everything here is reproducible from this
repository; every number names the file that produces it, and `docs/REPRODUCE.md`
gives the commands, the measured runtimes and the expected output. The whole
audit runs in **91 seconds**.

---

## 1. The result

A traced scroll surface cannot pass through itself. Any place where a published
`tifxyz` surface does is a defect in the representation, whatever the cause, and
no threshold is needed to say so.

Testing all 179 published segments of Scroll 1, Scroll 5, PHerc0139, PHerc0814
and PHerc1667 — 228 million triangles:

```
                            segments   with crossings    triangles    transverse pairs
  Scroll 1  (PHercParis4)       55           48            140.1 M          633,415
  Scroll 5  (PHerc0172)         53           52             56.9 M          822,072
  PHerc1667                     20           20             19.7 M          648,536
  PHerc0139                     38           29              8.8 M           17,385
  PHerc0814                     13           11              2.6 M           96,856
  ---------------------------------------------------------------------------------
  TOTAL                        179          160            228.1 M        2,218,264
```

**160 of 179 published segments contain at least one place where the traced sheet
passes through itself.**

Those crossings are not one phenomenon. Expressed as a fraction of each
segment's own revolution period, they fall into three bands with a large empty
gap between the first two, and each band is explained by how far around the
scroll its segments reach (§4).

The tool emits a per-trace certificate and an overlay that loads directly into
volume-cartographer's existing point-collection widget, so a result is something
you open rather than something you read.

---

## 1b. Which open problem this addresses

Vesuvius Challenge lists **sheet switching** among its current bottlenecks:

> *"mergers, where two nearby sheets are joined by mistake; holes, where the
> predicted surface disappears; sheet switches"*
>
> *"Meshes can jump from one wrap to another"*
>
> *"a small local error can send a traced mesh onto the wrong wrap entirely,
> with no easy way to recover"*
>
> *"label quality is now one of the main unwrapping bottlenecks"*
>
> *"automatic growth still needs human inspection and correction"*

A trace that switches sheets and later returns visits the same physical place at
two points far apart along its own parameter. That is exactly the wrap-scale band
in §4: two parts of a trace **1.8 revolutions or more apart in parameter, meeting
in 3D**. Seventeen of the 179 segments show it.

This is a detector for a specific, checkable *consequence* of sheet switching,
and the scope of that claim matters:

- A trace that switches wraps and never comes back need not self-intersect, so
  this is not a complete sheet-switch detector.
- A self-intersection is not by itself proof that a switch occurred — see §9,
  where three separate attempts to establish that failed.

What it does give is an automated, deterministic, corpus-scale answer to
"which published segments contain a place where the traced sheet meets itself a
full wrap or more away, and where exactly" — for meshes that today rely on human
inspection to find such things. On the five scrolls tested it runs in 91 seconds.

---

## 2. Why this measurement and not proximity

An earlier version of this tool measured **proximity** — how close a trace comes
to another part of itself. Two volume-cartographer maintainers pointed out,
independently and correctly, that this is undecidable:

> *"some wraps are in reality very tightly packed in places (~5um)"*
>
> *"part of it is just the optimization pushing it that way and part of it is the
> 20vx spacing in the quadmesh where in tightly packed regions with some 'wobble'
> … you can bilinearly interpolate to even closer positions"*

Both are right, and together they mean a small distance can indicate a defect,
genuinely tight packing, or an artifact of interpolating a coarse quad mesh. No
threshold separates them.

**Transverse self-intersection is different in kind.** An embedded surface cannot
pass through itself at any packing density, so the tight-packing objection does
not apply — and that is demonstrated rather than asserted. `tests/test_selfcross.py`
contains folds at 3.0 voxels and at **0.5 voxels** of separation, closer than the
~5 µm described above, and both report nothing while a planted crossing is found.

That is the whole reason for the change. It replaces an argument with a
measurement.

---

## 3. What is in this repository

| component | file |
|---|---|
| the predicate | `engines/selfcross.cpp` — Möller interval overlap, classifying TRANSVERSE / COPLANAR / GRAZING with penetration margin and crossing angle |
| corpus census | `bench/crossing_census.py` — both quad diagonals, records the exact mesh it measured |
| events and regions | `bench/crossing_analyse.py` |
| revolution period | `bench/revolution_diag.py`, `bench/period_cross_check.py` |
| four-corpus comparison | `bench/revolution_summary.py` |
| **external control** | `bench/winding_control.py` — period against published winding counts |
| triangulation invariance | `bench/triangulability.py` |
| mesh quality | `bench/crossing_quality.py`, `bench/twist_threshold.py` |
| library agreement | `bench/validate_fcl.py` |
| clustering sensitivity | `bench/clustering_sensitivity.py` |
| consumable output | `src/windcheck/certificate.py`, `bench/make_certificates.py` |
| figure | `bench/figure_strata.py` |
| tests | `tests/test_selfcross.py`, `tests/test_geometry.py` — 21 tests, no data needed |

No GPU. No labels, no ground truth, no model. The core measurement needs no
volume data at all — only the published surface meshes.

---

## 4. The three bands

Millimetres do not compare between scrolls: voxel size differs, scrolls differ in
diameter, and a crushed scroll has no single circumference. So separation is
expressed as a fraction of each segment's **own** revolution period, measured
from its own geometry. No axis is fitted; axis fitting failed its own positive
control on this data (residual sd 79.8 vx against 12–17 vx sheet spacing — these
scrolls are crushed, not round).

Sorting every segment by its maximum separation and splitting wherever
consecutive values differ by 1.4× **and** by at least 0.05 revolutions:

```
                                              covering span
   n     separation, revolutions              of the segment
  ----  ---------------------------          ----------------
   54    0.002  -  0.115                      0.65  -  1.42       local
                            <-- gap 9.56x
   35    0.899  -  1.234                      1.03  -  5.55       one revolution
                            <-- gap 1.47x
   14    1.810  -  5.044                      2.88  -  9.08       wrap-scale
```

**The largest gap is 9.56×, and the next largest anywhere in the data is 1.47×.**
No threshold is chosen and defended; the dominant feature of the distribution is
a gap, and that gap is the criterion.

The third column is what gives the bands meaning, and it was not placed by hand:

- **Local** overlaps occur at every covering span.
- **One-revolution** overlaps occur *only* in segments covering 1.03 revolutions
  or more — never below, where they are geometrically impossible. These are a
  segment's own two ends covering the same angular sector: **correct geometry,
  not a defect.** They appear consistently across three scrolls.
- **Wrap-scale** overlaps occur only in traces covering 2.88 revolutions or more.

So the useful output is not a verdict but a distance: how far apart, in turns of
the scroll, the two overlapping parts of a trace lie. That distinguishes a pinch
from a segment closing on itself from a trace that has jumped its wrap.

`bench/figure_strata.py` plots all of it, with the period appearing on one axis
only, so an error in the period estimate cannot manufacture the vertical
structure.

---

## 5. The period is checked against published winding counts

Everything in §4 divides by an estimated revolution period, so that estimate
needed an external check rather than an internal one.

Scroll 1 supplies it. Its segments are named by the range of windings they cover
— `w073-076`, `w116-117`, `w010-027` — a count produced by the publishers' own
tracing, not by us and not from the geometry we measure. `bench/winding_control.py`
compares the two on the 33 segments that carry both a range and a mesh on the
annotated volume:

```
   31 of 33 agree
      correlation declared vs measured    r = 0.9999
      mean absolute error                 0.033 windings
      declared range                      2 to 18 windings

      w116-117   declared  2   measured  2.022
      w089-091   declared  3   measured  2.984
      w038-045   declared  8   measured  7.923
      w010-027   declared 18   measured 18.033
```

Two disagree, and loudly — measured 1.010 and 0.257 against declared 18 and 10.
No cut could mistake 0.257 for 10. Two candidate explanations were tested and
**both rejected**: it is not proximity to the scroll core, because the same
winding ranges are published twice and the other trace of each measures
correctly; and it is not self-intersection density, because the densest segment
in the whole set passes. The cause is open and is not guessed at here.

---

## 6. Three further invariances

**Invariant to sampling density.** PHerc0814 publishes the same traces on three
volumes. Taking one identical physical band at 9.362 / 2.399 / 1.129 µm:

```
                 triangles     max span
   9.362 µm        118,914     215.9 mm
   2.399 µm      1,840,264     215.1 mm
   1.129 µm      8,400,324     215.4 mm
```

0.3% spread across a 70× change in triangle count. What this does *not* show is
stated too: these are one segment mapped into three registered volumes, not three
independent reconstructions, so persistence is expected if they share source
geometry.

**Invariant to triangulation.** A quad grid splits into triangles two ways per
quad. Testing all four combinations for each intersecting pair, **~90% cross
under all four**, and the argument for those needs no assumption about how a
global assignment couples pairs: any global choice gives some specific pair of
diagonals, and that combination is one of the four. Scope travels with the claim
— vertices fixed, connectivity fixed, no added vertices.

**A maintainer's hypothesis, supported.** It was suggested that quad-interior
"wobble" lets bilinear interpolation reach closer positions than the samples
support. An earlier check of ours used a discrete Laplacian and found nothing —
but a Laplacian is *blind to exactly this*, since a saddle `z = uv` has zero
Laplacian and nonzero twist. Measuring the bilinear cross-term
`|p00 − p10 − p01 + p11|` instead, local crossings show **9.40×** twist
enrichment (52 of 52 segments above 2×) and wrap-scale ones 1.21×. That is
association, not causation, and a related prediction of ours failed: we expected
a twist threshold above which a quad cannot be triangulated cleanly, and there is
none.

---

## 7. Verification

**FCL**, an established C++ collision library, in both directions:

```
   POSITIVES  (we say transverse)     FCL agrees  250/250  = 100.0%
   NEGATIVES  (nearby, unreported)    FCL agrees  249/250  =  99.6%
```

The single disagreement is explained and both tools are right: quads
(547, 2256) and (547, 2258) share a vertex **exactly**, so FCL reports contact
and we report no penetration. They touch without crossing. It was predicted in
the script before it ran, and it surfaced a separate mesh degeneracy worth
reporting on its own.

**Determinism.** Counts are invariant across thread counts and broad-phase cell
sizes — 10,907 transverse and 1,989 grazing across nine configurations, pinned by
two regression tests. This matters because it was not always true, and fixing it
found two real bugs.

**An independent implementation** agreed on 6,000 randomised float32 triangle
pairs with zero disagreements away from grazing cases.

---

## 8. What this is for

`bench/make_certificates.py` writes two things per trace.

A **certificate** in physical and revolution units, carrying its own caveats so
they cannot be separated from the number, and enough provenance to reproduce it.
It discriminates without a human reading anything.

An **overlay** in volume-cartographer's own `PointCollections` JSON schema (keys
read from `core/src/PointCollections.cpp`), so it opens in the existing widget
with no transform — coordinates are already volume voxels. One representative
point per crossing *event*, not per pair: 379 points for the worst Scroll 5 trace
rather than 117,445, which is the difference between an overlay someone uses and
one they close.

**The gap this fills.** `volume-cartographer` carries no self-contact,
injectivity or fold check on the mesh representation, and its one relevant
parameter, `resume_local_row_self_intersection_cell_factor`, has no consumer in
that tree. Related work exists elsewhere in the same monorepo and is not
superseded here: a vertex-proximity penalty in the fibre-registration optimiser,
a UV-injectivity check in the unwrap pipeline that detects the *dual* case (close
in UV, far in 3D), nearest-vertex overlap masks in the neural-tracing datasets,
and an animation-clearance certificate which states that intra-roll geometry is
"intersection-free at source" and is **"deliberately not measured"**
(`scroll-unwrap-pipeline/src/scrollkit/metrics/clearance.py:18`). That stated
assumption is what this measures.

---

## 9. What this does not establish

Stated in the same voice as the results, because a reader with the data will find
it.

- **It does not show that any crossing is a tracing error.** It shows the surface
  overlaps itself. Three separate attempts to connect crossings to sheet
  misassignment were made and none succeeded: CT adjudication (best effect
  Cohen's d = 0.245), branch identification (50%, chance), and per-voxel instance
  labels — where the pairwise design turned out to be circular, because two
  quads that intersect are at the same point and so receive the same label by
  construction, and where annotation covers only ~0.6% of any segment's cells.
- **It does not identify which branch of a doubled trace is wrong.**
- **It does not establish a cause.** Tracing, interpolation, reparametrisation
  and optimisation all remain possible; no evidence here attributes these to any
  pipeline stage.
- **One-revolution overlaps are correct geometry** and are reported as such, not
  counted as defects.
- **It is a deterministic floating-point validator, not exact predicates.** Plane
  tolerance is derived from the operands, `16 × FLT_EPSILON × max|coordinate|`,
  because the input is float32 running to ~1.3e4 voxels where one ULP is already
  ~1e-3 vx.
- **Effective sample is smaller than 179.** Segments within a scroll are
  correlated, and Scroll 5's nine multi-wrap files are partitions of one artifact.
- **Event counts are definition-dependent** (1.65× for local overlaps between two
  defensible clustering rules, 1.04× for wrap-scale, on the analysed subset);
  span is not.

---

## 10. Errors found and corrected along the way

Included rather than hidden, because every one of these moved our own numbers
**down**, which is the direction that suggests they are real.

| bug | symptom | effect |
|---|---|---|
| absolute-mm wrap-scale cut | held on one scroll, inverted on two others | headline rewritten in revolutions |
| one corpus read at the wrong resolution | an 11-segment cluster nobody could explain | the anomaly *was* the bug |
| per-thread deduplication | 11,211 pairs on one thread, 21,213 on twelve | every early count invalid |
| endpoint touch counted as crossing | intervals `[0,1]` and `[1,2]` reported transverse | semantics overstated |
| fixed `EPS = 1e-6` | finer than the float32 input's own resolution | −6.7% of our own headline |
| no AABB early-out | grazing tracked broad-phase cell size | 92–99% of grazing spurious |
| pair-weighted quality statistic | described as per-quad | split understated, 8.4× vs 9.4× |
| thread count of 0 spawned no workers | every query returned "nothing found" | caught before it reached a result |

Two were found by our own determinism checks rather than by review. A second
measure of a corrupted input agrees with the first, so the resolution-mismatch bug
was caught only by checking a result against the *provenance* of its input — the
census now records the exact mesh it measured, and the analysis refuses to run
against a different one.

---

## 11. Reproducing

```sh
uv sync --extra viz
clang++ -O3 -std=c++17 -pthread -o engines/selfcross   engines/selfcross.cpp
clang++ -O3 -std=c++17 -pthread -o engines/atlas_query engines/atlas_query.cpp
uv run pytest -q                                  # 21 tests, no data needed

uv run python -m windcheck.fetch --sample all     # 18.07 GB, hash-pinned
uv run python -m windcheck.fetch --sample all --verify

sh bench/time_pipeline.sh                         # the whole audit, timed
```

`docs/REPRODUCE.md` has every step against the section it produces, the measured
runtime of each, what each output file contains, and what a mismatch means.

All inputs are pinned by SHA-256 across four manifests: 1,724 files, 18.07 GB.
