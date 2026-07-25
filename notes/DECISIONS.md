# Decision log

Append-only. Every entry records what was decided, why, and what would reverse it.

## 2026-07-25 — target: cross-winding consistency

**Decided:** build a tool that detects sheet-switches from geometry alone, rather
than reporting individual data issues one at a time.

**Why:** a sheet-switch is a *systematic* failure — the same mechanism recurs
across traces and scrolls — so a check that runs over the whole published corpus
is worth more than any single finding, and it keeps working as new segments are
released. It also needs no ground truth, which is the binding constraint here:
there is no published set of known-bad segments to learn from.

**Reverses if:** the invariant turns out not to separate correct traces from
incorrect ones (see the entries below, where three stronger versions of the claim
were tested and dropped).

## 2026-07-25 — invariant chosen: winding index monotonicity

**Decided:** use the 44 contiguous single-winding segments (w052-w095) as a
reference atlas, and require the winding index along a trace to be monotone and
continuous.

**Why:** it is a discrete, certifiable violation rather than a learned score; it
needs no ground-truth labels, no GPU, and no access to the 1.27 TB volume.

**Reverses if:** the wNNN labels turn out to be derived from the trace (task #11),
which would make the check circular.

## 2026-07-25 — dropped: dual-scan differential oracle

**Decided:** do not use the 838/839 tifxyz pair as independent reconstructions.

**Why:** the catalog gives both `bbox_transformed` identical to the last decimal
and `overlap_ratio: 1.0` for both volumes. Same geometry, one coordinate frame.
Checked before relying on it; it was a bonus signal, not load-bearing.

## 2026-07-25 — GATE FAILED: the reference-atlas approach is dropped

**Result:** the pre-registered gate fails. Dropping the atlas-based detector.

**What was built and works:** `engines/atlas_query.cpp` (nearest point-to-surface
winding lookup, uniform-grid accelerated, 12-thread, ~27k queries/s), the tifxyz
reader, and the atlas assembly. None of this is in doubt.

**What the measurements established (all reproducible):**

1. True sheet spacing is **12.3 voxels (~97 um)**, not the ~30 that a
   nearest-*sample* lookup suggests. Point-to-surface interpolation was
   necessary, not a refinement.
2. The `wNNN` labels **are a true radial index** — median distance rises
   monotonically with |dw| (19 -> 25-30 -> 37-51 -> 51-71 -> 74-93) across three
   independent probes, with modal rank order w, w+-1, w+-1, w+-2, w+-2, ...
   This closes task #11 for the purpose at hand: the ordering is global
   structure, so using it to check an individual trace is not circular.
3. Point-wise winding assignment is only ~82% accurate under leave-one-out.
   Aggregating over v fixes most of that, but the known-clean baseline still
   drifts up to 6 windings with steps of 3 — larger than the +-1 signal a real
   sheet-switch produces.
4. **The killer:** the 44 reference windings occupy median radii 2015-2532
   voxels; the 9 auto_grown traces occupy 779-1970. **Radially disjoint.** The
   traces sit 65-417 voxels from the nearest reference surface against a
   12-voxel sheet spacing, and d1 falls monotonically with trace index
   (286, 276, 286, 417, 372, 314, 220, 136, 66) exactly as that geometry
   predicts. The 52/78/93 profile oscillation was noise, not signal.
5. Across all 311 published segments in 12 samples, **PHerc0172 is the only one
   with both labelled windings and auto_grown traces**, and its two families do
   not overlap. The approach cannot be tested on public data anywhere.

**Not a bug.** Every component checked out; the data simply does not contain the
overlap the method needs. Reported as a negative rather than rescued.

**The pivot this suggests (not yet adopted):** consecutive wraps of a spiral
cannot intersect. That **no-crossing invariant needs no reference atlas at
all** — it applies directly to the 61 auto_grown traces across 5 scrolls that
have no labels, which is exactly where the community has no ground truth. And
PHerc0139 (37 windings) + PHerc1667 (20 windings) provide 57 more labelled
sheets as an independent positive control for it. Same engine, same reader.

## 2026-07-25 — self-gap detector: both controls pass, signal found

**The invariant that survived.** Distance from a point on a trace to the nearest
*other part of the same trace*, excluding its own neighbourhood in u. Needs no
scroll axis and no reference atlas, so it is immune to both failures above.

**Why the axis-based version was dropped:** fitting radius-vs-winding on the 44
labelled sheets -- ground truth, known radial order -- gives residual sd 79.8 vx
against a 12-17 vx sheet spacing. A Herculaneum scroll is crushed, not round.
Without that control I would have reported traces 3/4/7/8 gaining 4-14 sheets
per revolution; it was entirely an artifact of assuming a circle.

**Controls.**
- Calibration: 1/2/3-sheet gaps measure 17.0 / 34.1 / 51.2 vx over 126 surface
  pairs -- ratios 1.00 / 2.00 / 3.01, dead linear. Per-POINT discrimination is
  only 77.8%, so nothing is claimed per point; the verdict is per region.
- NULL control: single-wrap labelled sheets have no previous wrap, and the
  detector returns 0.0-0.5% finite gaps on w060/w070/w078/w090. It stays quiet
  where there is nothing.

**Result on the 9 auto_grown traces.** Median gap 19.7-23.1 vx vs the calibrated
17.0, so the traces are mostly correct. The signal is in the tail and it is
spatially clustered:

    trace  flagged(<6vx)  largest blob  top-5 share  u-span
      0        0.9%             78          25%        36    scattered = noise
      1        2.1%            243          16%        87    scattered = noise
      2        3.1%          1,547          65%       240    localised
      7        8.3%          2,555          46%       381    localised
      4       10.7%          2,140          39%       354    localised
      5       14.5%         13,807          81%     1,533    half the trace

Inside the blobs the gap is ~0 (median 0.37-2.31 vx, 22-69% within one voxel),
i.e. the trace is coincident with itself. That rules out the compressed-region
explanation, where sheets touch but stay distinct.

**NOT YET RULED OUT -- the load-bearing open question.** Exclusion is by u-INDEX,
which assumes u advances monotonically with angle along the trace. If an
auto_grown surface is an irregular patch rather than a neat spiral, two points
60 columns apart in u could be physically adjacent on the SAME wrap, and the
detector would read that as self-overlap. The fix is to exclude by geodesic
distance on the surface instead of by u-index, and to verify that the excluded
band really is far in 3D. **No claim should be published until this is closed.**
This is the first thing an independent review should attack.

## 2026-07-25 — BLOCKER CLEARED: u-index exclusion is sound

The worry was that exclusion by |du| < 60 columns could leave same-wrap
neighbours in play, manufacturing the signal. Falsified directly: in self-gap
mode the engine returns the u-column of the nearest triangle, so the test is
just |u_nearest - u_query| for the flagged points.

    trace          0    1    2    3    4    5    6    7    8
    flagged |du| 271  313  359  436  480  532  562  576  618
    % within 120 of the 60-column cutoff: 0.0% everywhere except trace 7 (1.4%)

Nothing is being caught just past the cutoff. Better, those offsets match the
revolution periods measured independently by angular unwrapping:

    u per turn   236  309  365  418  467  507  547  593  616

Agreement across all nine traces by two unrelated routes. The flagged points sit
at exactly one revolution, so they are genuine wrap neighbours, and a gap of ~0
there means the trace came back onto the sheet it had already traced.

Note this also retro-validates the revolution-period estimate from the session
where the RADIUS fit failed: the angle measurement was fine, only the
radius-from-a-fixed-axis model was broken.

Geodesic exclusion (task #24) is no longer load-bearing. Keep it as a robustness
item, not a blocker.

## 2026-07-25 — CT spot-check does NOT confirm the finding (test was ill-designed)

Verified first that the tifxyz->volume mapping is sound: trace points sit at the
77th percentile of nonzero voxels in a 121^3 box, i.e. on dense material.

Then compared flagged (gap<6) against clean (gap 14-24) regions of trace 5:

                     n    CT intensity at trace   bright bands within +-28 vx
    FLAGGED         38          77.0               median 3.0  mean 2.97
    CLEAN           22          90.1               median 3.0  mean 3.23

Statistically identical. **The pre-registered PASS condition is not met.**

The design was wrong, not merely the result. "How many sheets are near this
point" cannot discriminate: both regions sit in normal papyrus at ~17 vx
spacing, so both yield ~3 bands in a +-28 window regardless of whether the
trace is correct. I measured a quantity that is invariant to the hypothesis.

The discriminating test is the MISSED-SHEET signature: if a trace doubled back
onto a wrap it had already traced, some nearby bright band is occupied by NO
part of the trace. At a correct point every neighbouring band is claimed by an
adjacent wrap; at a genuine doubling one is orphaned. The engine already reports
where the trace's other wraps sit, so this is buildable from what exists.

Weak hint, underpowered, not leaned on: flagged points sit on less dense
material (77 vs 90), consistent with the trace being slightly off-surface.

STATUS: the self-gap signal is real and reproducible as a geometric statement
(the trace IS coincident with itself at one revolution). What is NOT established
is that this constitutes a segmentation ERROR. No claim ships until the
missed-sheet test either confirms or refutes it.

## 2026-07-25 — CT confirmation FAILS on a third, properly-powered design

Missed-sheet test, occupancy measured as distance-to-surface via the engine
(sanity check passes: distance at t=0 is 0.00 in both arms, so the trace really
is where we think it is).

                  n    bands   orphans   orphan rate   %>=1 orphan
    FLAGGED      99     2.99     1.788       59.8%        92.9%
    CLEAN        90     2.92     1.567       53.6%        92.2%

    Mann-Whitney one-sided p = 0.0496,  Cohen's d = 0.245

The difference is in the predicted direction but the effect is tiny. Two things
follow, and the second matters more than the first:

1. Underpowered to establish the effect: d=0.245 needs n~250/arm for p<0.01, and
   the per-box neighbourhood filter capped us near 95.
2. **Even if real, d=0.245 is useless as a detector.** Certifying an individual
   region needs separation, not a shifted mean. This is not a threshold away
   from working.

THREE test designs have now failed to confirm: band-count, missed-sheet v1
(broken occupancy), missed-sheet v2 (correct occupancy, powered). The
pre-registered PASS condition -- "spot-checks in the CT volume confirm the trace
is genuinely on the wrong wrap" -- is NOT met.

**Standing conclusion.** The geometric finding is solid and reproducible: four of
nine traces are coincident with themselves at exactly one revolution, over large
contiguous regions, and that survived direct falsification. What is NOT
established, and should not be claimed, is that this constitutes a segmentation
error. It may be normal for these surfaces, an artifact of how they are meshed,
or a real defect the CT cannot resolve at 7.91 um.

Per the pre-registered gate, no claim ships. The honest next move is to ASK
rather than assert: show the community the geometry and ask whether
self-coincidence at one revolution is expected in auto_grown output.

## 2026-07-25 — out-of-sample: the signal is CONCENTRATION, not rate

PHerc0814's multi-wrap auto_grown traces, run at identical parameters (stride 4,
exclude_u 60, threshold 6 vx) to Scroll 5:

                       traces    flag% median    largest blob median / max
    Scroll 5 (PHerc0172)   9          3.3            276 / 7,607
    PHerc0814             10          3.5             32 /   130

Two separable facts, previously conflated:

1. Scattered self-coincidence at ~3% is a PIPELINE-WIDE property. It reproduces
   almost exactly on an independent scroll, so "flagged % is high" means nothing
   on its own. This kills the naive reading of the earlier result.
2. CONCENTRATION is anomalous. Scroll 5 traces 5 (7,607), 7 (1,445) and 4
   (1,110) exceed anything in the out-of-sample set by 10-59x.

**Revised deliverable.** Not "I found errors" -- the CT cannot support that. What
the evidence does support is a cheap, label-free statistic that RANKS published
traces by geometric anomaly, calibrated across two independent scrolls. That is
Open Problem #4 (route only uncertain regions to humans). No ground truth, no
GPU, minutes on a laptop, and it yields a triage order.

**Caveat to fix before anything ships.** One PHerc0814 trace reported |du| = 95
against an exclusion cutoff of 60. Close enough that its flags may be the
same-wrap artifact falsified earlier. |du| must become a per-trace VALIDITY
FILTER (reject a trace whose flagged |du| is not comfortably above the cutoff),
not a number checked by eye.

## 2026-07-25 — source audit: the finding is grounded in THEIR OWN constant

Independent source audit of volume-cartographer + lasagna. Verdict SHOULD-BE-PREVENTED.
I independently verified the load-bearing claims before accepting them.

**Verified in source:**
- `same_surface_th = 2.0` — GrowSurface.cpp:121, overridable GrowSurface.cpp:1429.
- The guard, GrowSurface.cpp:2919-2931: `pointTo(...)` against the whole
  already-grown surface, then `if (dist <= same_surface_th) { best_inliers = -1;
  best_ref_seed = false; }`. Its own debug string is
  "candidate rejected: nearest same-surface location had no valid state".
  So the authors treat a candidate landing within 2.0 vx of already-traced
  surface as a defect and try to stop it.
- `resume_local_row_self_intersection_cell_factor` exists ONLY in
  JsonProfilePresets.hpp:50 — a self-intersection guard with no consumer.
- vc_straighten.cpp:13-15 header: "pair each grid point with the nearest 3D
  point one winding over (its self-overlap)". NUANCE missed on the first pass:
  "self-overlap" in their vocabulary means the normal neighbouring wrap, so
  self-overlap per se is expected. Only the DISTANCE distinguishes normal
  (~17 vx) from defective (~0). This is why the sub-2.0 fraction is the number
  that matters, not the flagged fraction.

**The number that reframes the project.** Fraction of points within 2.0 vx of
the trace's own neighbouring wrap, i.e. inside the pipeline's own rejection
threshold (stride 3):

    trace      0     1     2     3     4     5     6     7     8
    % <2.0vx  0.36  0.83  1.97  1.71  5.42  6.60  2.13  5.70  2.71
    total: 48,304 points at stride 3 (~435k at full grid resolution)

**Why this is a better claim than anything before it.** It needs no CT
confirmation and no ground truth, because "defect" is defined by THEIR constant,
not by my judgement. The claim becomes: published auto_grown output violates the
same-surface criterion the grower enforces during growth.

**Honest caveat that must ship with it.** The guard runs at growth time on a
low-res grid (pitch 200 vx) via a stochastic search (QuadSurface.cpp:1443-1490,
greedy descent + up to 1000 random restarts, authors' own FIXMEs). The published
tifxyz is the POST-OPTIMISATION output: GrowSurface.cpp:1362 `points = points_out`
rewrites every point, and optimize_surface_mapping re-runs no such check, nor any
injectivity/fold test. So a violation in the final output does not prove the
guard "failed" -- it may have been introduced afterwards. Either way there is no
stage that re-checks the criterion after optimisation. That gap is the finding.

**A decisive attribution test is NOT available on published data:** it
proposed using the `generations` channel (inpainted cells are generation 0,
GrowSurface.cpp:960) to separate "grower missed it" from "inpainter created it".
Published tifxyz ships only x/y/z/meta.json, no generations channel.

## 2026-07-25 — validity filter enforced in code; rate is pipeline-wide, only concentration is not

Consolidated the ad-hoc scripts into `src/windcheck/selfgap.py` with the |du|
validity filter enforced rather than eyeballed: a trace is rejected when its
flagged |du| p10 falls within 2x the exclusion cutoff, since those partners may
be same-wrap neighbours rather than the next wrap.

**It rejected the two highest-scoring traces.** 0814-20251001060526 (11.12%
sub-2.0, |du| p10 = 70 vs cutoff 60) and 0814-20250925222237 (8.18%, p10 = 108).
Both would have been headline findings; both are artifacts. This is exactly why
the filter had to be code and not judgement.

**After filtering, 19 traces across 2 scrolls:**

    sub-2.0vx rate     Scroll 5  0.36-6.60%     PHerc0814  0.32-6.77%
    largest blob       Scroll 5  up to 13,807   PHerc0814  up to 240

The RATE is pipeline-wide and overlapping -- 0814's worst valid trace (6.77%)
now exceeds Scroll 5's worst (6.60%). So "% of points violating the grower's own
2.0 vx criterion" is NOT distinctive and must not be sold as if it were.

What remains distinctive is CONCENTRATION: Scroll 5 traces 5 / 7 / 4 carry
contiguous blobs of 13,807 / 2,555 / 2,140 cells against a maximum of 240
anywhere in the out-of-sample set, a 58x gap.

Third time this session the first reading dissolved under a control and a
narrower claim survived. The surviving claim is deliberately small:
**scattered violation is normal; concentrated violation is not.**

## 2026-07-25 — adversarial audit: machinery survives, three headline claims die

An independent review rebuilt the engine, re-derived measurements, and read both
codebases. I verified its load-bearing claims myself before accepting them.

### Validated (independent re-derivation)
- Calibration: its 16.75 / 33.76 / 49.87 vx vs my 17.0 / 34.1 / 51.2; ratios
  2.015 / 2.977. Solid.
- Trace 5 reproduced exactly: 6.596% <2vx, 14.541% <6vx, blob 13,807, |du| 532.
- `dist2_point_triangle` correct: max error 3.6e-6 vx vs an independent
  implementation over 250 randomised queries at three cell sizes.
- Shell early-exit bound correct (triangles are inserted into every cell their
  AABB touches).
- Periods re-derived from angular unwrapping with several centres:
  241/312/368/421/470/506/548/594/620 — matches. Trace 5 stable for exclude_u
  from 30 to 240.

### ENGINE BUG (real, mine, now fixed)
`max_dist` did not behave as its own comment claimed: exiting by exhausting
`max_shell` let a finite non-global result escape instead of becoming infinity.
This is exactly the "w078 0.5% finite at median 331 vx" I saw in the null
control and dismissed as far-field noise. Fixed by certifying only results
inside `max_dist - cell`. Null control 0.5% -> 0.079%; trace 5 headline
unchanged at 6.599% / 14.548%, confirming nearby hits were always certified.

### Claims that DIED
1. **"58x concentration anomaly" — dead, it mixed thresholds.** The rate came
   from the 2 vx (grower) threshold and the blob from my 6 vx threshold.
   Normalised by analysed points, the anomaly is single-digit: ~1.8x at 6 vx and
   at 2 vx it can even reverse. Raw 57.5x reproduces but is meaningless
   unnormalised. Must report blob fraction, physical area, and both 4- and
   8-connectivity.
2. **"9 out-of-sample traces" — dead. VERIFIED MYSELF.** All nine Scroll 5
   traces share the artifact id `auto_grown_20251115002740308`; they are
   partitions `_0`.._8` of ONE run. Effective n is 1 + 8 PHerc0814 segments, not
   19. This was visible in filenames I had already printed and read past.
3. **"Published output violates the grower's own criterion" — overstated.
   VERIFIED MYSELF.** The catalog carries 793 `s3://philodemos/...` private
   origins; Scroll 5 geometry passed through private artifacts, flattening,
   orientation, OBJ conversion and tifxyz reconstruction. Nothing establishes
   these meshes came from the audited GrowSurface.cpp revision or parameters.
   Also `same_surface_th` is an acceptance heuristic across several proximity
   semantics, not a documented final-mesh invariant. And the earlier reading's
   "2 vx threshold on a 200 vx grid" was wrong: `pointTo_` interpolates
   continuously and refines to ~0.1% of a cell.

   DEFENSIBLE: published meshes contain concentrated nonlocal
   self-near-coincidence, and current VC code runs no final self-contact check.
   NOT DEFENSIBLE: they violate the criterion of the grower that produced them.
