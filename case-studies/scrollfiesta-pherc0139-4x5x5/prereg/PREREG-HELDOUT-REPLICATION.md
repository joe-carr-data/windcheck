# Pre-registration: held-out replication of the unroll-export finding

Written 2026-08-05, BEFORE any carve, run, export or census of the
held-out window, per house rule. an external review requirement: a second grid is
required for any claim beyond the canonical fixture.

## Design

- HELD-OUT WINDOW (named now): the z-adjacent disjoint 4x5x5 window of
  the SAME sample and volume -- bbox 4864 5376 3072 3712 2560 3200,
  umbilicus (y,x) = (3405, 2878), same pred/raw zarrs, same carve
  script, same canonical flags at every stage. It shares NO voxels with
  the original window (z 4352-4864). The stronger 4x21x21 canonical run
  is explicitly deferred on operational cost; this is a disjoint 4x5x5
  window and will be described exactly as that.
- ONE pre-named fallback, used ONLY if the primary window is not
  traceable (empty/failed carve or zero meshed cubes): the x-adjacent
  window bbox 4352 4864 3072 3712 3200 3840. No further windows; a
  second failure is the recorded outcome.
- NEGATIVE CONTROL (declared, already published): the 84,316
  verified-patch census (windcheck v0.1.0-beta record: 84,311 clean,
  105x suppression vs corpus tiles) stands as the predicted-clean
  verified-output control. It is cited, not rerun.
- RESOLUTION CONTROL (descriptive only): one additional export of the
  ORIGINAL fixture at --tifxyz-du 2 --tifxyz-dv 2, censused the same
  way, to establish raster-resolution dependence of the counts. No
  claim hangs on it.

## Protocol (all parameters frozen; their canonical flags throughout)

1. Carve with their script (anon S3), same flags as the original except
   the bbox above.
2. grid_pipeline --halo 13 --exe build/cube_mesh; per-cube failures
   recorded and accepted as deviations (the original window had 1/100).
3. scroll_whole placement, then --reregister --audit. THEIR quality
   gate result is recorded and reported either way; a FAIL does not
   stop the experiment but is stated prominently with the result.
4. scroll_unroll --steps 12345 + --export-tifxyz + facemap tap
   (byte-neutrality already proven for the tap; not re-proven).
5. Census the export: corpus-frozen parameters, both diagonals.
   Classify rows by the frozen four-corner label rule; cluster events
   with bench/crossing_analyse.py's region-pair definition verbatim.
6. Apply the frozen transform through bench/fiesta_adapter.py; record
   retention, handover verdict, and the official-loader reload
   (merged vc_tifxyz_selfcross) on the adapter output.

## Endpoints

- PRIMARY (replication): >= 1 transverse event in the held-out final
  export. The finding replicates iff the primary holds.
- Secondary, descriptive: rows/events (d0/d1); label-disjoint split;
  their audit gate result; transform retention and 0/0; official
  reload; the du=dv=2 resolution-control counts on the original
  fixture.

## Outcomes (fixed now)

- REPLICATES: >= 1 event. Claim allowed: "the finding replicates on a
  disjoint window of the same sample" -- never more than that; other
  samples remain untested.
- DOES NOT REPLICATE: zero transverse events in the held-out export.
  Reported as such; the canonical-fixture result stands alone and every
  public statement is scoped to that fixture.
- NOT TRACEABLE: both named windows fail to produce a placed grid.
  Reported; no third window.

## Kill criteria

One day of wall clock for the whole chain. No parameter tuning, no
window shopping beyond the one pre-named fallback, no event-rule
changes. Every count reported including zeros.
