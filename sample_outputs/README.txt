sample_outputs -- real runs, reproducible from published data
=============================================================

Every file here was produced by the commands shown at the top of it, against
the Vesuvius Challenge open-data bucket. Nothing is hand-edited.

    calibration_PHerc0172.txt   what one, two and three sheets of separation
                                 measure on Scroll 5: 17.05 / 34.14 / 51.22
                                 voxels, ratios 1.000 / 2.003 / 3.004, over
                                 126 surface pairs
    status_PHerc0172.txt        the corpus: 44 contiguous labelled windings
                                 (w052-w095) and 9 auto_grown traces
    selfgap_PHerc0172.txt       the check over all 53 Scroll 5 surfaces
    selfgap_PHerc0814.txt       the same over PHerc0814, a different scroll
    *_report.json               the same tables as JSON

To regenerate:

    uv run python -m windcheck.fetch --sample PHerc0172
    uv run windcheck calibrate data/scroll5_tifxyz
    uv run windcheck selfgap   data/scroll5_tifxyz --json report.json


WHAT TO LOOK AT
---------------

1. The 44 labelled single-winding segments report 0.00% flagged.
   A one-wrap segment has no previous wrap, so a working check must find
   nothing on them. It does. This is the null control and it is in the output,
   not just asserted.

2. Four labelled segments are marked REJECT, not flagged.
   The validity filter caught them: their flagged partners sit at |du| of 60-67
   against an exclusion cutoff of 60, so those are same-wrap neighbours, not a
   real wrap away. The check refuses to report a row it cannot stand behind.

3. On PHerc0814 the filter removes the two highest-scoring traces
   (11.27% and 8.18%), both artifacts by the same test. Without it they would
   be the headline numbers.

4. Scattered flagging at roughly 3% appears on both scrolls.
   It is a property of this data, not a defect signature, and is not
   informative on its own.

5. Concentration is what separates traces. Compare blob% -- the largest
   contiguous flagged region as a fraction of analysed points:

       PHerc0172 auto_grown _5    7.02%
       PHerc0172 auto_grown _7    1.52%
       PHerc0172 auto_grown _0    0.06%
       PHerc0814 highest valid    2.11%

   Trace _5's largest single region spans roughly half its length.

Use blob%, not the raw cell count, for any comparison: raw counts scale with
trace size and with --stride.
