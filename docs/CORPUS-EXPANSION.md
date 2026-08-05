# The expansion corpus: 99 further traces, dispositioned and verified

Additive companion to `docs/CORPUS.md` (the 185 pinned segments). Same frozen transform policy, same census parameters, same independent verification. Nothing in the pinned corpus is restated.

## Reconciliation of the denominators

185 pinned segments (179 censusable) + 99 expansion traces (95 censusable) = 284 indexed artifacts, 274 censusable. The July full-corpus census reported 278 traces = 179 pinned censusable + 99 expansion traces processed; the 4 expansion refusals below the census validity floor and the 6 pinned not-censusable inputs account for 278 != 284 and 274 != 278.

## Retention

- Expansion inventory: **0.985177** area-weighted over 95 censusable traces.
- Worst single trace: `20251217234605-w2_20251217234605189` retains **0.8270** and its 99.9%-core gate is **FAIL** -- reported, not averaged away.
- Combined over both inventories (recomputed from per-segment areas, 274 censusable of 284 indexed): **0.994566**.
- areas are in canonical grid vx^2 of each trace's own volume; voxel sizes differ across scrolls, so the combined figure is a vx^2-weighted summary, not a physical-area claim

## Per-trace records

| trace | sample | disposition | retained | core gate | verified |
|---|---|---|---|---|---|
| `20240711124827-20240618142020` | PHerc0332 | already_clean | 1.000000 |  | yes |
| `20240828190516-20240716140050` | PHerc0332 | already_clean | 1.000000 |  | yes |
| `20250502180708-auto_grown_20250502160708188` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250502180748-auto_grown_20250502160748721` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250502182142-auto_grown_20250502161324419` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250502182456-auto_grown_20250502161202782` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250502183138-auto_grown_20250502162038685` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250502183421-auto_grown_20250502161744358` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250502184201-auto_grown_20250502163549332` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250502184658-auto_grown_20250502163923577` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250502184845-auto_grown_20250502164121265` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250502185519-auto_grown_20250502164303733` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250502205333-auto_grown_20250502181030065` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250510172639` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250510172804` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250511003658-tifxyz` | PHerc0343P | already_clean | 1.000000 |  | yes |
| `20250511200236` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250609204953-z_dbg_gen_00612` | PHerc0500P2 | transformed | 0.999914 | PASS | yes |
| `20250611132918-z_dbg_gen_00302_inp_hr` | PHerc0500P2 | transformed | 0.999993 | PASS | yes |
| `20250611133823-z_dbg_gen_00130_inp_hr` | PHerc0500P2 | transformed | 0.999966 | PASS | yes |
| `20250611134135-z_dbg_gen_00331_inp_hr` | PHerc0500P2 | transformed | 0.999991 | PASS | yes |
| `20250611153340-z_dbg_gen_00275` | PHerc0500P2 | transformed | 0.999774 | PASS | yes |
| `20250611153346-z_dbg_gen_00286_inp_hr` | PHerc0500P2 | transformed | 0.999968 | PASS | yes |
| `20250611153347-z_dbg_gen_00105_inp_hr` | PHerc0500P2 | transformed | 0.999963 | PASS | yes |
| `20250611154256-z_dbg_gen_00354_inp_hr` | PHerc0500P2 | transformed | 0.999991 | PASS | yes |
| `20250611161023-z_dbg_gen_00205_inp_hr` | PHerc0500P2 | transformed | 0.999946 | PASS | yes |
| `20250611161415-z_dbg_gen_00411_inp_hr` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250611170257-z_dbg_gen_00308_inp_hr` | PHerc0500P2 | transformed | 0.999544 | PASS | yes |
| `20250611171318-z_dbg_gen_00440_inp_hr` | PHerc0500P2 | transformed | 0.999919 | PASS | yes |
| `20250611171745-z_dbg_gen_00225_inp_hr` | PHerc0500P2 | transformed | 0.999993 | PASS | yes |
| `20250611172135-z_dbg_gen_00346_inp_hr` | PHerc0500P2 | transformed | 0.998954 | PASS | yes |
| `20250611173014-z_dbg_gen_00506_inp_hr` | PHerc0500P2 | transformed | 0.999936 | PASS | yes |
| `20250628074500-500P2_front` | PHerc0500P2 | transformed | 0.999980 | PASS | yes |
| `20250702235910-auto_grown_20250702235910292` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250703025628-auto_grown_20250703025628283` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20250703034159-auto_grown_20250703034159599` | PHerc1447 | transformed | 0.999960 | PASS | yes |
| `20250716055227-z_dbg_gen_00260_inp_hr` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250716055228-z_dbg_gen_00333_inp_hr` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250716055229-z_dbg_gen_00371_inp_hr` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250716055230-z_dbg_gen_00377_inp_hr` | PHerc0500P2 | transformed | 0.999904 | PASS | yes |
| `20250716055231-z_dbg_gen_00473_inp_hr` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250716055232-z_dbg_gen_00505_inp_hr` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250716055233-z_dbg_gen_00565_inp_hr` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250716055234-z_dbg_gen_00709_inp_hr` | PHerc0500P2 | transformed | 0.999995 | PASS | yes |
| `20250716055235-z_dbg_gen_00779_inp_hr` | PHerc0500P2 | transformed | 0.999981 | PASS | yes |
| `20250716055236-z_dbg_gen_00356_inp_hr` | PHerc0500P2 | transformed | 0.999959 | PASS | yes |
| `20250902170435--5_b2` | PHerc0343P | already_clean | 1.000000 |  | yes |
| `20250902170441--4_b2` | PHerc0343P | already_clean | 1.000000 |  | yes |
| `20250902170447--3_b2` | PHerc0343P | already_clean | 1.000000 |  | yes |
| `20250902171202--2_b2` | PHerc0343P | already_clean | 1.000000 |  | yes |
| `20250902171204--1_b2` | PHerc0343P | already_clean | 1.000000 |  | yes |
| `20250904233748-0_b2` | PHerc0343P | already_clean | 1.000000 |  | yes |
| `20250905172054-1_b2` | PHerc0343P | already_clean | 1.000000 |  | yes |
| `20250919064353-auto_grown_20250918234353791` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919123506-auto_grown_20250919053506178_inp_hr` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919124722-auto_grown_20250919054722721_inp_hr` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919124917-auto_grown_20250919054917419_inp_hr` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919125302-auto_grown_20250919055302061` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919125605-auto_grown_20250919055605407` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919125754-auto_grown_20250919055754487_inp_hr` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919130642-auto_grown_20250919060642061` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919131352-auto_grown_20250919061352722_inp_hr` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919131642-auto_grown_20250919061642215_inp_hr` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919132115-auto_grown_20250919062115592_inp_hr` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919133208-auto_grown_20250919063208154_inp_hr` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919135433-auto_grown_20250919065433578_inp_hr` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919135915-auto_grown_20250919065915092_inp_hr` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919140325-auto_grown_20250919070325767_inp_hr` | PHerc0009B | already_clean | 1.000000 |  | yes |
| `20250919184428-0500P2-wrap01_0919` | PHerc0500P2 | not_censusable | -- |  | yes |
| `20250919184429-0500P2-wrap02_0919` | PHerc0500P2 | not_censusable | -- |  | yes |
| `20250919184430-0500P2-wrap03_0919` | PHerc0500P2 | not_censusable | -- |  | yes |
| `20250919184431-0500P2-wrap04_0919` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250919184432-0500P2-wrap05_0919` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250919184433-0500P2-wrap06_0919` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250919184434-0500P2-wrap07_0919` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250919184435-0500P2-wrap08_0919` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250919184436-0500P2-wrap09_0919` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250920020222-0500P2-wrap10_0919` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250920020223-0500P2-wrap12_0919` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250920020224-0500P2-wrap13_0919` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20250922024644-0500P2-wrap11_0919` | PHerc0500P2 | already_clean | 1.000000 |  | yes |
| `20251028213516-auto_grown_20251028213516907` | PHerc0800 | already_clean | 1.000000 |  | yes |
| `20251028220042-auto_grown_20251028220042762` | PHerc0800 | already_clean | 1.000000 |  | yes |
| `20251028220955-auto_grown_20251028220955262` | PHerc0800 | already_clean | 1.000000 |  | yes |
| `20251028222030-auto_grown_20251028222030940` | PHerc0800 | already_clean | 1.000000 |  | yes |
| `20251028225813-auto_grown_20251028225813045` | PHerc0800 | already_clean | 1.000000 |  | yes |
| `20251029010146-auto_grown_20251029010146642` | PHerc0800 | not_censusable | -- |  | yes |
| `20251105093211-z_dbg_gen_00320` | PHerc1447 | already_clean | 1.000000 |  | yes |
| `20251217233843-w1_20251217233843496` | PHercMANBp | transformed | 0.971616 | PASS | yes |
| `20251217234605-w2_20251217234605189` | PHercMANBp | transformed | 0.827046 | FAIL | yes |
| `20251218010446-w0_20251218010446110` | PHercMANBp | transformed | 0.998685 | PASS | yes |
| `20251218211706-w3_20251218211706689` | PHercMANBp | transformed | 0.993967 | PASS | yes |
| `20251218212128-w4_20251218212128406` | PHercMANBp | transformed | 0.984272 | PASS | yes |
| `20251218212713-w5_20251218212713988` | PHercMANBp | transformed | 0.976401 | PASS | yes |
| `20251219211451-w8_20251219211451561` | PHercMANBp | transformed | 0.996339 | PASS | yes |
| `20251220012955-w7_20251220012955200` | PHercMANBp | transformed | 0.993203 | PASS | yes |
| `20251220015639-w6_20251220015639809` | PHercMANBp | transformed | 0.976774 | PASS | yes |
| `20251222204911-w6_20251222204911527` | PHercMANBp | already_clean | 1.000000 |  | yes |
| `20251222223312-w7_20251222223312675` | PHercMANBp | transformed | 0.998185 | PASS | yes |

Every number above is copied from the trace's own excision certificate or the independent verification record (`verification.json`); the four census refusals carry an explicit decline confirmation (`not_censusable_confirmed.json`).

