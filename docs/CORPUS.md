# Excised corpus

Generated 2026-07-31T14:37:59Z by `bench/build_release_index.py` from the pinned base manifest, the per-segment excision certificates and the independent recensus record. Every number on this page is copied from those files; nothing is hand-entered. The machine readable form of exactly this content is `out/release/index.json` (schema `release_index/v1`).

Throughout, **represented surface retained** means canonical quad area that survives the cut, divided by the canonical quad area of the stated denominator. It is a measure of the mesh and nothing more.

## Summary

| | |
|---|---|
| Pinned trace artifacts | 185 |
| Scrolls | 5 |
| Transformed | 154 |
| Already clean (no change needed) | 25 |
| Not censusable | 6 |
| Displacement-repaired bases | 103 |
| Original published bases | 82 |
| Unique base geometries | 185 |
| Exact duplicate bases | 0 |
| Represented surface retained, operational denominator, area-weighted | 99.5055% over 179 censusable segments |
| Represented surface retained, headline denominator, area-weighted | 99.5055% over 185 pinned artifacts |
| Lowest segment retention (headline) | 94.6785% (20251001060526-auto_grown_20251001060526760) |
| Runtime min / median / p90 / max | 0.07 s / 36.77 s / 139.97 s / 211.26 s |
| Runtime gate | 600 s, 0 over |
| Meshes re-hashed clean in a fresh workdir | 154/154 |
| Independently recensused | 179/179 clean under both canonical triangulations |
| Residual non-adjacent transverse contacts | 0 |
| Disagreements with the recorded census | 0 |
| Timeouts / errors | 0 / 0 |
| 99.9%-area core gate | 174 pass, 5 fail |

### Per scroll

| Scroll | Segments | Transformed | Already clean | Not censusable |
|---|---:|---:|---:|---:|
| PHerc0139 | 38 | 26 | 12 | 0 |
| PHerc0814 | 19 | 11 | 2 | 6 |
| PHerc1667 | 20 | 20 | 0 | 0 |
| Scroll 1 | 55 | 46 | 9 | 0 |
| Scroll 5 | 53 | 51 | 2 | 0 |

## What the columns mean

**Disposition.**
- `transformed` — the input carried transverse contacts; an excision was computed and an output mesh was emitted.
- `already_clean` — the input censused transverse-clean under both canonical triangulations, so no cut was made and no output mesh exists.
- `not_censusable` — the input carries no triangles, so no census, no cut and no cleanliness claim is defined on it.

**Input hash and output hash.** Each tifxyz mesh is three coordinate planes. The input hash column is the sha256 of the `x` plane of the pre-excision base that was actually cut; the output hash column is the sha256 of the `x` plane of the emitted mesh. `index.json` carries all three planes on both sides, plus the hashes of the original published mesh. The `y` and `z` planes are hashed the same way and are checked together, so the `x` plane alone is a convenient abbreviation, not the whole check.

**The two area denominators.** Two retention figures are recorded for every segment, because 103 of the 185 bases are displacement-repaired meshes rather than the published originals.

- **Operational retention** — measured against the PRE-EXCISION BASE, i.e. the mesh that was actually cut. Where that base is a displacement-repaired mesh rather than the published original, this denominator measures what the excision cost relative to the surface it was handed.
- **Headline retention** — measured against the ORIGINAL PUBLISHED mesh: the removed ORIGINAL quad indices are priced using the ORIGINAL published coordinates, and divided by the canonical area of the original mesh's retained quads. When the base IS the original mesh the two denominators coincide by construction.

The two summary figures are also evaluated over different populations, and each row above says which. **Operational retention** runs over the 179 censusable segments: every segment that was cut, plus every already-clean segment, each priced against the base it was actually handed. **Headline retention** runs over all 185 pinned artifacts — 179 of them carry a priced original-coordinate area, and the other 6 are the triangle-empty or invalid inputs, which carry no surface at all and so enter with zero area on both sides of the ratio, where they cannot move the figure.

A segment that was never cut is **not** dropped from the headline denominator. It enters at exactly 1.0, against the full canonical area of its original published mesh. Where the base IS that mesh, the area is on the certificate. Where the base was displacement-repaired, the certificate writes no area block at all — nothing was cut, so there was nothing to price — and the original area is instead recovered by re-reading the original published mesh and dividing its retained quads the same way every cut segment is divided; the result is memoised in `out/headline_original_areas.json`. The headline figure above, and the population it is taken over, are therefore the same ones the frozen decision rule (`bench/headline_decision.py`) reports: truncated to three decimals as that rule prints it — it never rounds a retention figure up — the figure is 99.505% over the same 185 artifacts.

In the per-segment tables below the two columns look identical, and that is a measured result rather than a copy. Of the 179 segments that carry both figures, 82 agree to the last bit — every segment cut from the original published mesh, where the two denominators are the same denominator, plus every already-clean segment, where nothing was removed. Across the rest, all of them displacement-repaired bases that were cut, the largest difference anywhere in the corpus is 1.425e-07 (`20260603042357-5753_-3`): the two figures agree to roughly seven significant figures. Pricing the removed quads on the repaired coordinates rather than the published ones barely moves the area. Both are still reported per segment, because that agreement is an observation about this corpus and not a guarantee.

**R_main.** R_main(C) = area(largest output descendant of input component C) / area(C), reported per INPUT component. It detects a cut that scores near-100% on area while shattering a component into pieces.

- **min R_main** is the minimum over *all* input components, with no size threshold of any kind.
- **area-weighted R_main** weights each input component's R_main by that component's share of input area.
- **core gate** — The 99.9%-AREA CORE is the smallest prefix of the input components, sorted by canonical area, whose cumulative area reaches 99.9% of the input retained area. The core gate passes when every component in that core has R_main >= 0.90.

**Runtime.** Wall-clock seconds for the segment under the frozen policy, against a 600 s gate.

## Inputs that could not be censused

6 of the 185 pinned artifacts carry no triangles at all. No census, no cut and no cleanliness claim is defined on them, so they hold no retention figure and no fragmentation figure. They are listed here rather than dropped, and they still have a certificate.

| Segment | Scroll | Grid | Valid vertices | Retained quads | Reason |
|---|---|---|---:|---:|---|
| `20250928235954-auto_grown_20250928235953722_copy` | PHerc0814 | 30x45 | 863 | 0 | the input carries no triangles, so no census, no cut and no cleanliness claim is defined on it |
| `20260225160055-auto_grown_20260225160055913` | PHerc0814 | 33x33 | 578 | 0 | the input carries no triangles, so no census, no cut and no cleanliness claim is defined on it |
| `20260226110425-auto_grown_20260226110425097_abf` | PHerc0814 | 30x48 | 868 | 0 | the input carries no triangles, so no census, no cut and no cleanliness claim is defined on it |
| `20260226113324-auto_grown_20260226113324624_abf` | PHerc0814 | 30x49 | 883 | 0 | the input carries no triangles, so no census, no cut and no cleanliness claim is defined on it |
| `20260226121302-auto_grown_20260226121302716_abf` | PHerc0814 | 30x49 | 870 | 0 | the input carries no triangles, so no census, no cut and no cleanliness claim is defined on it |
| `20260226123353-auto_grown_20260226123353106` | PHerc0814 | 23x56 | 878 | 0 | the input carries no triangles, so no census, no cut and no cleanliness claim is defined on it |

## Fragmentation

174 of the 179 segments where the core gate applies pass it. The following show material fragmentation: at least one component inside the 99.9%-area core has R_main below 0.90, meaning its largest surviving descendant holds less than 90% of that component's input area. Area retention alone does not show this — each of these segments still retains most of its area — which is why R_main is reported separately.

| Segment | Scroll | min R_main (core) | min R_main (all) | area-weighted R_main | core components below gate |
|---|---|---:|---:|---:|---:|
| `20250925204843-auto_grown_20250925204843950` | PHerc0814 | 0.426323 | 0.426323 | 0.951641 | 1 |
| `20250925212701-auto_grown_20250925212701145` | PHerc0814 | 0.033906 | 0.033906 | 0.976250 | 1 |
| `20250925222237-auto_grown_20250925222237459` | PHerc0814 | 0.396828 | 0.396828 | 0.955513 | 1 |
| `20250926051122-auto_grown_20250926051122829` | PHerc0814 | 0.769789 | 0.769789 | 0.827777 | 5 |
| `20251001060526-auto_grown_20251001060526760` | PHerc0814 | 0.265911 | 0.265911 | 0.928370 | 2 |

## Archive size on disk

Sizes are of the emitted corpus artifacts only: the excised meshes, the per-segment certificates and the run logs. Driver scratch directories are not release artifacts and are excluded. 'on disk' is allocated blocks, the number a disk-usage tool reports; 'bytes' is the sum of file sizes.

**Total emitted corpus artifacts: 1,685,368,143 bytes = 1.570 GiB (1.685 GB); 1,688,145,920 bytes allocated on disk.**

| Scroll | Segments | Emitted meshes | Mesh bytes | Certificate bytes | Log bytes | Total |
|---|---:|---:|---:|---:|---:|---:|
| PHerc0139 | 38 | 26 | 47.4 MiB | 0.9 MiB | 38 KiB | 48.3 MiB |
| PHerc0814 | 19 | 11 | 44.5 MiB | 0.7 MiB | 19 KiB | 45.2 MiB |
| PHerc1667 | 20 | 20 | 156.5 MiB | 1.1 MiB | 25 KiB | 157.7 MiB |
| Scroll 1 | 55 | 46 | 935.3 MiB | 3.8 MiB | 60 KiB | 939.1 MiB |
| Scroll 5 | 53 | 51 | 413.2 MiB | 3.8 MiB | 67 KiB | 417.0 MiB |

### Proposed split

One archive per scroll. The five scrolls are independently useful, the split keeps every archive under 1 GiB, and each archive is self-describing because it carries a copy of the shared manifests. No archive is created by this script. The sizes below are the uncompressed payload each archive would contain.

| Proposed archive | Scroll | Contents | Payload |
|---|---|---|---:|
| `windcheck-corpus-pherc0139.tar` | PHerc0139 | 26 excised tifxyz meshes; 38 excision certificates; 38 run logs; shared manifests (base manifest, driver summary, independent recensus record, release index) | 50.0 MiB (0.049 GiB) |
| `windcheck-corpus-pherc0814.tar` | PHerc0814 | 11 excised tifxyz meshes; 19 excision certificates; 19 run logs; shared manifests (base manifest, driver summary, independent recensus record, release index) | 46.8 MiB (0.046 GiB) |
| `windcheck-corpus-pherc1667.tar` | PHerc1667 | 20 excised tifxyz meshes; 20 excision certificates; 20 run logs; shared manifests (base manifest, driver summary, independent recensus record, release index) | 159.3 MiB (0.156 GiB) |
| `windcheck-corpus-scroll1.tar` | Scroll 1 | 46 excised tifxyz meshes; 55 excision certificates; 55 run logs; shared manifests (base manifest, driver summary, independent recensus record, release index) | 940.8 MiB (0.919 GiB) |
| `windcheck-corpus-scroll5.tar` | Scroll 5 | 51 excised tifxyz meshes; 53 excision certificates; 53 run logs; shared manifests (base manifest, driver summary, independent recensus record, release index) | 418.7 MiB (0.409 GiB) |

Each archive also carries the shared manifests, 1.7 MiB in total: `out/corpus_bases.json`, `out/excised/corpus/corpus_summary.jsonl`, `out/excised/corpus/verification.json`, `out/release/index.json`.

## Where the certificates live

One excision certificate per segment, at `out/excised/corpus/<segment>_excision_certificate.json`. Each one records the pinned input mesh and its plane hashes, the emitted output mesh and its plane hashes, the census before and after under both canonical triangulations, both area denominators, the per-component recovery table, the frozen selection policy, and the code provenance of the run.

Alongside them:

- `out/corpus_bases.json` — the pinned base manifest: which mesh each segment was cut from, and the hash that pins it.
- `out/excised/corpus/corpus_summary.jsonl` — one driver record per segment.
- `out/excised/corpus/verification.json` — the independent recensus: every emitted artifact re-hashed and re-censused in a fresh workdir. Every census here ran in a fresh temporary workdir created for that segment and deleted afterwards; no driver workdir, CSV or atlas was reused.
- `out/excised/corpus/logs/<segment>.log` — the per-segment run log.
- `out/release/index.json` — this document, machine-readable.

## Per-segment results

`op.` is operational retention, `headline` is headline retention; see the two denominators above. Hashes are the first 12 hex digits of the sha256 of the `x` plane. A dash means the certificate records no value for that segment, which is stated rather than filled in.

### PHerc0139 (38 segments)

| Segment | Disposition | Base | Input hash | Output hash | op. retained | headline retained | min R_main | area-wt R_main | core gate | Runtime (s) |
|---|---|---|---|---|---:|---:|---:|---:|:---:|---:|
| `20250108000000-w025_2025010863` | already_clean | repaired | `94716d64373d` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.22 |
| `20250108000001-w026_2025010854` | transformed | repaired | `b8a4a7f6bcba` | `6b11d52729a5` | 99.8334% | 99.8334% | 0.9983 | 0.9983 | pass | 24.14 |
| `20250108000002-w027_2025010845` | transformed | repaired | `0aeb136dc878` | `7c079a695b6d` | 99.9304% | 99.9304% | 0.9993 | 0.9993 | pass | 3.61 |
| `20250108000003-w028_2025010836` | transformed | repaired | `3d814d1f2828` | `e4357e21ca53` | 99.9879% | 99.9879% | 0.9999 | 0.9999 | pass | 1.17 |
| `20250108000004-w029_2025010827` | transformed | repaired | `9f8eee169e06` | `9b0ab2523e69` | 99.9963% | 99.9963% | 1.0000 | 1.0000 | pass | 0.72 |
| `20250108000005-w030_2025010818` | transformed | repaired | `65680c76bdb3` | `506327e58e30` | 99.9847% | 99.9847% | 0.9998 | 0.9998 | pass | 0.82 |
| `20250223000000-w059_2025022312` | transformed | repaired | `4e074f19ab1c` | `7d879e468030` | 99.9980% | 99.9980% | 1.0000 | 1.0000 | pass | 0.90 |
| `20250831000000-w040_2025083102` | transformed | repaired | `3836d0f60afc` | `9cdc87338fd8` | 99.9912% | 99.9912% | 0.9999 | 0.9999 | pass | 0.74 |
| `20251226000000-w055_2025122611` | transformed | repaired | `b3346f0ea516` | `bdc56ecc998b` | 99.7160% | 99.7160% | 0.9971 | 0.9971 | pass | 143.75 |
| `20260108000000-w041_2026010816` | transformed | repaired | `f6d823c613cb` | `7a0c8e4eed07` | 99.9840% | 99.9840% | 0.9998 | 0.9998 | pass | 1.23 |
| `20260112000000-w043_2026011217` | transformed | repaired | `3cd9b13ba15a` | `1bca323305ff` | 99.9945% | 99.9945% | 0.9999 | 0.9999 | pass | 0.85 |
| `20260115000000-w044_2026011522` | transformed | repaired | `3be889015abc` | `24c8d2426b17` | 99.9989% | 99.9989% | 1.0000 | 1.0000 | pass | 0.81 |
| `20260115000001-w056_2026011514` | transformed | repaired | `f62790117bca` | `674821881971` | 99.9577% | 99.9577% | 0.9996 | 0.9996 | pass | 13.41 |
| `20260126000000-w045_2026012619` | transformed | original | `152250af01f4` | `2cf30bf95757` | 99.9969% | 99.9969% | 1.0000 | 1.0000 | pass | 0.88 |
| `20260127000000-w057_2026012713` | transformed | repaired | `0c2561586afc` | `f5fc9902e083` | 99.9411% | 99.9411% | 0.9994 | 0.9994 | pass | 33.97 |
| `20260130000000-w031_2026013019` | transformed | repaired | `6bc81d7078fb` | `934d0d389030` | 99.8262% | 99.8262% | 0.9982 | 0.9982 | pass | 31.71 |
| `20260203000000-w032_2026020303` | transformed | repaired | `40c033837b35` | `abd3e634da25` | 99.7398% | 99.7398% | 0.9972 | 0.9972 | pass | 83.73 |
| `20260206000000-w042_2026020613` | already_clean | original | `a8f89d90cc0c` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.44 |
| `20260206000001-w047_2026020613` | already_clean | repaired | `db3a8a55d43d` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.33 |
| `20260206000002-w054_2026020617` | transformed | repaired | `f5a6a1189855` | `34689c767f6d` | 99.4534% | 99.4534% | 0.9945 | 0.9945 | pass | 134.06 |
| `20260210000000-w058_2026021020` | transformed | original | `caa2552e84c0` | `aee4c198aee7` | 99.9988% | 99.9988% | 1.0000 | 1.0000 | pass | 0.78 |
| `20260215000000-w024_2026021572` | already_clean | repaired | `e202b7990c71` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.23 |
| `20260219000000-w033_2026021913` | transformed | repaired | `3488357aad35` | `8608703f0803` | 99.8076% | 99.8076% | 0.9980 | 0.9980 | pass | 67.44 |
| `20260219000001-w048_2026021912` | already_clean | original | `095a1a7a38ac` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.19 |
| `20260220000000-w053_2026022021` | transformed | repaired | `17d29419ae84` | `91f35980bfbe` | 99.9527% | 99.9527% | 0.9995 | 0.9995 | pass | 7.61 |
| `20260227000000-w052_2026022705` | transformed | repaired | `5adaa195db5f` | `bb1814025420` | 99.7814% | 99.7814% | 0.9969 | 0.9969 | pass | 90.23 |
| `20260302000000-w039_2026030210` | already_clean | original | `b1915cd9722f` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.23 |
| `20260302000001-w050_2026030220` | already_clean | original | `bb3f5df9fe8c` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.35 |
| `20260303000000-w034_2026030317` | transformed | repaired | `27e21dbf0559` | `0f569faa024d` | 99.8462% | 99.8462% | 0.9981 | 0.9981 | pass | 54.35 |
| `20260306000000-w038_2026030608` | already_clean | original | `f5ff688f66a2` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.22 |
| `20260306000001-w051_2026030600` | transformed | original | `dff240d4c984` | `12172698d009` | 99.7527% | 99.7527% | 0.9975 | 0.9975 | pass | 90.95 |
| `20260310000000-w037_2026031015` | already_clean | original | `7a8817ced2b4` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.22 |
| `20260311000000-w023_202603112281` | already_clean | original | `ac664704d9b2` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.21 |
| `20260316000000-w036_2026031607` | already_clean | original | `733664739b47` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.20 |
| `20260317000000-w035_2026031718` | already_clean | original | `d827a5d2e510` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.19 |
| `20260318000000-w049_20260318` | transformed | repaired | `868d0e21c75e` | `797339b48f78` | 99.9162% | 99.9162% | 0.9990 | 0.9990 | pass | 15.50 |
| `20260325000000-w046_20260325` | transformed | original | `3ab134524dc4` | `b5cb7494f9a1` | 99.9922% | 99.9922% | 0.9999 | 0.9999 | pass | 1.00 |
| `20260422000000-title_2026042222_zmid_flatboi` | transformed | repaired | `cfc964f5f10c` | `0e1e5b7cb296` | 99.9103% | 99.9103% | 0.9990 | 0.9990 | pass | 2.17 |

### PHerc0814 (19 segments)

| Segment | Disposition | Base | Input hash | Output hash | op. retained | headline retained | min R_main | area-wt R_main | core gate | Runtime (s) |
|---|---|---|---|---|---:|---:|---:|---:|:---:|---:|
| `20250925161630-auto_grown_20250925161630635` | transformed | repaired | `54fd7fd3fc64` | `4d186e7d43c6` | 99.9848% | 99.9848% | 0.9985 | 0.9998 | pass | 1.07 |
| `20250925182632-auto_grown_20250925182632715` | transformed | repaired | `3b9805249577` | `2e8f9a32aa5b` | 98.9933% | 98.9933% | 0.9890 | 0.9893 | pass | 135.41 |
| `20250925204843-auto_grown_20250925204843950` | transformed | repaired | `ff2c07c3607d` | `ff6b12a22736` | 98.9331% | 98.9331% | 0.4263 | 0.9516 | FAIL | 66.20 |
| `20250925212701-auto_grown_20250925212701145` | transformed | original | `d80ab2ed2d6d` | `513af8e0278e` | 97.8985% | 97.8985% | 0.0339 | 0.9763 | FAIL | 40.80 |
| `20250925222237-auto_grown_20250925222237459` | transformed | repaired | `55fbdc6e89d7` | `c3a93d2be773` | 96.7141% | 96.7141% | 0.3968 | 0.9555 | FAIL | 38.02 |
| `20250926051122-auto_grown_20250926051122829` | transformed | original | `344005063f1f` | `308140334088` | 96.8136% | 96.8136% | 0.7698 | 0.8278 | FAIL | 134.47 |
| `20250926165636-auto_trace_20250926165636789` | transformed | repaired | `af112ca5b651` | `53f6b5348a60` | 99.8164% | 99.8164% | 0.9966 | 0.9981 | pass | 64.71 |
| `20250926203608-auto_grown_20250926203608932` | transformed | repaired | `acc16975cbc9` | `e60f7c8307d0` | 99.3221% | 99.3221% | 0.9081 | 0.9665 | pass | 54.16 |
| `20250928235953-auto_grown_20250928235953722` | transformed | repaired | `4a752c7aeb04` | `ebcda3c78894` | 99.5650% | 99.5650% | 0.9388 | 0.9926 | pass | 122.68 |
| `20250928235954-auto_grown_20250928235953722_copy` | not_censusable | original | `e9c6fe49fd53` | `—` | not recorded | not recorded | — | — | — | 0.07 |
| `20250929040153-auto_grown_20250929040153260` | already_clean | original | `bf640eb92993` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.15 |
| `20250930083742-auto_grown_20250930083742882` | transformed | repaired | `004366eed5a2` | `33c733256602` | 99.7325% | 99.7326% | 0.9427 | 0.9972 | pass | 37.43 |
| `20251001060526-auto_grown_20251001060526760` | transformed | original | `282a016ae487` | `f96bc785ad8b` | 94.6785% | 94.6785% | 0.2659 | 0.9284 | FAIL | 124.28 |
| `20260225160055-auto_grown_20260225160055913` | not_censusable | original | `f1466bdafaa7` | `—` | not recorded | not recorded | — | — | — | 0.10 |
| `20260226000000-46527_2um_try2` | already_clean | original | `1f0748d26fce` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.12 |
| `20260226110425-auto_grown_20260226110425097_abf` | not_censusable | original | `b6552bbe0bbb` | `—` | not recorded | not recorded | — | — | — | 0.09 |
| `20260226113324-auto_grown_20260226113324624_abf` | not_censusable | original | `3edefa20c5a8` | `—` | not recorded | not recorded | — | — | — | 0.08 |
| `20260226121302-auto_grown_20260226121302716_abf` | not_censusable | original | `b4923302cdc6` | `—` | not recorded | not recorded | — | — | — | 0.08 |
| `20260226123353-auto_grown_20260226123353106` | not_censusable | original | `e54b654a86ae` | `—` | not recorded | not recorded | — | — | — | 0.09 |

### PHerc1667 (20 segments)

| Segment | Disposition | Base | Input hash | Output hash | op. retained | headline retained | min R_main | area-wt R_main | core gate | Runtime (s) |
|---|---|---|---|---|---:|---:|---:|---:|:---:|---:|
| `20240304141531-w013_20240304141531_flatboi` | transformed | repaired | `16f75638245b` | `94fb78403055` | 99.9986% | 99.9986% | 1.0000 | 1.0000 | pass | 2.79 |
| `20240304144031-w018_20240304144031_flatboi` | transformed | repaired | `e86f86fd79d9` | `236f874cea09` | 99.9997% | 99.9997% | 1.0000 | 1.0000 | pass | 3.09 |
| `20240304161941-w023_20240304161941_flatboi` | transformed | repaired | `b60c17fd15ee` | `0e26e503d4d0` | 99.9987% | 99.9987% | 1.0000 | 1.0000 | pass | 3.09 |
| `20251206103305-w012_20251206103305555_flatboi` | transformed | repaired | `486f4111aaba` | `36a8dd2216ac` | 99.9589% | 99.9589% | 0.9995 | 0.9995 | pass | 3.35 |
| `20251208130119-w028_20251208130119156_flatboi` | transformed | repaired | `aa90cc226789` | `c6be5b6e7b18` | 99.9649% | 99.9649% | 0.9996 | 0.9996 | pass | 32.39 |
| `20251212185248-w029_20251212185248662_flatboi` | transformed | repaired | `281e870e2486` | `c7e9235dc85f` | 99.9876% | 99.9876% | 0.9999 | 0.9999 | pass | 2.08 |
| `20251220020000-w030_2025122002_flatboi` | transformed | repaired | `6a7dbc67bfc9` | `85c3c2d5644b` | 99.9478% | 99.9478% | 0.9994 | 0.9994 | pass | 21.61 |
| `20251223230000-w031_2025122323_flatboi` | transformed | repaired | `cb93c0cb4deb` | `ffa654878400` | 99.9407% | 99.9407% | 0.9994 | 0.9994 | pass | 33.40 |
| `20260105050000-w032_2026010505_flatboi` | transformed | repaired | `0726b2e748b7` | `9c9a7bd18719` | 99.8279% | 99.8279% | 0.9982 | 0.9982 | pass | 46.84 |
| `20260108140509-w011_20260108140509268_flatboi` | transformed | repaired | `3c91fba2be5f` | `2188a88d3024` | 99.7684% | 99.7684% | 0.9977 | 0.9977 | pass | 16.63 |
| `20260110010000-w033_2026011001_flatboi` | transformed | repaired | `5b076ddcb957` | `70f07f88c021` | 99.9815% | 99.9815% | 0.9998 | 0.9998 | pass | 1.72 |
| `20260115160000-w034_2026011516_flatboi` | transformed | repaired | `d531fc39824f` | `e6f18b69e9a6` | 99.9300% | 99.9300% | 0.9992 | 0.9992 | pass | 11.58 |
| `20260116230000-w035_2026011623_flatboi` | transformed | repaired | `3d316dc40f0c` | `155c7ef63dc5` | 99.6914% | 99.6914% | 0.9969 | 0.9969 | pass | 126.35 |
| `20260119120000-w036_2026011912_flatboi` | transformed | repaired | `58410f093dcc` | `1b2b39f5cb67` | 99.7202% | 99.7202% | 0.9972 | 0.9972 | pass | 105.09 |
| `20260123230000-w037_2026012323_flatboi` | transformed | repaired | `df2326bd9305` | `152cc6181959` | 99.9171% | 99.9171% | 0.9992 | 0.9992 | pass | 35.31 |
| `20260128140000-w038_2026012814_flatboi` | transformed | repaired | `c9d19d50f252` | `aa14bfe85aaa` | 99.9108% | 99.9108% | 0.9991 | 0.9991 | pass | 6.78 |
| `20260130150000-w039_2026013015_flatboi` | transformed | repaired | `4b4d2f4109aa` | `209c8be6c644` | 99.8610% | 99.8610% | 0.9986 | 0.9986 | pass | 25.66 |
| `20260203210000-w040_2026020321_flatboi` | transformed | repaired | `a52f69756bb4` | `26b0e5387bd1` | 99.9564% | 99.9564% | 0.9996 | 0.9996 | pass | 3.54 |
| `20260205070000-w041_2026020507_flatboi` | transformed | repaired | `294b307eae52` | `85e5ad46953c` | 99.7599% | 99.7599% | 0.9976 | 0.9976 | pass | 42.78 |
| `20260612121456-w011_20260108140509268_merged_v4_flatboi_straightened_v4` | transformed | original | `3b880bd1d1de` | `29e4488971f7` | 96.1655% | 96.1655% | 0.0000 | 0.9576 | pass | 197.33 |

### Scroll 1 (55 segments)

| Segment | Disposition | Base | Input hash | Output hash | op. retained | headline retained | min R_main | area-wt R_main | core gate | Runtime (s) |
|---|---|---|---|---|---:|---:|---:|---:|:---:|---:|
| `20230702185753` | already_clean | original | `4a4ce60cd6a4` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.78 |
| `20230929220926` | already_clean | original | `1f4ed8800744` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 1.04 |
| `20231005123336` | transformed | original | `d3e8b99d335f` | `18364329508f` | 99.9992% | 99.9992% | 1.0000 | 1.0000 | pass | 3.12 |
| `20231007101619` | already_clean | original | `9696d0abb7b8` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 2.09 |
| `20231012184424` | already_clean | repaired | `fe1f2f5ef802` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 1.22 |
| `20231016151002` | transformed | original | `ce832c2e630d` | `4e282d2a9a49` | 99.9886% | 99.9886% | 0.9999 | 0.9999 | pass | 2.87 |
| `20231022170901` | already_clean | original | `72e0594a586a` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 1.15 |
| `20231031143852` | already_clean | original | `f537a9030c48` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.64 |
| `20231106155351` | already_clean | repaired | `728197d292cc` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.61 |
| `20231210121321` | already_clean | original | `e7e250ca7a21` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.70 |
| `20231221180251` | already_clean | original | `869823176112` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.60 |
| `20260602204401-5753_-7` | transformed | repaired | `5df8aead3af7` | `900b165961ec` | 99.9547% | 99.9547% | 0.9995 | 0.9995 | pass | 2.25 |
| `20260602225659-5753_0` | transformed | repaired | `3029dc519b72` | `48742963c2f3` | 99.7785% | 99.7785% | 0.9977 | 0.9977 | pass | 108.31 |
| `20260603005223-5753_-1` | transformed | original | `f7247bbd3ef5` | `3d3e8702de65` | 99.6893% | 99.6893% | 0.9965 | 0.9965 | pass | 122.30 |
| `20260603024952-5753_-2` | transformed | repaired | `6ac955ebf6d1` | `ff9d049cd232` | 99.8651% | 99.8651% | 0.9985 | 0.9985 | pass | 77.08 |
| `20260603042357-5753_-3` | transformed | repaired | `8088ad6ea36c` | `4a13927eb3a4` | 99.8823% | 99.8823% | 0.9988 | 0.9988 | pass | 48.96 |
| `20260603145540-5753_-4` | transformed | original | `47d54bbf39e4` | `602cf9d8d898` | 99.5819% | 99.5819% | 0.9958 | 0.9958 | pass | 125.78 |
| `20260603185441-5753_-6` | transformed | original | `3fe958486d37` | `0a9123fe68b3` | 99.6323% | 99.6323% | 0.9962 | 0.9962 | pass | 108.04 |
| `20260603190005-5753_-5` | transformed | repaired | `4bdb7410c576` | `082bdbbe127f` | 99.7975% | 99.7975% | 0.9979 | 0.9979 | pass | 44.28 |
| `20260603193311-5753_-1_copy` | transformed | repaired | `5f24c7765569` | `022629c1fb90` | 99.9125% | 99.9125% | 0.9989 | 0.9989 | pass | 37.50 |
| `20260603222816-20231005123336_v2_flatboi` | transformed | repaired | `86adab33a90c` | `35bb4c4d5e5c` | 99.5670% | 99.5670% | 0.9948 | 0.9948 | pass | 144.07 |
| `20260604223808-20231210121321_v8` | transformed | repaired | `6d484e0da1a9` | `a72e994e78f6` | 99.7318% | 99.7318% | 0.0000 | 0.9969 | pass | 124.05 |
| `20260623141135-w046-052_jordi` | transformed | original | `f3f9f010bc11` | `66c5dfcf3bfb` | 99.9709% | 99.9709% | 0.9997 | 0.9997 | pass | 75.07 |
| `20260623141649-w053-058_jordi` | transformed | original | `26fdbeb9f03e` | `eaa39e2e72db` | 99.9411% | 99.9411% | 0.0000 | 0.9994 | pass | 137.35 |
| `20260623141924-w010-027` | transformed | original | `fbdbd52efc25` | `5cebf8dfd454` | 99.9714% | 99.9714% | 0.9997 | 0.9997 | pass | 110.12 |
| `20260623142658-w028-037` | transformed | repaired | `2b2168ce0fdb` | `49b3323fde42` | 99.9897% | 99.9897% | 0.9999 | 0.9999 | pass | 33.28 |
| `20260623143441-w038-045` | transformed | original | `f192e7d4db99` | `4cf93a06b7f9` | 99.9710% | 99.9710% | 0.9997 | 0.9997 | pass | 66.00 |
| `20260623144224-w046-052` | transformed | original | `7a60edfe74c4` | `aeb3a167c329` | 99.9585% | 99.9585% | 0.9996 | 0.9996 | pass | 120.27 |
| `20260623144957-w053-058` | transformed | repaired | `1183316c1e3a` | `89dc82eacd6d` | 99.9668% | 99.9668% | 0.4237 | 0.9996 | pass | 33.12 |
| `20260623145652-w059-063` | transformed | original | `73137efb6f50` | `1f2154869a59` | 99.9136% | 99.9136% | 0.0000 | 0.9990 | pass | 35.89 |
| `20260623150417-w064-068` | transformed | original | `5a67a56c987f` | `88ab6d2a0291` | 99.9790% | 99.9790% | 0.9998 | 0.9998 | pass | 7.26 |
| `20260623151041-w069-072` | transformed | repaired | `3b101312a575` | `c1fdd4bdb282` | 99.9973% | 99.9973% | 1.0000 | 1.0000 | pass | 6.69 |
| `20260623151729-w073-076` | transformed | repaired | `9b3f2ff759bf` | `71302a5aa91e` | 99.9935% | 99.9935% | 0.9999 | 0.9999 | pass | 13.12 |
| `20260623152443-w077-080` | transformed | original | `7f267e046fdf` | `7c22c260df4d` | 99.9911% | 99.9911% | 0.9999 | 0.9999 | pass | 36.77 |
| `20260623153216-w081-084` | transformed | repaired | `397268a7e43e` | `1b9bac84c985` | 99.9815% | 99.9815% | 0.9998 | 0.9998 | pass | 47.68 |
| `20260623154006-w085-088` | transformed | repaired | `8957cfd7729a` | `f0d8e7c9887d` | 99.9886% | 99.9886% | 0.9999 | 0.9999 | pass | 31.35 |
| `20260623154617-w089-091` | transformed | original | `408e81e00429` | `c92dbbfa004d` | 99.9876% | 99.9876% | 0.9999 | 0.9999 | pass | 18.66 |
| `20260623155240-w092-094` | transformed | repaired | `5b866d55cd15` | `1e4dea6a4cbf` | 99.9858% | 99.9858% | 0.9999 | 0.9999 | pass | 30.02 |
| `20260623155914-w095-097` | transformed | original | `ddb42664241d` | `1d11197f74e7` | 99.9816% | 99.9816% | 0.9998 | 0.9998 | pass | 48.83 |
| `20260623160554-w098-100` | transformed | repaired | `b7569cb69b45` | `cf0c35813582` | 99.9832% | 99.9832% | 0.9998 | 0.9998 | pass | 27.98 |
| `20260623161233-w101-103` | transformed | original | `b0a22412f8c7` | `abdc287f67d7` | 99.9645% | 99.9645% | 0.9996 | 0.9996 | pass | 60.50 |
| `20260623161921-w104-106` | transformed | original | `834d93fc6801` | `29d8814dfea7` | 99.9585% | 99.9585% | 0.9996 | 0.9996 | pass | 136.56 |
| `20260623162614-w107-109` | transformed | original | `610a59fed29a` | `989f1d5c3337` | 99.9300% | 99.9300% | 0.9993 | 0.9993 | pass | 122.53 |
| `20260623163339-w110-112` | transformed | original | `822574be19b5` | `e89911a2d09a` | 99.8794% | 99.8794% | 0.9988 | 0.9988 | pass | 100.32 |
| `20260623164141-w113-115` | transformed | original | `72e30042919f` | `556ea3efe561` | 99.8035% | 99.8035% | 0.9980 | 0.9980 | pass | 101.74 |
| `20260623164704-w116-117` | transformed | original | `9fc0456d629f` | `ef8ed96cb047` | 99.7929% | 99.7929% | 0.9979 | 0.9979 | pass | 85.16 |
| `20260623165230-w118-119` | transformed | original | `354306d22f94` | `5d6269171512` | 99.7177% | 99.7177% | 0.9971 | 0.9971 | pass | 130.47 |
| `20260623165742-w120-121` | transformed | original | `bb03be5dc8f9` | `4cebaa37a48e` | 99.6368% | 99.6368% | 0.9963 | 0.9963 | pass | 129.46 |
| `20260623170305-w122-123` | transformed | original | `7a9eae441164` | `e7cb1efc2222` | 99.5602% | 99.5602% | 0.9955 | 0.9955 | pass | 137.42 |
| `20260623170833-w124-125` | transformed | original | `06a59e1054fe` | `a5484fc0a160` | 99.4768% | 99.4768% | 0.0000 | 0.9945 | pass | 137.05 |
| `20260623171400-w126-127` | transformed | original | `060ca39b22c2` | `3a2b2780f9ba` | 98.6858% | 98.6858% | 0.0000 | 0.9853 | pass | 191.46 |
| `20260623171929-w128-129` | transformed | original | `fab16b1f282c` | `25c30de59703` | 96.1026% | 96.1026% | 0.0000 | 0.9497 | pass | 164.93 |
| `20260701183124-w010-027` | transformed | original | `47a5dbfd5af9` | `4edc67b9ed28` | 99.9448% | 99.9448% | 0.9994 | 0.9994 | pass | 130.63 |
| `20260701183125-w028-037` | transformed | repaired | `f19aea6fe109` | `7cac4c2a331f` | 99.9262% | 99.9262% | 0.9993 | 0.9993 | pass | 122.90 |
| `20260701183126-w038-045` | transformed | original | `68609f1e9d98` | `c0a1bab0cf99` | 99.8799% | 99.8799% | 0.9988 | 0.9988 | pass | 138.15 |

### Scroll 5 (53 segments)

| Segment | Disposition | Base | Input hash | Output hash | op. retained | headline retained | min R_main | area-wt R_main | core gate | Runtime (s) |
|---|---|---|---|---|---:|---:|---:|---:|:---:|---:|
| `20250917143559-w062_20250917143559205_flatboi` | already_clean | original | `168f56c764ab` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.39 |
| `20250926112011-w078_20250926112011918_flatboi` | transformed | original | `dd0ddf0dd8e2` | `92adad3d0c48` | 99.0070% | 99.0070% | 0.9900 | 0.9900 | pass | 131.04 |
| `20250926113336-w079_20250926113336891_flatboi` | transformed | original | `2a163dd87cd3` | `17c26c20de18` | 98.7633% | 98.7633% | 0.9875 | 0.9875 | pass | 151.22 |
| `20250926114310-w080_20250926114310705_flatboi` | transformed | repaired | `c3d5c0d25201` | `cca8c237402b` | 98.4900% | 98.4900% | 0.9847 | 0.9847 | pass | 130.05 |
| `20250926122938-w082_20250926122938619_flatboi` | transformed | original | `fff82ff301d9` | `f4fe53395078` | 97.2069% | 97.2069% | 0.9694 | 0.9694 | pass | 159.24 |
| `20250926132117-w083_20250926132117595_flatboi` | transformed | original | `fe21c00c5316` | `e84585287a9c` | 98.2309% | 98.2309% | 0.9822 | 0.9822 | pass | 136.24 |
| `20250926145427-w084_20250926145427439_flatboi` | transformed | repaired | `8aff22ef6ec9` | `11bba9d644e3` | 99.7029% | 99.7029% | 0.9970 | 0.9970 | pass | 142.70 |
| `20250926150345-w086_20250926150345176_flatboi` | transformed | repaired | `61c065b1582f` | `db3068466824` | 99.9004% | 99.9004% | 0.9990 | 0.9990 | pass | 61.66 |
| `20250926153406-w088_20250926153406704_flatboi` | transformed | repaired | `63ba269815c6` | `fc7f60bf8bb6` | 99.7712% | 99.7712% | 0.0000 | 0.9976 | pass | 122.38 |
| `20251106133152-w077_20251106133152767_flatboi` | transformed | repaired | `c822c673f629` | `f7691c474044` | 99.2531% | 99.2530% | 0.9924 | 0.9924 | pass | 125.28 |
| `20251106170358-w076_20251106170358232_flatboi` | transformed | repaired | `7a039ebd08d1` | `1ef6ad5413c3` | 98.7539% | 98.7539% | 0.9873 | 0.9873 | pass | 150.30 |
| `20251107063327-w063_20251107063327215_flatboi` | transformed | repaired | `c30f246f1e38` | `791108ab0242` | 99.9870% | 99.9870% | 0.9999 | 0.9999 | pass | 12.67 |
| `20251107110950-w064_20251107110950052_flatboi` | transformed | original | `39e8c4dd8e03` | `3fbc7cb4af77` | 99.9273% | 99.9273% | 0.9993 | 0.9993 | pass | 66.93 |
| `20251107175958-w061_20251107175958270_flatboi` | transformed | original | `6d7d4471ace0` | `a2c9a57d6ec9` | 99.9920% | 99.9920% | 0.9999 | 0.9999 | pass | 7.31 |
| `20251109232817-w065_20251109232817724_flatboi` | transformed | original | `6c424c838f2b` | `117fc0d4edff` | 99.9063% | 99.9063% | 0.9990 | 0.9990 | pass | 22.93 |
| `20251110110112-w066_20251110110112608_flatboi` | transformed | original | `54a78c07d17e` | `2a326a2c2c50` | 99.8972% | 99.8972% | 0.9989 | 0.9989 | pass | 41.04 |
| `20251110135146-w085_20251110135146457_flatboi` | transformed | original | `dc6c85ff9d78` | `544d396bf5b9` | 99.6158% | 99.6158% | 0.9959 | 0.9959 | pass | 142.16 |
| `20251110135803-w067_20251110135803677_flatboi` | transformed | repaired | `4ddcf6d1cf7d` | `019381a89e8d` | 99.9484% | 99.9484% | 0.9995 | 0.9995 | pass | 61.96 |
| `20251111010954-w068_20251111010954408_flatboi` | transformed | repaired | `a15272768639` | `ad74a870838e` | 99.9441% | 99.9441% | 0.9994 | 0.9994 | pass | 87.75 |
| `20251111151836-w069_20251111151836406_flatboi` | transformed | original | `5389f3983182` | `c045fb715b71` | 98.3263% | 98.3263% | 0.9829 | 0.9829 | pass | 152.99 |
| `20251111210734-w070_20251111210734722_flatboi` | transformed | repaired | `93d1ebbf2506` | `33557bd664c8` | 99.9363% | 99.9363% | 0.9994 | 0.9994 | pass | 9.33 |
| `20251112000002-w087_20251112000002214_flatboi` | transformed | repaired | `2fe747e90785` | `7b125d1ffbee` | 99.9557% | 99.9557% | 0.9995 | 0.9995 | pass | 41.19 |
| `20251112163527-w071_20251112163527854_flatboi` | transformed | repaired | `953790527a05` | `e4542117bd5b` | 99.9691% | 99.9691% | 0.9997 | 0.9997 | pass | 20.74 |
| `20251113021010-w072_20251113021010674_flatboi` | transformed | repaired | `c083bda8a6b8` | `69107b467017` | 99.9788% | 99.9788% | 0.9998 | 0.9998 | pass | 18.55 |
| `20251113130119-w073_20251113130119516_flatboi` | transformed | repaired | `ed513c5beb1b` | `15bfd9f97f4d` | 99.9798% | 99.9798% | 0.9998 | 0.9998 | pass | 32.80 |
| `20251113152136-w074_20251113152136889_flatboi` | transformed | repaired | `4c2e7bf2a3b5` | `236cc978d421` | 99.9953% | 99.9953% | 1.0000 | 1.0000 | pass | 1.56 |
| `20251113204233-w075_20251113204233755_flatboi` | transformed | repaired | `04fe0d840ec7` | `a8da14a7240b` | 99.9794% | 99.9794% | 0.9998 | 0.9998 | pass | 33.99 |
| `20251114202035-w081_20251114202035261_flatboi` | transformed | original | `b242a927a819` | `6b27db2b162d` | 99.4811% | 99.4811% | 0.9947 | 0.9947 | pass | 115.04 |
| `20251115002740-auto_grown_20251115002740308_0_flatboi` | transformed | original | `d211ddd26330` | `929924d2c06e` | 99.7350% | 99.7350% | 0.9972 | 0.9972 | pass | 136.82 |
| `20251115002741-auto_grown_20251115002740308_1_flatboi` | transformed | original | `b8b3f94fc8c6` | `f580d9273006` | 99.3888% | 99.3888% | 0.0000 | 0.9897 | pass | 140.65 |
| `20251115002742-auto_grown_20251115002740308_2_flatboi` | transformed | original | `bed937a16f0e` | `508d49cbfd45` | 99.2158% | 99.2158% | 0.0000 | 0.9917 | pass | 211.26 |
| `20251115002743-auto_grown_20251115002740308_3_flatboi` | transformed | original | `a1497be20629` | `af691f0c545e` | 99.4617% | 99.4617% | 0.9945 | 0.9945 | pass | 150.08 |
| `20251115002744-auto_grown_20251115002740308_4_flatboi` | transformed | original | `dfb697a92116` | `f8c4644807f4` | 97.6050% | 97.6050% | 0.9754 | 0.9754 | pass | 203.88 |
| `20251115002745-auto_grown_20251115002740308_5_flatboi` | transformed | original | `e31d05db540e` | `e4053897a6fc` | 97.9247% | 97.9247% | 0.9788 | 0.9788 | pass | 205.29 |
| `20251115002746-auto_grown_20251115002740308_6_flatboi` | transformed | original | `f3ecfe26f9d2` | `b58b9a587660` | 99.1892% | 99.1892% | 0.9916 | 0.9916 | pass | 148.94 |
| `20251115002747-auto_grown_20251115002740308_7_flatboi` | transformed | original | `c429597f6a24` | `bf241f411986` | 97.9415% | 97.9415% | 0.9785 | 0.9785 | pass | 195.94 |
| `20251115002748-auto_grown_20251115002740308_8_flatboi` | transformed | original | `634d5862857a` | `1bd5907720e6` | 99.0508% | 99.0508% | 0.9904 | 0.9904 | pass | 137.58 |
| `20251203230500-w089_20251203230500465_flatboi` | transformed | repaired | `106c1e8518be` | `6c211d0e8f00` | 99.7993% | 99.7993% | 0.9980 | 0.9980 | pass | 130.30 |
| `20251204060522-w090_20251204060522416_flatboi` | transformed | repaired | `8cb73989bfb2` | `c4fca9c2fa72` | 99.8397% | 99.8397% | 0.9984 | 0.9984 | pass | 141.41 |
| `20251204071911-w091_20251204071911070_flatboi` | transformed | repaired | `3f791c9d8adc` | `bde6c7ba395c` | 99.8531% | 99.8531% | 0.9985 | 0.9985 | pass | 138.96 |
| `20251204101314-w092_20251204101314959_flatboi` | transformed | repaired | `c1d88224e50a` | `5545a1585ef2` | 99.7983% | 99.7983% | 0.9979 | 0.9979 | pass | 133.42 |
| `20251204192708-w093_20251204192708445_flatboi` | transformed | repaired | `602c4277e77d` | `d50901233ff2` | 99.7562% | 99.7562% | 0.9975 | 0.9975 | pass | 116.94 |
| `20251205115859-w094_20251205115859448_flatboi` | already_clean | repaired | `84ed604a14ae` | `—` | 100.0000% | 100.0000% | 1 by constr. | 1 by constr. | pass | 0.36 |
| `20251206184729-w095_20251206184729906_flatboi` | transformed | repaired | `3638f21ea3fb` | `986a991d996a` | 99.5290% | 99.5290% | 0.9951 | 0.9951 | pass | 130.82 |
| `20251207123912-w060_20251207123912611_flatboi` | transformed | repaired | `62a675dc3999` | `e32bedf49e9f` | 99.9950% | 99.9950% | 0.9999 | 0.9999 | pass | 1.34 |
| `20251210142414-w059_20251210142414060_flatboi` | transformed | repaired | `328cd7296c05` | `a281de503e3d` | 99.9601% | 99.9601% | 0.9996 | 0.9996 | pass | 18.20 |
| `20251210153304-w058_20251210153304503_flatboi` | transformed | repaired | `82fe742a2c64` | `1abfe96394d5` | 99.9343% | 99.9343% | 0.9993 | 0.9993 | pass | 63.22 |
| `20251210160810-w057_20251210160810301_flatboi` | transformed | repaired | `e8f4b8e16e57` | `43873012c5fa` | 99.9373% | 99.9373% | 0.9994 | 0.9994 | pass | 61.59 |
| `20251211095724-w056_20251211095724922_flatboi` | transformed | repaired | `b9695d444d9c` | `3ee0bfa9345f` | 99.9332% | 99.9332% | 0.9993 | 0.9993 | pass | 61.37 |
| `20251211111451-w055_20251211111451620_flatboi` | transformed | repaired | `b381d983163a` | `4222f011a629` | 99.9057% | 99.9057% | 0.9990 | 0.9990 | pass | 61.63 |
| `20251211135303-w054_20251211135303394_flatboi` | transformed | repaired | `436b8e601ee5` | `f7097122c5e7` | 99.9532% | 99.9532% | 0.9995 | 0.9995 | pass | 52.42 |
| `20251211152543-w053_20251211152543204_flatboi` | transformed | repaired | `29cd56e215a3` | `be76ae607d14` | 99.9871% | 99.9871% | 0.9999 | 0.9999 | pass | 2.38 |
| `20251212105627-w052_20251212105627148_flatboi` | transformed | repaired | `de0c93086ba5` | `43238e3d5931` | 99.9817% | 99.9817% | 0.9998 | 0.9998 | pass | 4.95 |

### Reading the retention columns

A segment with disposition `already_clean` shows 100.0000% under both denominators because no cut was made. Its R_main is 1 by construction rather than by measurement, and the certificate records no per-component table for it; that is marked `1 by constr.` rather than reported as a measured value. A segment with disposition `not_censusable` shows a dash everywhere: no area and no component structure is defined on an input with no triangles.

6 already-clean segments sit on a displacement-repaired base, so their certificate records no original-coordinate area: nothing was cut, so nothing was priced. They are **not** dropped from the area-weighted headline denominator. The canonical area of each one's original published mesh was recomputed from that mesh, by the same route every cut segment's headline denominator is computed, and each segment enters at a retained fraction of exactly 1.0. The recomputed areas are memoised in `out/headline_original_areas.json`. They are: `20250108000000-w025_2025010863`, `20260206000001-w047_2026020613`, `20260215000000-w024_2026021572`, `20231012184424`, `20231106155351`, `20251205115859-w094_20251205115859448_flatboi`.

6 not-censusable inputs are inside the count the headline figure is reported over, carrying zero area on both sides of the ratio: they have no triangles, so there is no surface to retain and none to remove, and their presence cannot move the figure either way. They are: `20250928235954-auto_grown_20250928235953722_copy`, `20260225160055-auto_grown_20260225160055913`, `20260226110425-auto_grown_20260226110425097_abf`, `20260226113324-auto_grown_20260226113324624_abf`, `20260226121302-auto_grown_20260226121302716_abf`, `20260226123353-auto_grown_20260226123353106`.

Restricted to the 154 transformed segments alone — that is, excluding every already-clean segment from both numerator and denominator — the area-weighted headline retention is 99.4756%. That is the stricter reading: it prices only the segments that were actually cut.

