# Reproducing the ScrollFiesta case study

Everything below runs in ScrollFiesta's own build environment (a Linux
container with their `build-deps.sh` toolchain) plus this repository's
Python environment. Expected wall clock per fixture: ~2 min carve,
~21 min `grid_pipeline`, ~2 min `scroll_whole`, ~10 min for the five
staged exports, ~30 min for the censuses on a laptop.

## 0. One-time setup

```sh
git clone https://github.com/Hob3rMallow/scrollfiesta_public scrollfiesta
cd scrollfiesta
git apply /path/to/windcheck/case-studies/scrollfiesta-pherc0139-4x5x5/patches/scrollfiesta_instrumentation_and_portability.diff
./build-deps.sh
cmake -S . -B build -DSCROLLFIESTA_BUILD_TOOLS=ON -DCMAKE_BUILD_TYPE=Release -GNinja
cmake --build build --parallel
```

The patch (against their MIT-licensed source, attribution preserved)
adds two `#include <stdlib.h>` portability fixes and two env-gated,
read-only diagnostics that are byte-neutral when the variables are
unset: `SF_GG_DUMP=<path>` (registration-graph dump) and
`SF_TXZ_FACEMAP=<prefix>` (per-pixel first-cover face provenance in the
tifxyz export).

## 1. Carve the fixture (their script, anonymous S3)

Canonical window:

```sh
python python/scripts/carve_grid_tifs.py \
  --pred-zarr "s3://vesuvius-challenge-open-data/PHerc0139/representations/predictions/surfaces/20250728140407-surface-20260413222639-surface-m7-L0-th0.2.zarr" \
  --raw-zarr  "s3://vesuvius-challenge-open-data/PHerc0139/volumes/20250728140407-9.362um-1.2m-113keV-masked.zarr" \
  --bbox 4352 4864 3072 3712 2560 3200 --umbilicus 3405 2878 \
  --s3-anon yes --out PHerc0139-4x5x5
```

Held-out window: identical except `--bbox 4864 5376 3072 3712 2560 3200`
and `--out PHerc0139-heldout`.

## 2. Their pipeline, their canonical flags

```sh
./build/grid_pipeline PHerc0139-4x5x5 output/run --halo 13 --exe build/cube_mesh
./build/scroll_whole output/run/dump output/run_placed
./build/scroll_whole output/run_placed --reregister --audit
```

Expected: 99/100 cubes meshed (one cube crashes deterministically under
Linux/GCC; their Windows reference run meshes it), audit gate PASS at
turn-off 2.92% (canonical) / 3.48% (held-out).

## 3. Staged exports with the face-provenance tap

```sh
for S in 0 12 123 1234 12345; do
  SF_TXZ_FACEMAP=$PWD/output/unroll_s$S/run \
  ./build/scroll_unroll output/run_placed output/unroll_s$S \
      --raw PHerc0139-4x5x5/cubes_RAW --steps $S --id run \
      --export-tifxyz output/unroll_s$S/tifxyz_run
done
```

Controls: add `--apply-ownership` for the ownership-applied run;
`--export-atlas output/unroll_atlas/atlas_run` (without
`--export-tifxyz`) for the per-wrap atlas; `--tifxyz-du 2 --tifxyz-dv 2`
for the resolution control.

## 4. Census, classification, event clustering (this repository)

```sh
uv run python bench/fiesta_gate.py unroll \
  --export  <scrollfiesta>/output/unroll_s12345/tifxyz_run \
  --facemap-prefix <scrollfiesta>/output/unroll_s12345/run \
  --workdir out/fiesta/unroll --out out/fiesta/unroll_s12345.json \
  --tag s12345
```

Census parameters are the corpus-frozen set (both diagonals,
exclude=1, cell=40, maxedge=60, touch_tol=1e-3). Events are 8-connected
region pairs (`bench/crossing_analyse.py`). The planted control is
`bench/fiesta_gate.py plant`.

## 5. The transaction

```sh
uv run python bench/fiesta_adapter.py \
  --export  <scrollfiesta>/output/unroll_s12345/tifxyz_run \
  --workdir out/fiesta/adapter \
  --out-dir out/fiesta/adapter/tifxyz_run_clean \
  --report  out/fiesta/adapter/transaction.json
```

The adapter refuses to hand over unless retained pixels are
byte-identical in every band and preserved sidecar and the shipped
output censuses clean under both triangulations.

## 6. Official reload

Build `vc_tifxyz_selfcross` from volume-cartographer
(merged: ScrollPrize/villa#1303) and run it on both the original export
(expected: clean flag `false`, counts matching step 4 exactly) and the
adapter output (expected: 0/0, clean flag `true`):

```sh
vc_tifxyz_selfcross <export-or-clean-dir> -o report.json
```

## Checking against the released numbers

Every table in the case-study README regenerates from steps 4-6. The
full per-contact records and the original/transformed exports for both
fixtures are attached to the case-study release as compressed
artifacts, each with a `.sha256` sidecar.
