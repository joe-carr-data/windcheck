#!/bin/sh
# Measure every step of the published pipeline, end to end, on this machine.
#
# The provenance page quotes runtimes so a reviewer can tell a slow step from a
# hung one. Quoting them from memory is how they drift, so they come from here.
# Output is one line per step: elapsed seconds, exit status, command.
#
#   sh bench/time_pipeline.sh 2>&1 | tee out/timings.txt
set -u
cd "$(dirname "$0")/.." || exit 1

run() {
    label="$1"; shift
    start=$(python3 -c 'import time; print(time.time())')
    "$@" > "out/timing_$label.log" 2>&1
    status=$?
    end=$(python3 -c 'import time; print(time.time())')
    python3 -c "print(f'{'"'"'$label'"'"':28s} {$end - $start:8.1f}s  exit $status')"
}

mkdir -p out
echo "step                            elapsed  status"
echo "------------------------------------------------"

run build_selfcross   clang++ -O3 -std=c++17 -pthread -o engines/selfcross engines/selfcross.cpp
run build_atlasquery  clang++ -O3 -std=c++17 -pthread -o engines/atlas_query engines/atlas_query.cpp
run tests             uv run pytest -q

run census_scroll5    uv run python bench/crossing_census.py \
                          --root data/scroll5_tifxyz --volume 20241024131839 \
                          --json out/crossing/census_v3.json --work out/crossing
run analyse           uv run python bench/crossing_analyse.py
run physical          uv run python bench/physical_report.py
run revolution_diag   uv run python bench/revolution_diag.py \
                          --root data/scroll5_tifxyz --volume 20241024131839 \
                          --dir out/crossing --json out/crossing/revdiag.json
run revolution_summary uv run python bench/revolution_summary.py
run certificates      uv run python bench/make_certificates.py

run quality           uv run python bench/crossing_quality.py
run triangulability   uv run python bench/triangulability.py --sample 200
run clustering        uv run python bench/clustering_sensitivity.py
run period_crosscheck uv run python bench/period_cross_check.py \
                          --root data/PHerc1667_tifxyz --work out/period_1667 \
                          --json out/crossing_1667/period.json
run fcl               uv run --with trimesh --with python-fcl --with rtree \
                          python bench/validate_fcl.py --n 250

echo "------------------------------------------------"
echo "per-step logs in out/timing_<step>.log"
