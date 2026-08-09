#!/bin/bash
# benchmark_repro.sh — reproduce ONE benchmark transformation, end to end,
# in about a minute, and score the result the way a candidate is scored.
#
#   ./bench/benchmark_repro.sh                 # the default trace
#   SEGMENT=<name> ./bench/benchmark_repro.sh  # any pinned trace
#
# What it does, in order:
#
#   1. runs the frozen excision policy on the published input;
#   2. censuses the derivative it just produced with the OFFICIAL
#      vc_tifxyz_selfcross -- not our engine, because the only claim worth
#      making is that the artifact satisfies the invariant under the tool
#      the pipeline actually uses;
#   3. scores it with bench/topology_eval.py against the input, alongside
#      the reference derivative this project published.
#
# The default trace is deliberately one whose input IS the published
# original: no repair step stands between what you can download and what
# is cut, so the whole chain is visible in one run.
#
# Nothing here writes into out/excised: the run goes to its own directory.
set -euo pipefail

SEGMENT="${SEGMENT:-20260210000000-w058_2026021020}"
# NOT out/benchmark: that directory is the publishable package, and the
# official validator records the ABSOLUTE path of every surface it
# censuses. Local paths do not belong in a public artifact.
WORK="${WORK:-out/benchmark-evidence/repro/$SEGMENT}"
VALIDATOR="${VALIDATOR:-${VC_SELFCROSS:-vc_tifxyz_selfcross}}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO"

INDEX=out/release/index.json
python3 - "$SEGMENT" "$INDEX" <<'PY' > /tmp/.repro_paths.$$ || exit 1
import json, sys
seg, idx = sys.argv[1], sys.argv[2]
e = json.load(open(idx))["segments"].get(seg)
if e is None:
    sys.exit(f"{seg} is not in {idx}")
if e["disposition"] != "transformed":
    sys.exit(f"{seg} was not transformed ({e['disposition']}); pick a trace "
             "that actually required a cut")
print(e["input_mesh"]); print(e["output_mesh"]); print(e["certificate"])
PY
INPUT=$(sed -n 1p /tmp/.repro_paths.$$)
REFERENCE=$(sed -n 2p /tmp/.repro_paths.$$)
CERT=$(sed -n 3p /tmp/.repro_paths.$$)
rm -f /tmp/.repro_paths.$$

echo "segment   : $SEGMENT"
echo "input     : $INPUT"
echo "reference : $REFERENCE"
echo

rm -rf "$WORK"; mkdir -p "$WORK"

echo "== 1/3  running the frozen excision policy"
uv run python bench/excise_shadow.py --segment "$SEGMENT" --certificate \
  --out-root "$WORK" --base-manifest out/corpus_bases.json

CAND="$WORK/${SEGMENT}_excised.tifxyz"
[ -d "$CAND" ] || { echo "no derivative was emitted at $CAND" >&2; exit 1; }

echo
echo "== 2/3  censusing the derivative with the OFFICIAL validator"
REPORT="$WORK/census.json"
"$VALIDATOR" "$(cd "$CAND" && pwd)" -o "$(cd "$WORK" && pwd)/census.json"

echo
echo "== 3/3  scoring it as a candidate"
uv run python bench/topology_eval.py \
  --input "$INPUT" --candidate "$CAND" --reference "$REFERENCE" \
  --census-report "$REPORT" --json "$WORK/score.json"

echo
echo "== reproduction check: same cut as the published reference?"
uv run python - "$CAND" "$REFERENCE" "$CERT" <<'PY'
import json, sys
from pathlib import Path
sys.path.insert(0, "bench")
import numpy as np
from topology_eval import load_surface
cand, ref, cert = Path(sys.argv[1]), Path(sys.argv[2]), json.load(open(sys.argv[3]))
_, _, _, vc = load_surface(cand)
xr, yr, zr, vr = load_surface(ref)
same_cells = bool(np.array_equal(vc, vr))
print(f"retained cell set identical to the published reference: {same_cells}")
print(f"certificate says the cut removed "
      f"{cert['excision']['n_invalidated_vertices']} vertices, "
      f"{cert['excision']['n_removed_quads']} quads")
if not same_cells:
    print("\nA DIFFERENT cut is not a failure of the policy: the selection "
          "phase is budgeted and a solver may return a different optimal "
          "set of the same area. Compare the areas in score.json before "
          "concluding anything.")
PY
echo
echo "artifacts in $WORK"
