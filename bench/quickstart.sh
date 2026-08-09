#!/bin/bash
# quickstart.sh — see the benchmark work, on one small trace, in a minute.
#
#   ./tools/quickstart.sh                      # fetches the input for you
#   INPUT=/path/to/trace.tifxyz ./tools/quickstart.sh   # if you already have it
#
# It runs the two demonstrations that define the evaluation contract:
#
#   A. the ORIGINAL trace, scored as if it were a candidate
#      -> rejected: it still self-intersects, so it is NOT scored on cost
#   B. the REFERENCE derivative
#      -> clean, with its full retention / fragmentation / fidelity vector
#
# and then shows exactly where you substitute your own derivative.
#
# The reference is not shipped and is not downloaded: it is RECONSTRUCTED
# from the input and this package's own `masks.npz`, which is also the
# fastest way to see that the reconstruction claim in the README is true.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
PKG="$(cd "$HERE/.." && pwd)"
SEGMENT="${SEGMENT:-20260126000000-w045_2026012619}"
TRACE="$PKG/traces/$SEGMENT"
WORK="${WORK:-$PKG/quickstart-work}"
VALIDATOR="${VALIDATOR:-${VC_SELFCROSS:-vc_tifxyz_selfcross}}"
PY="${PYTHON:-python3}"

[ -d "$TRACE" ] || { echo "no such trace in this package: $SEGMENT" >&2; exit 2; }
command -v "$VALIDATOR" >/dev/null 2>&1 || [ -x "$VALIDATOR" ] || {
  echo "REFUSING: cannot find the validator '$VALIDATOR'." >&2
  echo "It is volume-cartographer's vc_tifxyz_selfcross. Set VC_SELFCROSS" >&2
  echo "to it, or to a launcher script that runs it." >&2
  exit 2; }

# The reconstruction below and the evaluator itself need three packages,
# and $PY defaults to the system python3, which on a fresh machine has
# none of them. Without this the first command in the README dies on a
# raw `ModuleNotFoundError` traceback partway through -- after the S3
# fetch has already run. Say what is missing, before spending anything.
"$PY" - <<'PY' || exit 2
import sys
# IMPORTED, not merely located. importlib.util.find_spec() finds a
# package whose compiled extension is broken or built for another
# architecture, and the failure then surfaces mid-run as a traceback.
need, missing = ("numpy", "tifffile", "scipy"), []
for m in need:
    try:
        __import__(m)
    except Exception:                                         # noqa: BLE001
        missing.append(m)
if missing:
    print(f"REFUSING: {sys.executable} cannot import: {', '.join(missing)}",
          file=sys.stderr)
    print(f"  pip install {' '.join(missing)}", file=sys.stderr)
    print("  ...or re-run with PYTHON=<an interpreter that has them>.",
          file=sys.stderr)
    sys.exit(2)
PY

echo "=== segment: $SEGMENT"
mkdir -p "$WORK"

# ---- 1. the input -----------------------------------------------------
if [ -n "${INPUT:-}" ]; then
  echo "using the input you supplied: $INPUT"
else
  INPUT="$WORK/input.tifxyz"
  if [ ! -d "$INPUT" ]; then
    URI=$("$PY" - "$PKG" "$SEGMENT" <<'PY'
import json, sys
loc = json.load(open(sys.argv[1] + "/locators.json"))["traces"][sys.argv[2]]
up = (loc.get("input") or {}).get("upstream") or {}
print(up.get("uri") or "")
PY
)
    [ -n "$URI" ] || { echo "no upstream locator for $SEGMENT; set INPUT=..." >&2; exit 2; }
    command -v aws >/dev/null 2>&1 || {
      echo "REFUSING: the AWS CLI is not installed, so I cannot fetch" >&2
      echo "  $URI" >&2
      echo "Fetch it yourself and re-run with INPUT=<that directory>." >&2
      exit 2; }
    echo "fetching the published input (a few MB):"
    echo "  aws s3 sync --no-sign-request $URI $INPUT/"
    aws s3 sync --no-sign-request "$URI" "$INPUT/" >/dev/null
  fi
fi

# ---- 2. verify it is the right bytes, then rebuild the reference -------
"$PY" - "$TRACE" "$INPUT" "$WORK/reference.tifxyz" <<'PY'
import hashlib, json, shutil, sys
from pathlib import Path
import numpy as np, tifffile

trace, inp, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
prov = json.loads((trace / "provenance.json").read_text())

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

want = (prov.get("input") or {}).get("hashes") or {}
bad = [a for a in "xyz" if want.get(a) and sha(inp / f"{a}.tif") != want[a]]
if bad:
    sys.exit(f"REFUSING: the input does not match provenance.json on {bad}. "
             "This is not the trace the benchmark measured.")
print(f"input verified by hash against provenance.json  ({', '.join('xyz')} ok)")

removed = np.load(trace / "masks.npz")["removed"].astype(bool)
if out.exists():
    shutil.rmtree(out)
out.mkdir(parents=True)
for a in "xyz":
    arr = np.asarray(tifffile.imread(inp / f"{a}.tif"))
    tifffile.imwrite(out / f"{a}.tif", np.where(removed, arr.dtype.type(-1), arr))
tifffile.imwrite(out / "mask.tif",
                 np.where(removed, np.uint8(0), np.uint8(255)))
for extra in ("meta.json",):
    if (inp / extra).is_file():
        shutil.copy(inp / extra, out / extra)
print(f"reference rebuilt from the input and masks.npz: "
      f"{int(removed.sum())} cells removed  -> {out}")
PY

# ---- the gate --------------------------------------------------------
# Both demonstrations are CHECKED, not narrated. The reasoning, and every
# condition, is in tools/quickstart_gate.py -- it lives in a module
# rather than here so the failure modes can be regression-tested by
# feeding it reports, instead of by arranging for a validator to break.
gate() {  # gate <dirty|clean> <report> <input> <candidate> <exit-status> [extra...]
  local expect="$1" report="$2" in="$3" cand="$4" rc="$5"; shift 5
  "$PY" "$HERE/quickstart_gate.py" --expect "$expect" --report "$report" \
    --input "$in" --candidate "$cand" --exit-status "$rc" \
    --package "$PKG" --workdir "$WORK" "$@"
}

echo
echo "=== A. the ORIGINAL trace, scored as a candidate"
echo "    (expected: rejected -- it still self-intersects)"
# Removed before the run, never after: a report left by a previous
# invocation is otherwise read as this one's result.
rm -f "$WORK/original-score.json"
set +e
"$PY" "$HERE/topology_eval.py" --input "$INPUT" --candidate "$INPUT" \
  --provenance "$TRACE/provenance.json" --validator "$VALIDATOR" \
  --json "$WORK/original-score.json"
RC=$?
set -e
gate dirty "$WORK/original-score.json" "$INPUT" "$INPUT" "$RC" || exit 1

echo
echo "=== B. the REFERENCE derivative, scored the same way"
echo "    (expected: clean, with its cost vector)"
rm -f "$WORK/reference-score.json"
set +e
"$PY" "$HERE/topology_eval.py" --input "$INPUT" \
  --candidate "$WORK/reference.tifxyz" \
  --provenance "$TRACE/provenance.json" --validator "$VALIDATOR" \
  --json "$WORK/reference-score.json"
RC=$?
set -e
# The extra bindings make this "the packaged reference", not merely
# "some clean surface": the reconstruction must hash to the reference
# derivative provenance.json names, and its cost vector must match the
# certificate figures shipped in fields.json.
gate clean "$WORK/reference-score.json" "$INPUT" "$WORK/reference.tifxyz" \
  "$RC" --provenance "$TRACE/provenance.json" \
  --fields "$TRACE/fields.json" || exit 1

cat <<EOF

=== now score YOUR derivative

  $PY tools/topology_eval.py \\
      --input      $INPUT \\
      --candidate  <YOUR derivative for this trace> \\
      --provenance $TRACE/provenance.json \\
      --validator  $VALIDATOR

Your candidate must be a tifxyz on the SAME grid as the input, with
mask.tif using 255 for retained cells. Beating the reference means
staying clean while retaining more area -- the reference is one
admissible answer, not the target.

Where the defects are, for this trace: traces/$SEGMENT/witness.json
EOF
