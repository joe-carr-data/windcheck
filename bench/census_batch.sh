#!/bin/bash
# census_batch.sh — census an arbitrary list of surfaces with the OFFICIAL
# validator, in ONE container.
#
# Generalises verify_corpus_official_batch.sh, which could only sweep a
# directory of *_excised.tifxyz. The benchmark's reference surfaces are not
# all in one directory: for a trace that needed no cut, the reference IS
# the published input.
#
# ONE container is not a style choice. The per-surface version spawned a
# container per surface and wedged the Docker daemon at 3 of 154, with no
# way to distinguish "slow" from "stuck".
#
# MANIFEST is a TSV of  <name>\t<surface-dir>  with paths that exist at the
# SAME absolute path inside the container, which is why the mount point and
# the host path are identical. A report written anywhere else -- macOS
# $TMPDIR, say -- succeeds and then lands nowhere the caller can read it.
set -uo pipefail

MANIFEST="${MANIFEST:?set MANIFEST (TSV: name<TAB>surface-dir)}"
REPORTS="${REPORTS:?set REPORTS (output directory)}"
IMAGE="${IMAGE:-vc-builder:latest}"
BIN="${BIN:-/villa/build/bin/vc_tifxyz_selfcross}"
VILLA="${VILLA:?set VILLA to a volume-cartographer checkout with build/bin}"
MOUNT="${MOUNT:-$PWD}"

[ -f "$MANIFEST" ] || { echo "manifest not found: $MANIFEST" >&2; exit 2; }
case "$MANIFEST" in "$MOUNT"/*) ;; *)
  echo "REFUSING: manifest $MANIFEST is outside the mount $MOUNT, so the" >&2
  echo "container cannot read it." >&2; exit 2 ;;
esac
case "$REPORTS" in "$MOUNT"/*) ;; *)
  echo "REFUSING: reports dir $REPORTS is outside the mount $MOUNT, so the" >&2
  echo "reports would land nowhere the caller can read." >&2; exit 2 ;;
esac
# A non-empty reports directory is refused. Reusing one silently mixes
# fresh reports with stale ones from an earlier operand, and every
# downstream check would then be reading a census of something else.
if [ -d "$REPORTS" ] && [ -n "$(ls -A "$REPORTS" 2>/dev/null)" ]; then
  echo "REFUSING: $REPORTS is not empty. Reports from a previous run would" >&2
  echo "be indistinguishable from this one's. Remove it or pick another." >&2
  exit 2
fi
mkdir -p "$REPORTS"

# The image tag is recorded by DIGEST as well as name: `latest` moves, and
# a report whose validator cannot be identified binds nothing.
IMAGE_DIGEST=$(docker image inspect "$IMAGE" --format '{{.Id}}' 2>/dev/null || echo "UNKNOWN")

set +e
docker run --rm --platform linux/arm64 --user "$(id -u):$(id -g)" \
  -v "$VILLA:/villa" -v "$MOUNT:$MOUNT" \
  -e MANIFEST="$MANIFEST" -e REPORTS="$REPORTS" -e BIN="$BIN" \
  --entrypoint /bin/bash "$IMAGE" -c '
set -uo pipefail
n=0; fail=0
while IFS="$(printf "\t")" read -r name dir; do
  [ -n "${name:-}" ] || continue
  n=$((n+1))
  if [ ! -d "$dir" ]; then
    fail=$((fail+1)); printf "%4d %-52s NO SUCH DIRECTORY\n" "$n" "$name"
    continue
  fi
  if "$BIN" "$dir" -o "$REPORTS/$name.json" > "$REPORTS/$name.stdout" 2>&1; then
    tri=$(grep -o "diagonal 0: [0-9]* triangles" "$REPORTS/$name.stdout" \
          | head -1 | awk "{print \$3}")
    printf "%4d %-52s triangles=%s\n" "$n" "$name" "${tri:-?}"
  else
    fail=$((fail+1)); printf "%4d %-52s CENSUS FAILED\n" "$n" "$name"
  fi
done < "$MANIFEST"
echo "done: $n censused, $fail failed"
[ "$fail" -eq 0 ]
'
RC=$?
set -e

# A sweep that "records failures but ends successfully" is a sweep whose
# caller will act on partial results.
if [ "$RC" -ne 0 ]; then
  echo "CENSUS SWEEP FAILED (exit $RC): at least one surface did not census" >&2
  exit "$RC"
fi

# Bind each report to the operand bytes it describes and to the validator
# that produced it, so a report cannot later be paired with a different
# surface of the same shape.
python3 "$(dirname "$0")/bind_census_reports.py" \
  --manifest "$MANIFEST" --reports "$REPORTS" \
  --image "$IMAGE" --image-digest "$IMAGE_DIGEST" --bin "$BIN" \
  --villa "$VILLA"
