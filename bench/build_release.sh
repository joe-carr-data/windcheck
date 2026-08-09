#!/bin/bash
# Assemble the topology-benchmark release archive, reproducibly.
#
#   VERSION=v0.3.2 ./bench/build_release.sh
#
# The tools directory used to be assembled by hand. That is how v0.3.0
# shipped tools built from UNCOMMITTED working-tree files -- caught in
# review, but only just. Every tool here is therefore taken from the
# committed tree with `git show`, and then verified byte-identical
# against the working copy, so the archive can only contain code that is
# in the repository at the recorded commit.
#
# It also does the three things macOS will otherwise break:
# deletes AppleDouble members, tars with COPYFILE_DISABLE, and audits the
# resulting member list. The first v0.3.0 candidate had 3,148 members,
# exactly half of them `._*` metadata.
set -euo pipefail

VERSION="${VERSION:?set VERSION, e.g. VERSION=v0.3.2}"
SRC="${SRC:-out/benchmark}"
STAGE="${STAGE:-out/release-stage}"
NAME="windcheck-${VERSION}-topology-benchmark"
REF="${REF:-HEAD}"

TOOLS=(
  topology_eval.py
  quickstart.sh
  quickstart_gate.py
  benchmark_card.py
  benchmark_accept.py
  benchmark_mutations.py
  benchmark_selfcheck.py
  benchmark_repro.sh
  bind_census_reports.py
  bind_witness.py
  build_locators.py
  census_batch.sh
  audit_archive.py
  pack_benchmark.py
)

# Payload files that must be present, or the archive is incomplete in a
# way a member count would not reveal.
REQUIRED=(
  index.json locators.json README.md CARD.md reconciliation.json
  B7-acceptance.json selfcheck.json witness-binding.json
  negative-controls.json
)

echo "=== $NAME from $SRC at $(git rev-parse --short "$REF")"

[ -d "$SRC" ] || { echo "no such payload directory: $SRC" >&2; exit 2; }
if [ -n "$(git status --porcelain -- bench/)" ]; then
  echo "REFUSING: bench/ has uncommitted changes. The archive records a" >&2
  echo "commit; shipping tools that are not in it is exactly the defect" >&2
  echo "this script exists to prevent." >&2
  git status --porcelain -- bench/ >&2
  exit 2
fi

rm -rf "$STAGE"
mkdir -p "$STAGE/$NAME"
# -a would carry macOS extended attributes into the archive.
cp -R "$SRC"/. "$STAGE/$NAME"/
mkdir -p "$STAGE/$NAME/tools"

for f in "${TOOLS[@]}"; do
  git show "$REF:bench/$f" > "$STAGE/$NAME/tools/$f" \
    || { echo "REFUSING: bench/$f is not in $REF" >&2; exit 2; }
  if ! cmp -s "$STAGE/$NAME/tools/$f" "bench/$f"; then
    echo "REFUSING: bench/$f differs between $REF and the working tree" >&2
    exit 2
  fi
done
chmod +x "$STAGE/$NAME/tools/"*.sh

for f in "${REQUIRED[@]}"; do
  [ -f "$STAGE/$NAME/$f" ] || {
    echo "REFUSING: the payload has no $f" >&2; exit 2; }
done

find "$STAGE/$NAME" \( -name '._*' -o -name '.DS_Store' \) -delete

# Byte-reproducible, not merely content-reproducible. Two builds of the
# same commit previously produced archives that were file-for-file
# identical but 33 bytes apart, because `cp -R` stamps fresh mtimes into
# the tar headers and gzip writes its own timestamp. A benchmark whose
# pitch is auditability should not have a digest that changes when
# nothing did: a reviewer rebuilding it must get the same sha256, or the
# hash proves nothing they can check.
#
# The fixed timestamp is the COMMIT's author date, so it is a property of
# the code being shipped rather than of the machine that shipped it.
STAMP="$(git show -s --format=%ct "$REF")"
TZ=UTC find "$STAGE/$NAME" -exec touch -h -t \
  "$(TZ=UTC date -r "$STAMP" +%Y%m%d%H%M.%S)" {} +
COPYFILE_DISABLE=1 tar --uid 0 --gid 0 --numeric-owner \
  -cf "$STAGE/$NAME.tar" -C "$STAGE" "$NAME"
gzip -n -9 -c "$STAGE/$NAME.tar" > "$STAGE/$NAME.tar.gz"
rm -f "$STAGE/$NAME.tar"

uv run python bench/audit_archive.py --archive "$STAGE/$NAME.tar.gz" \
  --expect-top "$NAME"

echo
echo "sha256  $(shasum -a 256 "$STAGE/$NAME.tar.gz" | cut -d' ' -f1)"
echo "bytes   $(wc -c < "$STAGE/$NAME.tar.gz" | tr -d ' ')"
echo "members $(tar tzf "$STAGE/$NAME.tar.gz" | wc -l | tr -d ' ')"
echo "tools   $(tar tzf "$STAGE/$NAME.tar.gz" | grep -c '/tools/[^/]*$')"
echo "commit  $(git rev-parse "$REF")"
echo "archive $STAGE/$NAME.tar.gz"
