#!/bin/sh
# Clean-checkout end-to-end smoke for the topology transaction.
#
# Proves the one command works from a pristine copy of this repository:
# no repo-local state, no pre-built engine, no ambient environment. It
# exports HEAD to a temporary directory, builds the census engine there,
# creates the environment with uv, synthesizes a clean segment, and runs
# a full transaction, asserting the committed output and its in-output
# authoritative report.
#
# Requires: git, uv, and a C++17 compiler (clang++ or c++).
set -eu

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

git archive HEAD | tar -x -C "$tmp/"
cd "$tmp"

if command -v clang++ >/dev/null 2>&1; then CXX=clang++; else CXX=c++; fi
"$CXX" -O3 -std=c++17 -pthread -o engines/selfcross engines/selfcross.cpp
uv sync --frozen --quiet

uv run python - <<'EOF'
import json
import numpy as np
import tifffile
from pathlib import Path

d = Path("smoke_seg")
d.mkdir()
H, W = 80, 80
x, y = np.meshgrid(np.arange(W, dtype=np.float32),
                   np.arange(H, dtype=np.float32))
tifffile.imwrite(d / "x.tif", x)
tifffile.imwrite(d / "y.tif", y)
tifffile.imwrite(d / "z.tif", np.full((H, W), 5.0, np.float32))
(d / "meta.json").write_text(json.dumps(
    {"scale": [1.0, 1.0], "uuid": "smoke", "type": "seg",
     "format": "tifxyz"}))
EOF

uv run windcheck transaction smoke_seg --out smoke_out \
    --report smoke_report.json
test -d smoke_out
test -f smoke_out/windcheck_transaction/certificate.json
test -f smoke_out/windcheck_transaction/transaction_report.json

uv run python - <<'EOF'
import json
rep = json.load(open("smoke_report.json"))
assert rep["exit_code"] == 0, rep
assert rep["committed"] is True, rep
inner = json.load(open(
    "smoke_out/windcheck_transaction/transaction_report.json"))
# the post-commit copy must equal the committed authoritative report
assert inner == rep, (inner, rep)
cert = json.load(open(
    "smoke_out/windcheck_transaction/certificate.json"))
# provenance is the recomputable source-tree digest, not a commit SHA
assert cert.get("source_tree_digest"), "certificate lacks a tree digest"
print("clean-checkout smoke OK: committed, exit 0, reports identical, "
      "certificate carries a source-tree digest")
EOF
