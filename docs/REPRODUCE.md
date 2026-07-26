# Reproducing every number in the report

This page exists so that a reviewer with a machine and no context can rerun the
audit and get the same answers, or find out precisely where they diverge. Every
figure quoted in `docs/submission.md` is produced by a command below.

Two conventions throughout:

- **Inputs are pinned by SHA-256.** Nothing under `data/` is committed. The
  manifests record the exact S3 key, byte count and hash of every file a result
  depends on, so an upstream file that is silently replaced shows up as a hash
  mismatch rather than as a mysteriously different number.
- **Runtimes are measured, not estimated.** The table below comes from
  `bench/time_pipeline.sh`, which writes `out/timings.txt` and a per-step log.
  Quote it from there; runtimes quoted from memory drift.

---

## 1. Machine the published timings were measured on

```
Apple M3 Pro, 12 cores, 36 GB
macOS 26.5.2
Apple clang 17.0.0 (clang-1700.6.3.2)
Python 3.13.2, numpy 2.5.1
```

Nothing requires a GPU. The core measurement needs no volume data at all — only
the published surface meshes.

---

## 2. Setup

```sh
git clone <repo> && cd windcheck
uv sync --extra viz

clang++ -O3 -std=c++17 -pthread -o engines/selfcross    engines/selfcross.cpp
clang++ -O3 -std=c++17 -pthread -o engines/atlas_query  engines/atlas_query.cpp

uv run pytest -q          # 21 tests, no data needed
```

The tests are the first checkpoint and need nothing downloaded. If they do not
pass, stop here — every later number depends on the predicate they pin, and two
of them exist because a determinism check found a real bug rather than a
reproducibility problem.

---

## 3. Data

Four published corpora. The S3 bucket is public and needs no credentials.

```sh
uv run python -m windcheck.fetch --sample all      # download + hash, 18.07 GB
uv run python -m windcheck.fetch --sample all --verify
```

`--verify` re-hashes every file against the manifest and exits nonzero on any
mismatch. To pin an already-downloaded tree without re-fetching, add
`--skip-download`.

| sample | local directory | manifest | files | size | role in the argument |
|---|---|---|---|---|---|
| `PHerc0172` | `data/scroll5_tifxyz` | `data/MANIFEST.json` | 424 | 0.64 GB | Scroll 5: 44 labelled segments + 9 multi-wrap partitions |
| `PHerc0814` | `data/PHerc0814_tifxyz` | `data/MANIFEST-PHerc0814.json` | 228 | 3.77 GB | second scroll; also the three-resolution sampling test |
| `PHerc0139` | `data/PHerc0139_tifxyz` | `data/MANIFEST-PHerc0139.json` | 756 | 7.10 GB | control corpus — its labelled windings cover 1.01–1.47 revolutions |
| `PHerc1667` | `data/PHerc1667_tifxyz` | `data/MANIFEST-PHerc1667.json` | 316 | 6.57 GB | control corpus — the one that broke the millimetre criterion |

Total 1,724 files, 18.07 GB.

The last two are the reason the headline is stated in revolutions rather than
millimetres. They are part of the evidence, not optional extras: see
`docs/submission.md` §4.1.

**A `.tifxyz` surface is a directory** — three float32 TIFFs plus `meta.json` —
not a file. Code that treats it as a file raises `IsADirectoryError`.

---

## 4. The pipeline

Run it all, with timings:

```sh
sh bench/time_pipeline.sh 2>&1 | tee out/timings.txt
```

Or step by step. Each line names the section of the report it produces.

```sh
# --- the census and the headline -----------------------------------------
uv run python bench/crossing_census.py \
      --root data/scroll5_tifxyz --volume 20241024131839 \
      --json out/crossing/census_v3.json --work out/crossing        # §3
uv run python bench/crossing_analyse.py                             # §4.1 table
uv run python bench/physical_report.py                              # §4.1 mm
uv run python bench/revolution_diag.py \
      --root data/scroll5_tifxyz --volume 20241024131839 \
      --dir out/crossing --json out/crossing/revdiag.json           # §4.1 rev
uv run python bench/revolution_summary.py                           # §4.1 pooled
uv run python bench/make_certificates.py                            # §6

# --- verification ---------------------------------------------------------
uv run python bench/crossing_quality.py                             # §4.4
uv run python bench/triangulability.py --sample 200                 # §4.3
uv run python bench/clustering_sensitivity.py                       # §7
uv run python bench/period_cross_check.py \
      --root data/PHerc1667_tifxyz --work out/period_1667 \
      --json out/crossing_1667/period.json                          # §4.1
uv run --with trimesh --with python-fcl --with rtree \
       python bench/validate_fcl.py --n 250                         # §5
```

### The four-corpus comparison

`revolution_summary.py` reads one `revdiag.json` per corpus, so each corpus has
to be censused and diagnosed first. Note the **separate `--work` directory per
corpus**: an earlier run wrote four scrolls' CSVs into one directory, and a
script that globbed it reported the mixture as a single corpus.

```sh
for c in "PHerc0814_tifxyz 20250804134230 0814" \
         "PHerc0139_tifxyz '' 0139" \
         "PHerc1667_tifxyz '' 1667"; do
  set -- $c
  uv run python bench/crossing_census.py --root data/$1 --volume "$2" \
        --json out/crossing_$3/census.json --work out/crossing_$3
  uv run python bench/revolution_diag.py --root data/$1 --volume "$2" \
        --dir out/crossing_$3 --json out/crossing_$3/revdiag.json
done
uv run python bench/revolution_summary.py
```

---

## 5. Measured runtimes

From `out/timings.txt` on the machine in §1. Per-step logs land in
`out/timing_<step>.log`.

```
step                            elapsed  status
------------------------------------------------
build_selfcross                    4.2s  exit 0
build_atlasquery                   0.6s  exit 0
tests                              3.9s  exit 0     21 tests
census_scroll5                    20.1s  exit 0     53 surfaces, both diagonals
analyse                            8.1s  exit 0
physical                           5.1s  exit 0
revolution_diag                    2.5s  exit 0     53 segments
revolution_summary                 0.3s  exit 0
certificates                       1.5s  exit 0
quality                            8.9s  exit 0
triangulability                    4.9s  exit 0     --sample 200
clustering                         6.9s  exit 0
period_crosscheck                 20.9s  exit 0     PHerc1667, 20 segments
fcl                                3.5s  exit 0     --n 250
------------------------------------------------
total                             91.4s
```

The whole audit reproduces in **under two minutes** once the data is on disk.
That is deliberate: a check nobody can afford to run is a check nobody runs.
Hashing all four corpora with `--verify` takes considerably longer than the
audit itself — 18.07 GB of SHA-256.

---

## 6. What each output file is

| file | contains |
|---|---|
| `out/crossing/census_v3.json` | per-surface triangle, pair and verdict counts, both quad diagonals, with the engine parameters that produced them |
| `out/crossing/*_d0.csv`, `*_d1.csv` | one row per intersecting triangle pair: the two quad indices, the verdict, penetration margin and crossing angle |
| `out/crossing/analysis.json` | pairs clustered into regions and events, setwise diagonal comparison |
| `out/crossing/physical.json` | per-trace areas and spans in mm |
| `out/crossing/revdiag.json` | covering span and separation in revolutions, per segment |
| `out/revolution_summary.json` | the pooled four-corpus table |
| `out/certificates/*_certificate.json` | the per-trace certificate |
| `out/certificates/*_points.json` | VC3D point collection, one point per event |

---

## 7. Expected results, and what a mismatch means

**Counts must be bit-identical across runs, thread counts and broad-phase cell
sizes.** 10,907 transverse and 1,989 grazing on `auto_grown_20251115002745`
across nine configurations. Two regression tests pin this. If your counts move
between runs of the same command, that is a bug and not a tolerance issue — it
was twice, and both times the determinism check found a real defect.

**The pooled four-corpus table** should read:

```
segments covering <= 1 revolution     n=56   max separation 0.115 rev
segments covering  > 2 revolutions    n=25   min nonzero    0.440 rev
```

**FCL agreement**: 250/250 on positives, 249/250 on negatives. The single
disagreement is expected and explained — quads (547,2256) and (547,2258) share
a vertex exactly, so FCL reports contact and we report no penetration. Both are
correct. If you see a *different* single disagreement, that is worth reporting.

**Timings** varying by 2–3× with hardware is normal. A step taking 100× longer
usually means the engine binary is missing and a Python fallback is running, or
that `--work` is pointing at a directory holding another corpus's CSVs.

---

## 8. Known limits of this provenance

Stated because a reviewer will find them.

- The certificate's `sha256_head` is a **partial** fingerprint: member names,
  sizes, and a bounded head of each file. Two surfaces differing only past the
  head collide. The manifests do hash full file contents; the certificate field
  is a cheap identifier, not a proof.
- `clustering_sensitivity.py` excludes files above 40,000 pairs because its
  product-space union is quadratic. It names every file it drops. Five of the
  six Scroll 5 exclusions are multi-wrap partitions, so that caveat is
  load-bearing and travels with the number.
- The revolution period is measured from the surface's own geometry and carries
  its own error. Two independent estimators agree to within about 15% where
  both are measurable; segments where they disagree are dropped, and the count
  of drops is printed.
- Runtimes were measured with the data already on local disk. First-run download
  time depends entirely on your link to S3.
