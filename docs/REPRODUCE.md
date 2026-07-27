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

Five published corpora. The S3 bucket is public and needs no credentials.

```sh
uv run python -m windcheck.fetch --sample all           # 740 files, 1.55 GB
uv run python -m windcheck.fetch --sample all --verify
```

`--verify` re-hashes every file against its manifest and exits nonzero on any
mismatch. To pin an already-downloaded tree without re-fetching, add
`--skip-download`.

| sample | local directory | volume the audit reads | files | size |
|---|---|---|---|---|
| `PHercParis4` | `data/scroll1_tifxyz` | `20230205180739` (7.91 µm) | 220 | 0.94 GB |
| `PHerc0172` | `data/scroll5_tifxyz` | `20241024131839` (7.91 µm) | 212 | 0.41 GB |
| `PHerc1667` | `data/PHerc1667_tifxyz` | `20231117161658` (7.91 µm) | 80 | 0.09 GB |
| `PHerc0139` | `data/PHerc0139_tifxyz` | `20250728140407` (9.362 µm) | 152 | 0.06 GB |
| `PHerc0814` | `data/PHerc0814_tifxyz` | `20250804134230` (9.362 µm) | 76 | 0.04 GB |

**Each corpus is pinned to one volume, and that is not a convenience.** These
scrolls are published at up to four resolutions. The crossing CSVs hold `(v, u)`
indices into the grid the census measured, so reading a different resolution of
the same segment silently reinterprets them against the wrong grid — which
happened here across all 20 PHerc1667 segments and produced a plausible,
entirely fictional cluster of results. The manifests therefore record exactly
the files a result depends on, and the fetcher downloads exactly that set:
1.55 GB rather than the 18 GB an unfiltered pull would bring.

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

### The five-corpus comparison

`revolution_summary.py` reads one `revdiag.json` per corpus, so each corpus has
to be censused and diagnosed first. Note the **separate `--work` directory** and
the **explicit `--volume`** on every line: an earlier run wrote four scrolls'
CSVs into one directory, and a later one read the wrong resolution back. Both
produced wrong numbers that looked entirely reasonable.

```sh
run() {   # corpus  volume  tag
  uv run python bench/crossing_census.py --root "data/$1" --volume "$2" \
        --json "out/crossing_$3/census.json" --work "out/crossing_$3"
  uv run python bench/revolution_diag.py --root "data/$1" --volume "$2" \
        --dir "out/crossing_$3" --json "out/crossing_$3/revdiag.json"
  uv run python bench/period_cross_check.py --root "data/$1" --volume "$2" \
        --work "out/period_$3" --json "out/crossing_$3/period.json"
}
run scroll1_tifxyz   20230205180739 s1
run scroll5_tifxyz   20241024131839 s5
run PHerc0139_tifxyz 20250728140407 0139
run PHerc0814_tifxyz 20250804134230 0814
run PHerc1667_tifxyz 20231117161658 1667

uv run python bench/revolution_summary.py
```

`revolution_diag.py` refuses to run when the mesh it reads has a different grid
shape from the census beside it, so a wrong `--volume` fails loudly rather than
producing a plausible answer.

Every segment, with a certificate and overlay each:

```sh
uv run python bench/precompute_all.py --out results
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

**Scope of that 91.4 s: it is the Scroll 5 pipeline only** — one corpus of 53
surfaces, plus the verification steps. Censusing all five corpora, 179 segments
and 228 M triangles, takes **82.7 s** on the same machine; `precompute_all.py`,
which additionally computes a second period estimate and writes a certificate
and overlay per segment, takes several minutes.

The point stands either way: a check nobody can afford to run is a check nobody
runs. Hashing all five corpora with `--verify` takes longer than the audit
itself.

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

**The pooled five-corpus strata** should read:

```
  no exclusion at all           <=1-period n=66  max 0.417   >2-period n=59  min nz 0.823
  period agreement enforced     <=1-period n=38  max 0.039   >2-period n=53  min nz 0.823
```

and `precompute_all.py` should classify 179 segments as: crossing present 160 /
none 19; period agreed 96 / disagreed 38 / unavailable 45.

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
- Each corpus is pinned to a single published volume. Results at other
  resolutions of the same segments are not covered by these manifests, and the
  analysis will refuse to mix them.
