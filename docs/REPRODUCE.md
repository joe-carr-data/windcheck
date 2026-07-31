# Reproducing the results

Everything published here is reproducible from the `tifxyz` surfaces alone. No
volume data, no GPU and no model weights are needed for any command on this
page.

There are four levels, in increasing cost:

| level | what it establishes | cost |
|---|---|---|
| 1. Spot reproduction | one segment, end to end, from published input to certified output | seconds |
| 2. Independent recensus | that an emitted artifact really is clean, without trusting this tool | seconds |
| 3. Corpus verification | that all 185 records are well-formed, hash-consistent and clean | ~20 min |
| 4. Full corpus pass | regenerating every output from scratch | ~50 min |

---

## 0. Setup

```sh
uv sync
clang++ -O3 -std=c++17 -pthread -o engines/selfcross engines/selfcross.cpp
uv run pytest -q
```

The compile takes about 4 s. The test suite needs no data.

Timings on this page were measured on an Apple-silicon laptop, single machine,
no cluster. The census kernel is multi-threaded; the default thread count is
chosen per command and recorded in every certificate, so a timing difference
never changes a count.

### Data

Segment surfaces come from the open data buckets:

```sh
uv run python -m windcheck.fetch --help
```

The fetcher writes a SHA-256 manifest next to what it downloads, and
`data/MANIFEST*.json` is committed. Every base mesh used in the corpus pass is
pinned by hash in `out/corpus_bases.json`, so a reproduction that starts from
different bytes is detected rather than silently accepted.

---

## 1. Spot reproduction of a single segment

One command, one segment, no drivers:

```sh
uv run windcheck transform data/scroll1_tifxyz/20231005123336 \
    --out out/spot/20231005123336
```

Expected, in about 2.4 s:

```
out/spot/20231005123336/
  20231005123336_transformed.tifxyz/   the emitted surface (x/y/z planes + mask)
  20231005123336_transform_certificate.json
  work_<hash>/                         scratch, safe to delete
```

The printed summary should read 20 transverse contacts under each diagonal
before, 0 and 0 after, 6 removed quads, 99.9992% retained area, `STATUS clean`. `transform` exits 0 only if the mesh it just wrote,
read back from disk, censuses at 0/0 under both canonical triangulations. A
non-zero exit means the claim failed, not that the tool crashed.

The second worked example exercises the other operator — bounded displacement
rather than excision:

```sh
uv run windcheck transform \
    data/scroll5_tifxyz/20251205115859-w094_20251205115859448_flatboi \
    --volume 20241024131839 --out out/spot/w094
```

About 0.25 s. This segment publishes two meshes, one per volume, so `--volume`
is required to disambiguate. The certificate records
`displacement.applied: true` and a terminal disposition of `already_clean`:
after the bounded displacement base is taken, nothing needs cutting.

To see the same segment before any transformation, use the report-only command,
which never writes geometry:

```sh
uv run windcheck check \
    data/scroll5_tifxyz/20251205115859-w094_20251205115859448_flatboi \
    --volume 20241024131839 --out out/spot/w094-check
```

About 0.2 s, and it should print a **not clean** verdict with 4 transverse
contacts under diagonal 0 and 7 under diagonal 1. It writes
`<segment>_check_certificate.json` and `<segment>_points.json`, the latter being
a viewer-loadable `PointCollection` placing one marker at the deepest contact of
each crossing event.

To force the excision path on the published mesh instead of the displaced base,
add `--no-displacement` (about 0.8 s, 4/7 contacts to 0/0, 6 quads removed,
99.9989% area retained).

---

## 2. Recensusing an emitted artifact independently

This is the check that matters, because it does not have to trust anything
above. Run the census engine directly on an emitted mesh:

```sh
uv run windcheck check out/spot/20231005123336/20231005123336_transformed.tifxyz \
    --out out/spot/recensus
```

Expected, in about 0.6 s: `clean`, 0 transverse contacts under diagonal 0 and 0
under diagonal 1, and 0 crossing events. The certificate this writes records the engine binary's own SHA-256 and every
census parameter, so a reader can confirm which kernel produced the verdict.

Two properties make this a real check rather than a restatement:

- the census runs on bytes read back from disk, in a fresh working directory,
  with no reuse of any intermediate the transform produced;
- it runs under **both** canonical quad triangulations. A cut that only
  resolves contacts under one diagonal fails here.

If you would rather not use this tool at all, the emitted `tifxyz` is a standard
surface: load it in your own viewer or intersection checker. Coordinates of
retained cells are bit-identical to the input, and excised cells are marked
invalid, so any independent tool sees exactly the same geometry.

---

## 3. Verifying the whole corpus

This re-derives every claim in the release from the artifacts on disk, without
trusting the driver that produced them: certificate well-formedness, roster
coverage in both directions, re-hashing of every emitted mesh, and a fresh
both-diagonal recensus of every one of them in a throwaway workdir.

```sh
uv run python bench/verify_corpus.py --jobs 4 --json out/verification.json
```

Expected result, roughly 20 minutes:

- 185/185 certificates well-formed, with base hashes still matching
  `out/corpus_bases.json`;
- 154/154 emitted meshes re-hash to their recorded output hashes;
- 179/179 recensus at 0 transverse contacts on both diagonals;
- 0 disagreements against the recorded census on any field;
- the 25 already-clean segments re-measured from their base mesh, also 0/0.

Then recompute the published headline and the corpus index from those artifacts:

```sh
uv run python bench/headline_decision.py
uv run python bench/build_release_index.py
```

`headline_decision.py` evaluates the eight pre-registered conditions and prints
the sentence with the numbers filled in from the data. It is deliberately not
editable prose: the only way to change what it prints is to change the corpus.
Expect six of eight conditions to pass, a corpus retention of 99.505%, and the
area and fragmentation qualifications printed separately below the sentence, not
folded into it. Exit status is the gate — 0 for the strong sentence, 1 for the
fallback.

`build_release_index.py` regenerates `out/release/index.json` and
`docs/CORPUS.md` from the same certificates. Both should be byte-stable apart
from the generation timestamp.

---

## 4. Regenerating the corpus from scratch

```sh
uv run python bench/excise_corpus.py --jobs 4
```

About 50 minutes on four workers. The base manifest, not a directory glob,
decides which bytes are processed, and a segment is skipped only when a valid
certificate already exists whose base hashes, policy hash and source-tree
digest all match. Use `--force` to ignore checkpoints.

Per-segment runtimes to expect: minimum 0.07 s, median 36.8 s, 90th percentile
140 s, maximum 211 s, against a hard 600 s per-segment process limit. No segment
in the published pass came within a factor of two of that limit.

Smaller slices are useful for a smoke test:

```sh
uv run python bench/excise_corpus.py --corpus "Scroll 1" --limit 5
uv run python bench/verify_corpus.py --limit 3
```

---

## 5. What each output file is

| path | contents |
|---|---|
| `out/corpus_bases.json` | the base manifest: per segment, the pinned input bytes, their content manifest and directory digest, and whether the base is the original mesh or a displacement-repaired one |
| `out/excised/corpus/<segment>_excision_certificate.json` | the terminal record for one segment |
| `out/excised/corpus/corpus_summary.jsonl` | one line per segment, as the driver saw it |
| `out/verification.json` | the independent verification record from level 3 |
| `out/release/index.json` | the machine-readable release index: input hash to output hash, disposition, certificate path, both area denominators, fragmentation metrics, runtime |
| `docs/CORPUS.md` | the same index, human-readable, including the six non-censusable inputs |

---

## 6. Verifying provenance without a repository

Nothing published cites a commit sha, because a reader of a release cannot
resolve one. Certificates and `out/release/index.json` instead record the code
version, the frozen policy version and hash, and a source-tree digest over
every published file that can change a result. Recompute it from the release:

```sh
uv run python -m windcheck.provenance
uv run python -m windcheck.provenance --verify <source_tree_digest>
uv run python -m windcheck.provenance --manifest
```

The serialisation is frozen and stated in `src/windcheck/manifest.py`, and it
is the same rule the per-mesh content manifests use, so an independent
implementation can check every digest in the release.

---

## 7. When a number does not match

Contact counts come from a deterministic floating-point triangle-intersection
predicate with scale-aware tolerances, cross-validated against independent
implementations. The predicate is deterministic, not exact: the same binary on
the same bytes reproduces the same counts, and a different build may not.
If counts do not match, check in this order:

1. **The engine binary.** Every certificate records the engine's source and
   binary SHA-256. A different compiler or optimisation level produces a
   different binary hash; the counts should still agree, and if they do not,
   that is a genuine finding worth reporting.
2. **The input bytes.** Compare against `original_hashes` in
   `out/corpus_bases.json`. Published meshes have been revised before.
3. **The diagonal.** A count taken under one triangulation is not comparable to
   a count taken under both.
4. **Adjacency.** The census excludes contacts between grid-adjacent cells by
   construction; a tool that does not will report far more.

Wall-clock timings will differ from this page and are not part of any claim,
except for the 600 s per-segment limit, which is a policy bound and is recorded
per segment.

---

## 8. Known limits of this reproduction

- The six triangle-empty or invalid inputs cannot be censused at all. They have
  terminal records, not verdicts, and no level of this page will produce a
  cleanliness claim for them.
- Bounded displacement was applied before excision on 103 of the 185 segments.
  Reproducing those bases from the original published meshes is a separate and
  longer computation than the excision pass; the corpus pass starts from the
  hash-verified repaired meshes recorded in the base manifest.
- The census certifies non-adjacent transverse contacts. Coplanar overlap and
  grazing contact are outside its scope, so a clean verdict is not a claim that
  no two parts of the surface touch at all.
