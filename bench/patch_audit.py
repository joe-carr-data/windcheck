"""Census the published `verified_patches` set for self-intersection.

The patches are a different population from the pinned corpus: tens of
thousands of small surfaces rather than 185 large ones, published on
`dl.ash2txt.org` rather than the S3 open-data bucket, and asserted by
their author to be hard to find errors in. That assertion is exactly the
kind this tool can answer without a threshold to argue about, so the
result is worth having whichever way it falls: crossings found are a bug
report on a flagship dataset, and none found is the broadest positive
control the predicate has ever had.

Two deliberate departures from the corpus driver, both recorded in every
output row so no reader has to infer them:

* `pipeline.census` refuses any surface under 5,000 valid cells, because
  a verdict on a sliver of a large trace is not meaningful. Here the
  whole population is small by construction, so that floor would refuse
  the dataset rather than measure it. The engine is called directly and
  the valid-cell count is reported per patch instead.
* Sampling is deterministic given `--seed` and the sorted index, and the
  exact patch list plus per-file SHA-256 are written alongside the
  results, so a run can be repeated or contested file by file.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import http.client
import json
import os
import random
import re
import shutil
import ssl
import sys
import threading
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402

from windcheck import pipeline, tifxyz  # noqa: E402

HOST = "dl.ash2txt.org"
PREFIX = "/datasets/spiral_datasets/PHercParis4/verified_patches"
BASE = f"https://{HOST}{PREFIX}"
FILES = ("meta.json", "x.tif", "y.tif", "z.tif")
INDEX = Path("out/patches/index.txt")

# One kept-alive HTTPS connection per worker thread. Four small files per
# patch over a fresh connection each costs ~2.6 s in TLS handshakes alone;
# reusing one connection cuts that to ~1.0 s, and it is also the polite
# thing to do to a community file server. Measured, not assumed: see
# notes/DECISIONS.md for the profile that motivated it.
_local = threading.local()


def _conn() -> http.client.HTTPSConnection:
    c = getattr(_local, "conn", None)
    if c is None:
        c = http.client.HTTPSConnection(
            HOST, timeout=120, context=ssl.create_default_context())
        _local.conn = c
    return c


def _drop_conn() -> None:
    c = getattr(_local, "conn", None)
    if c is not None:
        try:
            c.close()
        except Exception:      # noqa: BLE001
            pass
        _local.conn = None


def get(path: str, retries: int = 3) -> bytes | None:
    """GET one file over this thread's connection, reconnecting on error."""
    for attempt in range(retries):
        try:
            c = _conn()
            c.request("GET", path)
            r = c.getresponse()
            body = r.read()
            if r.status == 200:
                return body
            _drop_conn()
        except Exception:      # noqa: BLE001
            _drop_conn()
        time.sleep(0.5 * (attempt + 1))
    return None


def fetch_index(force: bool = False) -> list[str]:
    """The sorted list of patch directory names, cached on disk."""
    if INDEX.exists() and not force:
        return INDEX.read_text().split()
    with urllib.request.urlopen(BASE + "/", timeout=600) as r:
        html = r.read().decode("utf-8", "replace")
    names = sorted(set(re.findall(r'href="(band-[^"]+\.tifxyz)/"', html)))
    INDEX.parent.mkdir(parents=True, exist_ok=True)
    INDEX.write_text("\n".join(names) + "\n")
    return names


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        while b := fh.read(1 << 20):
            h.update(b)
    return h.hexdigest()


def have_patch(d: Path) -> bool:
    """A patch is cached only when EVERY file is present.

    `d.is_dir()` is not this check. Under a thread pool one worker can
    create the directory while another is still writing into it, and a
    reader that trusts the directory gets a half-written patch. That race
    killed a 40,000-pair run at 3,500.
    """
    return all((d / f).exists() for f in FILES)


def download(name: str, dest: Path) -> dict[str, str] | None:
    """Fetch one patch. Returns per-file hashes, or None if incomplete.

    Each file lands under a unique temporary name and is renamed into
    place, so a concurrent reader sees either the finished file or no
    file, never a partial one.
    """
    d = dest / name
    d.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for f in FILES:
        target = d / f
        if not target.exists():
            body = get(f"{PREFIX}/{name}/{f}")
            if body is None:
                return None
            tmp = d / f".{f}.{threading.get_ident():x}.part"
            tmp.write_bytes(body)
            os.replace(tmp, target)
        hashes[f] = sha256(target)
    return hashes


def audit_one(mesh: Path, work: Path, tag: str) -> dict:
    """Both-diagonal census of one patch, bypassing the corpus size floor."""
    s = tifxyz.read(mesh)
    nv, nu = s.shape
    out = {"grid": [int(nv), int(nu)],
           "n_valid": s.n_valid,
           "n_valid_upstream_rule": s.n_valid_pipeline,
           "z_floor_cells": s.z_floor_cells}
    per = {}
    for d in (0, 1):
        _, counts = pipeline.run_engine(mesh, tag, work, d,
                                        {"threads": 1})
        per[f"d{d}"] = counts
    out["transverse_d0"] = int(per["d0"]["transverse"])
    out["transverse_d1"] = int(per["d1"]["transverse"])
    out["transverse_both"] = out["transverse_d0"] + out["transverse_d1"]
    out["coplanar_both"] = int(per["d0"].get("coplanar", 0)
                               + per["d1"].get("coplanar", 0))
    out["clean"] = out["transverse_both"] == 0
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--work", default="out/patches/work")
    ap.add_argument("--data", default="out/patches/data")
    ap.add_argument("--out", default="out/patches/audit.jsonl")
    ap.add_argument("--keep", action="store_true",
                    help="keep downloaded patches instead of deleting them")
    ap.add_argument("--refresh-index", action="store_true")
    a = ap.parse_args()

    names = fetch_index(a.refresh_index)
    print(f"index: {len(names)} patches")
    rng = random.Random(a.seed)
    chosen = names if a.sample <= 0 else rng.sample(names, min(a.sample,
                                                               len(names)))
    chosen = sorted(chosen)
    print(f"sample: {len(chosen)} (seed {a.seed})")

    data, work = Path(a.data), Path(a.work)
    out_path = Path(a.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if out_path.exists():
        done = {json.loads(l)["patch"] for l in out_path.read_text().splitlines()
                if l.strip()}
        print(f"resuming: {len(done)} already recorded")

    todo = [n for n in chosen if n not in done]
    t0 = time.time()
    n_clean = n_dirty = n_failed = 0

    def process(name: str) -> dict:
        """Download AND census one patch in the same worker.

        Both stages live in one task so they overlap across workers: the
        run is network-latency-bound (~1.0 s of round trips per patch
        against ~38 ms of census), so a serial census in the main thread
        becomes the ceiling as soon as concurrency passes about 26.
        """
        rec = {"patch": name, "source": f"{BASE}/{name}"}
        hashes = download(name, data)
        if hashes is None:
            rec["status"] = "download_failed"
            return rec
        rec["sha256"] = hashes
        try:
            rec.update(audit_one(data / name, work, name[:40]))
            rec["status"] = "ok"
        except Exception as e:                              # noqa: BLE001
            rec["status"] = "census_failed"
            rec["error"] = f"{type(e).__name__}: {e}"
        if not a.keep:
            shutil.rmtree(data / name, ignore_errors=True)
        return rec

    with out_path.open("a") as sink, \
            cf.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for rec in pool.map(process, todo):
            if rec["status"] != "ok":
                n_failed += 1
            elif rec["clean"]:
                n_clean += 1
            else:
                n_dirty += 1
            sink.write(json.dumps(rec) + "\n")
            sink.flush()
            n = n_clean + n_dirty + n_failed
            if n % 100 == 0:
                el = time.time() - t0
                rate = n / el if el else 0
                left = (len(todo) - n) / rate if rate else 0
                print(f"  {n}/{len(todo)}  clean={n_clean} "
                      f"with_crossings={n_dirty} failed={n_failed}  "
                      f"{rate:.1f}/s  eta {left/60:.0f}m", flush=True)

    print(f"\npatches censused        {n_clean + n_dirty}")
    print(f"  transverse-clean      {n_clean}")
    print(f"  with crossings        {n_dirty}")
    print(f"  failed                {n_failed}")
    print(f"wall                    {time.time() - t0:.0f}s")
    print(f"records                 {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
