"""Census the published samples the pinned corpus does not cover.

The August target says "every published Herculaneum surface trace". The
pinned corpus is five samples and 185 traces; the open-data bucket
publishes 315 segment directories across fourteen samples. This closes
the other nine.

Each sample publishes its meshes under one or more volume ids, and a
census is only meaningful against ONE of them: the grid indices a census
emits are indices into the volume it read, and mixing two resolutions of
the same trace is the exact bug that produced a fictitious cluster in
July. So the volume is chosen per sample by majority coverage, recorded,
and every downstream record names it.

Nothing here changes the pinned corpus or any published figure. It is an
additive second corpus with its own manifest, so the two can be reported
separately or together without either being restated.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from windcheck import pipeline, tifxyz  # noqa: E402

BUCKET = "vesuvius-challenge-open-data"
KNOWN = {"PHercParis4", "PHerc0172", "PHerc0814", "PHerc0139", "PHerc1667"}
ROOT = Path("data/expand")
OUT = Path("out/expand")


def s3(*args: str) -> str:
    return subprocess.run(["aws", "s3", *args, "--no-sign-request"],
                          capture_output=True, text=True).stdout


def samples_with_segments() -> list[str]:
    out = s3("ls", f"s3://{BUCKET}/")
    names = [l.split()[-1].rstrip("/") for l in out.splitlines() if "PRE" in l]
    return [n for n in names if n not in KNOWN]


def plan_sample(sample: str) -> dict:
    """Every censusable mesh in one sample, with its S3 prefix.

    Two published layouts, both volume-space traced geometry:

      <seg>/mesh/<name>-on-<volume>-<res>.tifxyz/   five pinned samples
      <seg>/mesh/intermediate/tifxyz_original/      nine others

    `tifxyz_normalized` beside the second is flattening output rather
    than traced geometry, and is excluded. Where a volume id appears in
    the path it is parsed and the majority one chosen, because grid
    indices are indices into the volume that produced them and mixing two
    resolutions of one trace is a bug we have already paid for once.
    """
    listing = s3("ls", f"s3://{BUCKET}/{sample}/segments/", "--recursive")
    direct: dict[str, dict[str, str]] = {}
    inter: dict[str, str] = {}
    vols: Counter = Counter()
    for line in listing.splitlines():
        parts = line.split()
        if not parts:
            continue
        key = parts[-1]
        if not key.endswith("/x.tif"):
            continue
        meshdir = key.rsplit("/", 1)[0]
        seg = key.split("/segments/")[1].split("/")[0]
        if meshdir.endswith("/intermediate/tifxyz_original"):
            inter[seg] = meshdir
        elif ".tifxyz" in meshdir:
            m = re.search(r"-on-(\d+)-", meshdir)
            if m:
                vols[m.group(1)] += 1
                direct.setdefault(seg, {})[m.group(1)] = meshdir

    meshes = {}
    volume = ""
    if direct:
        volume = vols.most_common(1)[0][0]
        for seg, byvol in direct.items():
            if volume in byvol:
                meshes[seg] = byvol[volume]
    for seg, d in inter.items():
        meshes.setdefault(seg, d)
    layout = ("direct" if direct and not inter
              else "intermediate" if inter and not direct else "mixed")
    return {"volume": volume, "layout": layout, "n_meshes": len(meshes),
            "volumes_seen": dict(vols.most_common()), "meshes": meshes}


def fetch_mesh(prefix: str, dest: Path) -> bool:
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run(["aws", "s3", "cp", f"s3://{BUCKET}/{prefix}/", str(dest),
                    "--no-sign-request", "--recursive", "--only-show-errors"],
                   check=False)
    return (dest / "x.tif").exists()


def census_one(mesh: Path, work: Path, tag: str) -> dict:
    s = tifxyz.read(mesh)
    rec = {"grid": [int(x) for x in s.shape], "n_valid": s.n_valid,
           "n_valid_upstream_rule": s.n_valid_pipeline,
           "z_floor_cells": s.z_floor_cells}
    tot = 0
    for d in (0, 1):
        _, c = pipeline.run_engine(mesh, tag, work, d)
        rec[f"transverse_d{d}"] = int(c["transverse"])
        tot += int(c["transverse"])
    rec["transverse_both"] = tot
    rec["clean"] = tot == 0
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["plan", "run"])
    ap.add_argument("--out", default=str(OUT / "expand.jsonl"))
    ap.add_argument("--plan-file", default=str(OUT / "plan.json"))
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    if a.stage == "plan":
        plan = {}
        for s in samples_with_segments():
            info = plan_sample(s)
            if info["n_meshes"]:
                plan[s] = info
                print(f"{s:16s} {info['layout']:12s} {info['n_meshes']:3d} traces"
                      + (f"  volume {info['volume']}" if info["volume"] else ""))
        Path(a.plan_file).write_text(json.dumps(plan, indent=1) + "\n")
        print(f"\n{len(plan)} samples, "
              f"{sum(p['n_meshes'] for p in plan.values())} traces")
        print(f"wrote {a.plan_file}")
        return 0

    plan = json.loads(Path(a.plan_file).read_text())
    out_path = Path(a.out)
    done = set()
    if out_path.exists():
        done = {json.loads(l)["segment"] for l in out_path.read_text().splitlines()
                if l.strip()}
    work = OUT / "work"
    work.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    with out_path.open("a") as sink:
        for sample, p in plan.items():
            for seg, prefix in sorted(p["meshes"].items()):
                if seg in done:
                    continue
                dest = ROOT / sample / seg
                rec = {"sample": sample, "segment": seg,
                       "volume": p["volume"], "layout": p["layout"],
                       "s3_prefix": prefix}
                if not fetch_mesh(prefix, dest):
                    rec["status"] = "fetch_failed"
                else:
                    try:
                        rec.update(census_one(dest, work, seg[:40]))
                        rec["status"] = "ok"
                    except Exception as e:                  # noqa: BLE001
                        rec["status"] = f"{type(e).__name__}: {e}"
                sink.write(json.dumps(rec) + "\n")
                sink.flush()
                v = rec.get("transverse_both")
                print(f"  {sample:14s} {seg[:40]:40s} "
                      f"{'clean' if rec.get('clean') else (v if v is not None else rec['status'][:22])}",
                      flush=True)

    recs = [json.loads(l) for l in out_path.read_text().splitlines() if l.strip()]
    ok = [r for r in recs if r.get("status") == "ok"]
    print(f"\ncensused {len(ok)} of {len(recs)} records in "
          f"{time.time()-t0:.0f}s")
    print(f"clean {sum(1 for r in ok if r['clean'])}, "
          f"self-intersecting {sum(1 for r in ok if not r['clean'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
