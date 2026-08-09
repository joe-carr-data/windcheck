"""Give every operand a locator someone else can actually fetch.

The package references its inputs and reference derivatives by hash and
by the path they had when they were captured. A hash proves you have the
right bytes; it does not tell you where to get them. Without a locator,
"these are published already" is true and useless.

Per trace this emits:

  input      where the mesh that was CUT comes from -- an upstream S3 URI
             when it is the published original, or the release asset that
             carries it when it is a certified displacement-repaired base
  reference  the release tag, asset, member path inside the archive, and
             the archive's sha256
  fetch      a command that gets it

Coverage is reported, not assumed. A trace whose operand could not be
located says so.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

OPEN_DATA_BUCKET = "vesuvius-challenge-open-data"
RELEASE_URL = "https://github.com/{repo}/releases/download/{tag}/{asset}"


def upstream_index() -> dict:
    """tail path -> (bucket, s3_key) from the pinned data manifests."""
    out = {}
    for f in sorted(glob.glob("data/MANIFEST*.json")):
        d = json.loads(Path(f).read_text())
        b = d.get("bucket")
        for e in d.get("files", []):
            out[e["path"]] = (b, e["s3_key"])
    return out


def find_upstream(mesh: Path, idx: dict) -> dict | None:
    # data/<sample>_tifxyz/<segment>/mesh/<name>.tifxyz -> the manifest
    # records <segment>/mesh/<name>.tifxyz/<file>
    tail = "/".join(mesh.parts[2:])
    for path, (bucket, key) in idx.items():
        if path.startswith(tail + "/"):
            prefix = key[: key.rindex("/")]
            return {"bucket": bucket, "s3_prefix": prefix,
                    "uri": f"s3://{bucket}/{prefix}/"}
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", default="out/benchmark")
    ap.add_argument("--repo", default="joe-carr-data/windcheck")
    ap.add_argument("--repack-report", default="notes/release-repack-report.json")
    ap.add_argument("--correction-manifest",
                    default="notes/mask-correction-manifest.json")
    ap.add_argument("--indices", nargs="+",
                    default=["out/release/index.json",
                             "out/release/expansion_index.json"])
    a = ap.parse_args()

    entries: dict = {}
    for p in a.indices:
        entries.update(json.loads(Path(p).read_text())["segments"])
    bench = json.loads((Path(a.benchmark) / "index.json").read_text())

    up = upstream_index()

    # member path inside a release archive, per output mesh
    member: dict = {}
    archive_sha: dict = {}
    rr = Path(a.repack_report)
    if rr.is_file():
        for asset in json.loads(rr.read_text())["assets"]:
            if asset.get("sha256_after"):
                archive_sha[asset["asset"]] = asset["sha256_after"]
            for c in asset.get("changes") or []:
                surf = str(Path(c["member"]).parent)
                member[Path(surf).name] = {"tag": asset["tag"],
                                           "asset": asset["asset"],
                                           "member": surf}
    cm = Path(a.correction_manifest)
    if cm.is_file():
        for ar in json.loads(cm.read_text()).get("archives", []):
            if ar.get("sha256_after_correction"):
                archive_sha[ar["asset"]] = ar["sha256_after_correction"]

    rows, n_in, n_ref = {}, 0, 0
    for t in bench["traces"]:
        if not t.get("scored"):
            continue
        seg = t["segment"]
        e = entries[seg]
        rec: dict = {"segment": seg}

        # ---- the input that was cut ----------------------------------
        base_kind = e.get("base_kind")
        inp = Path(e["input_mesh"])
        if base_kind == "original":
            loc = (find_upstream(Path(e.get("original_mesh") or inp), up)
                   or (({"bucket": OPEN_DATA_BUCKET,
                         "s3_prefix": e["s3_prefix"],
                         "uri": f"s3://{OPEN_DATA_BUCKET}/{e['s3_prefix']}/"}
                        if e.get("s3_prefix") else None)))
            rec["input"] = {
                "role": "the published original trace",
                "upstream": loc,
                "hashes": e.get("input_hashes"),
                "fetch": (f"aws s3 sync --no-sign-request {loc['uri']} "
                          f"./{seg}.tifxyz/" if loc else None)}
        else:
            key = f"{seg}_repaired.tifxyz"
            rec["input"] = {
                "role": ("a certified displacement-repaired base, NOT the "
                         "published original: this project produced it and "
                         "publishes it in the repaired-bases release asset"),
                "release_asset": "repaired-bases",
                "member_hint": key,
                "hashes": e.get("input_hashes"),
                "original_upstream": find_upstream(
                    Path(e["original_mesh"]), up) if e.get("original_mesh")
                    else None}
        if rec["input"].get("upstream") or rec["input"].get("release_asset"):
            n_in += 1

        # ---- the reference derivative --------------------------------
        om = e.get("output_mesh")
        if om:
            m = member.get(Path(om).name)
            if m:
                url = RELEASE_URL.format(repo=a.repo, tag=m["tag"],
                                         asset=m["asset"])
                rec["reference"] = {
                    "release_tag": m["tag"], "asset": m["asset"],
                    "url": url,
                    "archive_sha256": archive_sha.get(m["asset"]),
                    "member": m["member"],
                    "hashes": e.get("output_hashes"),
                    "fetch": (f"curl -L -O {url} && tar xzf {m['asset']} "
                              f"{m['member']}")}
                n_ref += 1
            else:
                rec["reference"] = {
                    "unlocated": ("this derivative is published but the "
                                  "repack report does not name its archive "
                                  "member")}
        else:
            rec["reference"] = {
                "identical_to_input": True,
                "note": "no cut was made; the reference IS the input"}
            n_ref += 1
        rows[seg] = rec

    out = {
        "schema": "windcheck_topology_benchmark/v1#locators",
        "what": ("where to obtain every operand this package references "
                 "but does not re-host"),
        "release_repo": a.repo,
        "coverage": {"n_traces": len(rows),
                     "inputs_located": n_in,
                     "references_located": n_ref},
        "caveat": ("a locator is a pointer, not a guarantee: verify the "
                   "hashes in each trace's provenance.json after fetching"),
        "traces": rows,
    }
    p = Path(a.benchmark) / "locators.json"
    p.write_text(json.dumps(out, indent=1) + "\n")
    print(f"{n_in}/{len(rows)} inputs and {n_ref}/{len(rows)} references "
          f"located -> {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
