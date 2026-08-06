"""Drop-in post-export topology transaction for a Scroll Fiesta export.

External review, rank 3: the scientific transform (frozen round-28 policy) is
proven; what was missing is a transaction a pipeline could actually
consume. This adapter:

  1. Accepts the export DIRECTLY, mask and all. Validity is
     sentinel-valid AND mask >= 255 (VC3D's own keep rule). If the mask
     invalidates pixels the sentinel does not, they are baked to the
     (-1,-1,-1) sentinel in the working copy fed to the operator, and
     the count is reported.
  2. Runs the UNCHANGED frozen transform (bench/excise_shadow.py
     --certificate) on the maskless working copy via a generated
     one-entry base manifest.
  3. Reassembles a drop-in output directory next to the input:
     transformed x/y/z; mask.tif rebuilt (0 where excised or originally
     invalid, 255 retained); Fiesta's grid-aligned sidecars preserved at
     retained pixels and invalidated at excised ones -- provenance.tif
     set to 0, winding.tif set to 0, an SF_TXZ_FACEMAP sidecar (when
     present next to the export) set to -1; meta.json carries the
     export's own VC3D fields with the operator's recomputed bbox and a
     transaction note.
  4. Verifies, and refuses to hand over on any failure:
     - retained pixels BYTE-IDENTICAL to the input in every band and
       every preserved sidecar;
     - excised pixels sentinel-invalid, mask 0, provenance 0;
     - reload census of the OUTPUT DIRECTORY AS SHIPPED (mask and all):
       zero transverse under both triangulations;
     - sha256 of every emitted file recorded in the transaction report.

The operator itself is never touched: same validator, same thresholds,
same selection policy, same budgets.

    uv run python bench/fiesta_adapter.py --export <tifxyz dir> \
        --workdir out/fiesta/adapter --out-dir <output tifxyz dir> \
        --report <report.json>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import tifffile

REPO_ROOT = Path(__file__).resolve().parents[1]
for _p in (REPO_ROOT / "src", REPO_ROOT / "bench"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from corpus_bases import geometry_key, head_commit               # noqa: E402
from windcheck.manifest import mesh_manifest                     # noqa: E402


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_bands(d: Path):
    x = np.asarray(tifffile.imread(d / "x.tif"))
    y = np.asarray(tifffile.imread(d / "y.tif"))
    z = np.asarray(tifffile.imread(d / "z.tif"))
    return x, y, z


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--export", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--segment-name", default="fiesta_adapter_input")
    args = ap.parse_args(argv)

    export = Path(args.export)
    work = Path(args.workdir)
    out_dir = Path(args.out_dir)
    work.mkdir(parents=True, exist_ok=True)

    report: dict = {
        "schema": "fiesta_adapter_transaction/v1",
        "generated_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "input": str(export),
        "operator": "bench/excise_shadow.py --certificate (frozen "
                    "round-28 policy, unchanged)",
    }

    # ---- 1. validity: sentinel AND mask ---------------------------------
    x, y, z = read_bands(export)
    H, W = x.shape
    sentinel_valid = ~((x == -1) & (y == -1) & (z == -1))
    mask_p = export / "mask.tif"
    mask = (np.asarray(tifffile.imread(mask_p)) if mask_p.exists()
            else np.full((H, W), 255, np.uint8))
    valid = sentinel_valid & (mask >= 255)
    baked = int((sentinel_valid & (mask < 255)).sum())
    report["input_grid"] = [int(H), int(W)]
    report["input_valid_px"] = int(valid.sum())
    report["mask_only_invalid_px_baked_to_sentinel"] = baked

    # working copy: maskless, mask baked into the sentinel
    workmesh = work / "input_baked.tifxyz"
    if workmesh.exists():
        shutil.rmtree(workmesh)
    workmesh.mkdir(parents=True)
    for band, name in ((x, "x.tif"), (y, "y.tif"), (z, "z.tif")):
        b = band.copy()
        b[~valid] = -1
        tifffile.imwrite(workmesh / name, b)
    shutil.copy2(export / "meta.json", workmesh / "meta.json")

    # ---- 2. one-entry manifest + the frozen transform -------------------
    h = mesh_manifest(workmesh)
    planes = {a: h.get(a) for a in ("x", "y", "z")}
    planes["mask"] = h.get("mask")
    manifest = {
        "schema": "corpus_bases/v1",
        "corpus_note": "fiesta_adapter one-entry manifest",
        "generated_utc": report["generated_utc"],
        "code_commit": head_commit(),
        "n_segments": 1, "n_canonical": 1, "n_duplicate_aliases": 0,
        "n_displacement_repaired": 0, "n_original": 1,
        "entries": [{
            "segment": args.segment_name,
            "corpus": "fiesta_adapter", "volume": None, "voxel_um": None,
            "voxel_um_note": "adapter input; physical units not claimed",
            "original_mesh": str(workmesh), "original_hashes": planes,
            "published_manifest_hashes": None,
            "published_manifest_source": None,
            "published_manifest_agrees_with_disk": None,
            "base_kind": "original", "base_mesh": str(workmesh),
            "base_hashes": planes,
            "repair_certificate": None, "repair_certificate_sha256": None,
            "repair_cert_output_hashes_verified": None,
            "repair_verification_note": ("adapter working copy: input "
                                         "mask baked into the sentinel "
                                         f"({baked} px beyond it)"),
            "geometry_key": geometry_key(planes),
            "original_geometry_key": geometry_key(planes),
            "duplicate_of": None, "is_canonical": True,
        }],
    }
    mpath = work / "adapter_base.json"
    mpath.write_text(json.dumps(manifest, indent=1))

    troot = work / "transform"
    r = subprocess.run(
        [sys.executable, str(REPO_ROOT / "bench" / "excise_shadow.py"),
         "--segment", args.segment_name, "--certificate",
         "--base-manifest", str(mpath), "--out-root", str(troot)],
        capture_output=True, text=True, cwd=REPO_ROOT)
    cert_p = troot / f"{args.segment_name}_excision_certificate.json"
    if not cert_p.exists():
        raise SystemExit(f"transform produced no certificate; stderr tail: "
                         f"{r.stderr[-500:]}")
    cert = json.loads(cert_p.read_text())
    report["certificate"] = str(cert_p)
    report["certificate_sha256"] = sha(cert_p)
    report["terminal_disposition"] = cert.get("terminal_disposition")
    report["retained_fraction"] = cert.get("headline_retained_fraction")
    if cert.get("terminal_disposition") not in ("transformed",
                                                "already_clean"):
        report["handed_over"] = False
        Path(args.report).write_text(json.dumps(report, indent=1))
        raise SystemExit(f"transform disposition "
                         f"{cert.get('terminal_disposition')}; NOT handing "
                         "over an output")

    if cert["terminal_disposition"] == "already_clean":
        tx, ty, tz = x, y, z
        excised = np.zeros((H, W), bool)
    else:
        tmesh = troot / f"{args.segment_name}_excised.tifxyz"
        tx, ty, tz = read_bands(tmesh)
        t_valid = ~((tx == -1) & (ty == -1) & (tz == -1))
        # the operator may add its own mask; honor it
        tm = tmesh / "mask.tif"
        if tm.exists():
            t_valid &= np.asarray(tifffile.imread(tm)).astype(bool)
        excised = valid & ~t_valid
    report["excised_px"] = int(excised.sum())
    retained = valid & ~excised

    # ---- 3. assemble the drop-in output --------------------------------
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    for band, name in ((x, "x.tif"), (y, "y.tif"), (z, "z.tif")):
        b = band.copy()
        b[~retained] = -1
        tifffile.imwrite(out_dir / name, b)
    new_mask = np.where(retained, np.uint8(255), np.uint8(0))
    tifffile.imwrite(out_dir / "mask.tif", new_mask)

    preserved = []
    for name in ("provenance.tif", "winding.tif"):
        p = export / name
        if not p.exists():
            continue
        s = np.asarray(tifffile.imread(p)).copy()
        s[excised] = 0
        tifffile.imwrite(out_dir / name, s)
        preserved.append(name)
    fm_in = None
    fms = sorted(export.glob("*_facemap.i32"))
    if len(fms) > 1:
        report["handed_over"] = False
        report["problems"] = [f"multiple facemaps in input: "
                              f"{[f.name for f in fms]}"]
        Path(args.report).write_text(json.dumps(report, indent=1))
        raise SystemExit("multiple facemaps; refusing")
    if fms:
        fm_in = fms[0]
        if fm_in.stat().st_size != H * W * 4:
            report["handed_over"] = False
            report["problems"] = ["facemap size mismatch"]
            Path(args.report).write_text(json.dumps(report, indent=1))
            raise SystemExit("facemap size mismatch; refusing")
    if fm_in is not None:
        fm = np.fromfile(fm_in, dtype="<i4")
        if len(fm) == H * W:
            fm = fm.reshape(H, W).copy()
            fm[excised] = -1
            (out_dir / fm_in.name).write_bytes(
                fm.astype("<i4").tobytes())
            preserved.append(fm_in.name)
    report["sidecars_preserved_and_invalidated"] = preserved

    meta = json.loads((export / "meta.json").read_text())
    cert_meta = None
    if cert["terminal_disposition"] == "transformed":
        cert_meta = json.loads(
            (troot / f"{args.segment_name}_excised.tifxyz" /
             "meta.json").read_text())
        if cert_meta.get("bbox"):
            meta["bbox"] = cert_meta["bbox"]
    meta["windcheck_transaction"] = {
        "certificate_sha256": report["certificate_sha256"],
        "excised_px": report["excised_px"],
        "note": ("validator-clean derivative under the stated validator; "
                 "retained coordinates bit-identical to the input export"),
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=4))

    # ---- 4. verification: refuse to hand over on ANY failure ------------
    problems = []
    ox, oy, oz = read_bands(out_dir)
    for name, inp, outp in (("x", x, ox), ("y", y, oy), ("z", z, oz)):
        if not np.array_equal(inp[retained], outp[retained]):
            problems.append(f"{name}.tif retained px not byte-identical")
        if not ((outp[~retained] == -1).all()):
            problems.append(f"{name}.tif non-retained px not sentinel")
    om = np.asarray(tifffile.imread(out_dir / "mask.tif"))
    if not ((om == 255) == retained).all():
        problems.append("mask.tif does not equal the retained set")
    for scn, zero in (("provenance.tif", 0), ("winding.tif", 0)):
        if scn in preserved:
            pin = np.asarray(tifffile.imread(export / scn))
            pout = np.asarray(tifffile.imread(out_dir / scn))
            if not np.array_equal(pin[retained], pout[retained]):
                problems.append(f"{scn} retained px altered")
            if not (pout[excised] == zero).all():
                problems.append(f"{scn} excised px not {zero}")
    if fm_in is not None and (out_dir / fm_in.name).exists():
        fin = np.fromfile(fm_in, dtype="<i4").reshape(H, W)
        fout = np.fromfile(out_dir / fm_in.name,
                           dtype="<i4").reshape(H, W)
        if not np.array_equal(fin[retained], fout[retained]):
            problems.append("facemap retained px altered")
        if not (fout[excised] == -1).all():
            problems.append("facemap excised px not -1")

    # reload census of the OUTPUT AS SHIPPED (with its mask)
    from crossing_census import census_one
    from excise_segment import CENSUS
    with tempfile.TemporaryDirectory() as td:
        row = census_one(out_dir, "adapter_out", CENSUS["exclude"],
                         CENSUS["cell"], CENSUS["threads"],
                         CENSUS["maxedge"], Path(td))
    if row is None:
        problems.append("output below census validity floor")
    else:
        report["output_census"] = {"d0": row["d0"], "d1": row["d1"]}
        if row["d0"]["transverse"] or row["d1"]["transverse"]:
            problems.append(
                f"output NOT transverse-clean: "
                f"{row['d0']['transverse']}/{row['d1']['transverse']}")

    report["file_sha256"] = {p.name: sha(p)
                             for p in sorted(out_dir.iterdir())}
    report["problems"] = problems
    report["handed_over"] = not problems
    Path(args.report).write_text(json.dumps(report, indent=1))
    if problems:
        raise SystemExit(f"TRANSACTION REFUSED: {problems}")
    print(f"transaction complete: {out_dir}  excised {report['excised_px']}"
          f" px of {report['input_valid_px']}; retained fraction "
          f"{report['retained_fraction']}; census 0/0; "
          f"{len(report['file_sha256'])} files hashed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
