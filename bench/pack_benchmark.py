"""Pack the paired corpus as a surface-topology benchmark (PREREG B3).

What this produces, per trace, is everything a candidate method needs to
be measured the same way ours was, and nothing that is already public
elsewhere:

  witness.json / witness.npz   where the input self-intersects
  masks.npz                    removed / retained / changed, on the input grid
  fields.json                  retention, fragmentation, displacement --
                               copied VERBATIM from the certificate
  provenance.json              hashes and policy binding for all of it

The inputs and the reference derivatives are NOT re-hosted. The inputs are
published upstream and the derivatives are published in this project's own
release; both are referenced by hash. Re-hosting 1.7 GB of derivative that
already has a public URL would create a second artifact that can drift
from the first. Instead the packer PROVES the derivative is reconstructible
from the input and the shipped mask (`reference_reconstructible`), so the
package is self-sufficient without duplicating anything.

Every number that a reader could check is checked here rather than
asserted:

  - the certificate file is re-hashed against the release index;
  - the input and output coordinate planes are re-hashed against the index;
  - the witness contact counts must equal the certificate's own
    `census_before` transverse counts, per diagonal -- this is what binds
    a workdir to the trace it claims to belong to;
  - the reference is reconstructed from input + removed mask and required
    to be array-identical to the published derivative.

Any failure REFUSES that trace and records why. A refused trace is listed
in the index with its reason; it is never silently dropped.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from topology_eval import CENSUS, load_surface  # noqa: E402

SCHEMA = "windcheck_topology_benchmark/v1"

# The terminology block is reproduced into every emitted index and README.
# PREREG B2: this wording binds every file, README and figure.
TERMINOLOGY = {
    "reference_derivative": (
        "The paired output for each trace is a REFERENCE DERIVATIVE: a "
        "topology-certified surface produced under a published, frozen "
        "policy. It is NOT ground truth. Nothing here establishes that the "
        "removed geometry was wrong papyrus, that the retained geometry is "
        "right, or that a different cut would be worse."),
    "one_admissible_answer": (
        "The reference is ONE ADMISSIBLE ANSWER, not the target. A "
        "candidate that stays clean and unfragmented while retaining MORE "
        "area is better, and the evaluator says so."),
    "no_single_score": (
        "No single scalar score is published. A leaderboard number would "
        "invite trading fragmentation against retention. The evaluator "
        "emits the vector and refuses to collapse it."),
    "banned": (
        "'ground truth' is not used of these derivatives anywhere in this "
        "package."),
}


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def plane_hashes(d: Path) -> dict:
    out = {}
    for a in "xyz":
        p = d / f"{a}.tif"
        out[a] = sha256_file(p) if p.is_file() else None
    return out


def read_witness(workdir: Path, tag: str) -> tuple[dict, dict]:
    """Transverse contacts from the census the transform consumed.

    Returns (arrays, per-diagonal counts). Only `transverse` rows are
    kept: coplanar and grazing contacts are recorded in the certificate
    but are not the defect a candidate has to fix, and mixing them into
    the witness would misstate the task.
    """
    cols = {k: [] for k in ("diagonal", "v1", "u1", "v2", "u2",
                            "tri1", "tri2")}
    counts = {}
    for d in (0, 1):
        f = workdir / f"{tag}_before_d{d}.csv"
        n = 0
        if not f.is_file():
            counts[f"d{d}"] = None
            continue
        with open(f, newline="") as fh:
            r = csv.DictReader(fh)
            need = {"v1", "u1", "v2", "u2", "verdict", "tri1", "tri2"}
            if not need <= set(r.fieldnames or []):
                raise ValueError(f"{f}: columns {r.fieldnames} lack {need}")
            for row in r:
                if row["verdict"] != "transverse":
                    continue
                n += 1
                cols["diagonal"].append(d)
                for k in ("v1", "u1", "v2", "u2", "tri1", "tri2"):
                    cols[k].append(int(row[k]))
        counts[f"d{d}"] = n
    arrays = {
        "diagonal": np.asarray(cols["diagonal"], np.uint8),
        "v1": np.asarray(cols["v1"], np.int32),
        "u1": np.asarray(cols["u1"], np.int32),
        "v2": np.asarray(cols["v2"], np.int32),
        "u2": np.asarray(cols["u2"], np.int32),
        "tri1": np.asarray(cols["tri1"], np.uint8),
        "tri2": np.asarray(cols["tri2"], np.uint8),
    }
    return arrays, counts


def fields_from_certificate(cert: dict) -> dict:
    """PREREG B3: retention, fragmentation and displacement figures copied
    VERBATIM from the certificate. Nothing is recomputed here, so this file
    cannot disagree with the published certificate."""
    cr = cert.get("component_recovery") or {}
    return {
        "copied_verbatim_from": "the segment's excision certificate",
        "grid_shape": cert.get("grid_shape"),
        "n_valid_vertices": cert.get("n_valid"),
        "n_retained_quads": cert.get("n_retained_quads"),
        "input_area_canonical": cert.get("input_area_canonical"),
        "retention": {
            "operational_retained_fraction":
                cert.get("operational_retained_fraction"),
            "headline_retained_fraction":
                cert.get("headline_retained_fraction"),
            "area": cert.get("area"),
            "headline_area": cert.get("headline_area"),
            "denominator_note": (cert.get("area") or {}).get("denominator"),
        },
        "excision": cert.get("excision"),
        "fragmentation": {
            "core_gate": cr.get("core_gate"),
            "unthresholded": cr.get("unthresholded"),
            "definition": cr.get("definition"),
            "core_gate_pass": cert.get("core_gate_pass"),
            "core_gate_note": cert.get("core_gate_note"),
        },
        "displacement": {
            "base_kind": cert.get("base_kind"),
            "coordinate_guarantee":
                (cert.get("emission") or {}).get("coordinate_guarantee"),
            "reload_checks": cert.get("reload_checks"),
            "note": ("`base_kind: displacement_repaired` means the mesh "
                     "that was cut is itself a repair of the published "
                     "original; the operational denominator is that "
                     "repaired mesh and the headline denominator is the "
                     "original. Both are carried above."),
        },
        "census_before": cert.get("census_before"),
        "census_after": cert.get("census_after"),
        "emptiness_guard": cert.get("emptiness_guard"),
    }


def pack_trace(entry: dict, cert_path: Path, out_dir: Path,
               roots: dict) -> dict:
    """Pack one trace, or refuse it with a reason."""
    seg = entry["segment"]
    row: dict = {"segment": seg, "disposition": entry.get("disposition"),
                 "scored": False, "refused": None}

    # ---- certificate identity, against the release index ---------------
    if not cert_path.is_file():
        row["refused"] = f"certificate not found: {cert_path}"
        return row
    got = sha256_file(cert_path)
    want = entry.get("certificate_sha256")
    if want and got != want:
        row["refused"] = (f"certificate sha256 {got[:16]} does not match "
                          f"the release index {want[:16]}")
        return row
    cert = json.loads(cert_path.read_text())
    if cert.get("segment") != seg:
        row["refused"] = (f"certificate names segment "
                          f"{cert.get('segment')!r}, index says {seg!r}")
        return row

    # ---- the operands, by hash -----------------------------------------
    inp = Path(entry["input_mesh"])
    if not inp.is_dir():
        row["refused"] = f"input mesh not on disk: {inp}"
        return row
    ih = plane_hashes(inp)
    for a in "xyz":
        w = (entry.get("input_hashes") or {}).get(a)
        if w and ih[a] != w:
            row["refused"] = (f"input {a}.tif is {ih[a][:16]}, the release "
                              f"index says {w[:16]}")
            return row

    out_mesh = entry.get("output_mesh")
    ref = Path(out_mesh) if out_mesh else inp
    if out_mesh:
        if not ref.is_dir():
            row["refused"] = f"reference derivative not on disk: {ref}"
            return row
        oh = plane_hashes(ref)
        for a in "xyz":
            w = (entry.get("output_hashes") or {}).get(a)
            if w and oh[a] != w:
                row["refused"] = (f"reference {a}.tif is {oh[a][:16]}, the "
                                  f"release index says {w[:16]}")
                return row
    else:
        oh = ih

    # ---- arrays ---------------------------------------------------------
    xi, yi, zi, vi = load_surface(inp)
    xr, yr, zr, vr = load_surface(ref)
    if vi.shape != vr.shape:
        row["refused"] = (f"reference grid {list(vr.shape)} != input grid "
                          f"{list(vi.shape)}")
        return row
    if cert.get("grid_shape") and list(vi.shape) != list(cert["grid_shape"]):
        row["refused"] = (f"loaded grid {list(vi.shape)} != certificate "
                          f"grid {cert['grid_shape']}")
        return row

    removed = vi & ~vr
    retained = vi & vr
    added = ~vi & vr
    changed = np.zeros_like(retained)
    if retained.any():
        changed[retained] = ~((xi[retained] == xr[retained])
                              & (yi[retained] == yr[retained])
                              & (zi[retained] == zr[retained]))

    # ---- the reference must be reconstructible from input + mask -------
    # This is what lets the package ship masks instead of a second copy of
    # a published 1.7 GB artifact. It is checked, not assumed.
    # COMPLETE arrays, not just the cells either side calls valid. The
    # first version compared `rebuilt[vr | removed]`, which left every
    # cell invalid in BOTH operands unexamined -- 29 million of them
    # across this corpus. Arbitrary values could have sat there and the
    # derivative would still have been called reconstructible.
    recon_ok = True
    for src, dst in ((xi, xr), (yi, yr), (zi, zr)):
        if src.shape != dst.shape or src.dtype != dst.dtype:
            recon_ok = False
            break
        rebuilt = np.where(removed, src.dtype.type(-1.0), src)
        if not np.array_equal(rebuilt, dst):
            recon_ok = False
            break
    recon_ok = bool(recon_ok and not added.any() and not changed.any())
    row["reference_reconstructible"] = recon_ok
    if not recon_ok:
        # Not fatal for a MOVING method, but under the frozen excision
        # policy it means something is wrong with this pair, so refuse.
        row["refused"] = ("the reference is not reconstructible from the "
                          f"input and the removed mask (added="
                          f"{int(added.sum())}, changed="
                          f"{int(changed.sum())})")
        return row

    # ---- witness, bound to the certificate's own census ----------------
    tag = hashlib.sha256((seg + "certificate").encode()).hexdigest()[:12]
    wdir = roots[entry["_inventory"]] / f"work_{tag}"
    if not wdir.is_dir():
        row["refused"] = f"census workdir not found: {wdir}"
        return row
    try:
        warr, wcounts = read_witness(wdir, tag)
    except ValueError as e:
        row["refused"] = str(e)
        return row
    cb = cert.get("census_before") or {}
    for d in (0, 1):
        want_n = (cb.get(f"d{d}") or {}).get("transverse")
        if want_n is None:
            continue
        if wcounts.get(f"d{d}") != want_n:
            row["refused"] = (
                f"witness d{d} carries {wcounts.get(f'd{d}')} transverse "
                f"contacts, the certificate's census_before says {want_n}; "
                "the workdir does not belong to this trace")
            return row
    n_contacts = int(warr["diagonal"].size)

    # ---- write ----------------------------------------------------------
    td = out_dir / seg
    td.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(td / "masks.npz",
                        removed=removed.astype(np.uint8),
                        retained=retained.astype(np.uint8),
                        changed=changed.astype(np.uint8))

    witness_npz = None
    if n_contacts:
        np.savez_compressed(td / "witness.npz", **warr)
        witness_npz = {"file": "witness.npz",
                       "sha256": sha256_file(td / "witness.npz"),
                       "arrays": {k: [int(v.size), str(v.dtype)]
                                  for k, v in warr.items()}}
    sample = [
        {"diagonal": int(warr["diagonal"][i]),
         "quad_a": [int(warr["v1"][i]), int(warr["u1"][i])],
         "tri_a": int(warr["tri1"][i]),
         "quad_b": [int(warr["v2"][i]), int(warr["u2"][i])],
         "tri_b": int(warr["tri2"][i])}
        for i in range(min(n_contacts, 50))
    ]
    (td / "witness.json").write_text(json.dumps({
        "schema": f"{SCHEMA}#witness",
        "segment": seg,
        "what": ("the non-adjacent transverse self-intersections found in "
                 "the INPUT, from the census the transform consumed. This "
                 "is what a candidate method has to fix."),
        "census_parameters": CENSUS,
        "contact_identity": ("a contact is an unordered pair of quad "
                             "origins WITH triangle indices; quad origin "
                             "(v, u) is the top-left lattice cell of the "
                             "quad, matching the emitting engine's "
                             "convention"),
        "not_included": ("coplanar and grazing contacts are recorded in "
                         "the certificate but are not part of the witness: "
                         "the invariant is transverse-clean"),
        "n_contacts_total": n_contacts,
        "n_contacts_per_diagonal": wcounts,
        "bound_to_certificate_census_before": cb,
        "payload": witness_npz,
        "payload_note": (
            "the full contact list is a compressed array file, not inline "
            "JSON: the corpus carries millions of contacts and an inline "
            "encoding would be unreadable and an order of magnitude larger. "
            "See PREREG-TOPOLOGY-BENCHMARK ADDENDUM 1."
            if witness_npz else
            "no payload file: this trace has no transverse contacts"),
        "sample_first_50": sample,
    }, indent=1) + "\n")

    (td / "fields.json").write_text(
        json.dumps({"schema": f"{SCHEMA}#fields", "segment": seg,
                    **fields_from_certificate(cert)}, indent=1) + "\n")

    prov = {
        "schema": f"{SCHEMA}#provenance",
        "segment": seg,
        "inventory": entry["_inventory"],
        "scroll_or_sample": entry.get("scroll") or entry.get("sample"),
        "voxel_um": entry.get("voxel_um"),
        "voxel_um_note": entry.get("voxel_um_note"),
        "disposition": entry.get("disposition"),
        "base_kind": entry.get("base_kind") or cert.get("base_kind"),
        "input": {
            "role": ("the mesh that was actually cut; NOT re-hosted here, "
                     "it is published upstream and referenced by hash"),
            "path_at_capture": str(inp),
            "hashes": ih,
            "s3_prefix": entry.get("s3_prefix"),
        },
        "original": {
            "path_at_capture": entry.get("original_mesh")
                               or cert.get("original_mesh"),
            "hashes": entry.get("original_hashes")
                      or cert.get("original_hashes"),
            "input_is_original":
                (cert.get("headline_area") or {}).get("input_is_original"),
        },
        "reference_derivative": {
            "role": ("published in this project's release; NOT re-hosted "
                     "here. It is reconstructible from the input and "
                     "masks.npz, which this packer verified."),
            "path_at_capture": out_mesh,
            "hashes": oh if out_mesh else None,
            "identical_to_input": out_mesh is None,
            "reconstructible_from_input_and_mask": recon_ok,
        },
        "certificate": {
            "path_at_capture": str(cert_path),
            "sha256": got,
        },
        "policy": {
            "version": cert.get("policy_version"),
            "hash": cert.get("policy_hash"),
            "invalidation": cert.get("invalidation"),
        },
        "census_parameters": CENSUS,
        "code_provenance": cert.get("code_provenance"),
        "terminology": TERMINOLOGY,
    }
    (td / "provenance.json").write_text(json.dumps(prov, indent=1) + "\n")

    row.update(
        scored=True,
        grid=list(vi.shape),
        n_valid_input=int(vi.sum()),
        n_removed_cells=int(removed.sum()),
        n_changed_cells=int(changed.sum()),
        n_contacts=n_contacts,
        n_contacts_per_diagonal=wcounts,
        retained_fraction_certificate=cert.get("operational_retained_fraction"),
        files={f.name: {"bytes": f.stat().st_size, "sha256": sha256_file(f)}
               for f in sorted(td.iterdir()) if f.is_file()},
    )
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pinned-index", default="out/release/index.json")
    ap.add_argument("--expansion-index",
                    default="out/release/expansion_index.json")
    ap.add_argument("--pinned-root", default="out/excised/corpus")
    ap.add_argument("--expansion-root", default="out/excised/expand")
    ap.add_argument("--out", default="out/benchmark")
    ap.add_argument("--limit", type=int, default=0,
                    help="pack only the first N traces (development)")
    ap.add_argument("--only", default=None, help="pack one named segment")
    a = ap.parse_args()

    roots = {"pinned": Path(a.pinned_root),
             "expansion": Path(a.expansion_root)}
    out_dir = Path(a.out) / "traces"
    out_dir.mkdir(parents=True, exist_ok=True)

    inventories = []
    for name, path in (("pinned", a.pinned_index),
                       ("expansion", a.expansion_index)):
        p = Path(path)
        if not p.is_file():
            print(f"REFUSING: index not found: {p}")
            return 2
        inventories.append((name, p, json.loads(p.read_text()),
                            sha256_file(p)))

    rows, t0 = [], time.time()
    for name, path, idx, idx_sha in inventories:
        order = idx["order"]
        for i, seg in enumerate(order):
            if a.only and seg != a.only:
                continue
            if a.limit and len(rows) >= a.limit:
                break
            entry = dict(idx["segments"][seg])
            entry["_inventory"] = name
            disp = entry.get("disposition")
            if disp == "not_censusable":
                # PREREG B6: listed with the disposition, excluded from
                # scoring, never hidden.
                # The two inventories spell this differently: the pinned
                # index has flat `not_censusable_reason` /
                # `not_censusable_evidence`, the expansion index nests
                # both under `not_censusable`. Reading only the first
                # spelling made four exclusions look evidence-free while
                # their valid-cell counts were sitting right there.
                nested = entry.get("not_censusable") or {}
                reason = (entry.get("not_censusable_reason")
                          or nested.get("reason"))
                evidence = entry.get("not_censusable_evidence") or (
                    {k: v for k, v in nested.items() if k != "reason"}
                    or None)
                if not reason or not evidence:
                    print(f"REFUSING: {seg} is excluded but carries "
                          f"{'no reason' if not reason else 'no evidence'}; "
                          "an exclusion without evidence is a hidden trace")
                    return 2
                rows.append({
                    "segment": seg, "inventory": name,
                    "disposition": disp, "scored": False,
                    "excluded_reason": reason,
                    "excluded_evidence": evidence,
                })
                continue
            cert_path = Path(entry["certificate"])
            row = pack_trace(entry, cert_path, out_dir, roots)
            row["inventory"] = name
            rows.append(row)
            mark = ("ok" if row.get("scored")
                    else f"REFUSED: {row.get('refused')}")
            print(f"[{len(rows):3d}] {name:9s} {seg:48s} {mark}", flush=True)

    n_ok = sum(1 for r in rows if r.get("scored"))
    n_ref = sum(1 for r in rows
                if not r.get("scored") and r.get("refused"))
    n_exc = sum(1 for r in rows if r.get("excluded_reason"))
    index = {
        "schema": SCHEMA,
        "what": ("A surface-topology benchmark: paired published Herculaneum "
                 "surface traces and topology-certified reference "
                 "derivatives, with the self-intersection witness for every "
                 "input and per-cell masks for every transformation."),
        "terminology": TERMINOLOGY,
        "census_parameters": CENSUS,
        "evaluator": "bench/topology_eval.py",
        "sources": {
            name: {"index": str(path), "sha256": sha,
                   "schema": idx.get("schema"),
                   "generated_utc": idx.get("generated_utc")}
            for name, path, idx, sha in inventories
        },
        "summary": {
            "n_traces_indexed": len(rows),
            "n_packaged": n_ok,
            "n_excluded_not_censusable": n_exc,
            "n_refused": n_ref,
            "n_contacts_total": sum(r.get("n_contacts") or 0 for r in rows),
            "n_removed_cells_total": sum(r.get("n_removed_cells") or 0
                                         for r in rows),
        },
        "pack_seconds": round(time.time() - t0, 1),
        "traces": rows,
    }
    Path(a.out).mkdir(parents=True, exist_ok=True)
    (Path(a.out) / "index.json").write_text(json.dumps(index, indent=1) + "\n")

    print(f"\npackaged {n_ok}, excluded {n_exc} (not censusable), "
          f"refused {n_ref}, of {len(rows)} indexed")
    print(f"witness contacts: {index['summary']['n_contacts_total']:,}")
    if n_ref:
        print("\nREFUSALS:")
        for r in rows:
            if r.get("refused"):
                print(f"  {r['segment']}: {r['refused']}")
    return 0 if n_ref == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
