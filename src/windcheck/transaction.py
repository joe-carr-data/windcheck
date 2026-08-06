"""The topology transaction: stage, optionally transform, verify, commit.

One command takes a candidate tifxyz, applies the FROZEN transform
policy only if the input is not already transverse-clean, verifies the
result the hard way (retained pixels byte-identical, sidecars preserved
at retained pixels and invalidated at excised ones, reload census 0/0
under both triangulations, optionally the official volume-cartographer
validator), writes a hash-bound certificate, and atomically promotes
the result to the output path. The input is never modified, and nothing
is promoted unless every transformation and verification gate succeeds.

The authoritative transaction report is committed INSIDE the output
(windcheck_transaction/transaction_report.json), written before
promotion so it exists iff the output does. `--report` is a post-commit
copy: if that copy fails, a warning is printed and the exit code of the
committed transaction stands -- a promoted output is never rolled back
over a report-copy failure.

Exit codes (stable, for CI):
    0   committed: input was already transverse-clean; no
        retained coordinate changed (mask-only-invalid coordinates
        may be normalized to the -1 sentinel)
    10  committed, transformed and verified clean
    3   refused: transform could not produce a verified-clean result
        within the frozen policy, or a verification check failed
    2   invalid or unverifiable input (missing bands, unknown sidecars
        without an adapter that declares them, bad mask semantics)
    1   internal failure

Adapters declare sidecar semantics. `generic` accepts only the core
tifxyz files; any unknown file in the input is a refusal (exit 2) --
the transaction never guesses how to invalidate what it does not know.
`scrollfiesta` additionally handles provenance.tif, winding.tif and
*_facemap.i32 (preserved at retained pixels, zeroed/-1 at excised).

No physical-correctness claim is made anywhere: the output is a
validator-clean derivative under the stated validator.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import tifffile

REPO_ROOT = Path(__file__).resolve().parents[2]

CORE_FILES = {"x.tif", "y.tif", "z.tif", "mask.tif",
              "meta.json"}
ADAPTER_FILES = {
    "generic": set(),
    "scrollfiesta": {"provenance.tif", "winding.tif"},
}


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Refusal(Exception):
    def __init__(self, code: int, msg: str):
        self.code = code
        super().__init__(msg)


def run(argv_ns) -> int:
    src = Path(argv_ns.input)
    out = Path(argv_ns.out)
    adapter = argv_ns.adapter
    report_p = Path(argv_ns.report) if argv_ns.report else None
    rep: dict = {"schema": "windcheck_transaction/v1",
                 "generated_utc": datetime.now(timezone.utc)
                 .strftime("%Y-%m-%dT%H:%M:%SZ"),
                 "input": str(src), "adapter": adapter,
                 "policy": "frozen round-28 excision policy, unchanged",
                 "committed": False}

    def finish(code: int, note: str) -> int:
        rep["exit_code"] = code
        rep["note"] = note
        if report_p:
            # Never let the report copy change the outcome: on a refusal
            # it would mask the refusal code, and after promotion it would
            # contradict an output that already exists.
            try:
                report_p.parent.mkdir(parents=True, exist_ok=True)
                report_p.write_text(json.dumps(rep, indent=1))
            except OSError as e:
                print(f"transaction: WARNING: report copy failed: {e}",
                      file=sys.stderr)
        print(f"transaction: {note} (exit {code})")
        return code

    try:
        if not src.is_dir():
            raise Refusal(2, "input is not a tifxyz directory")
        srcr = src.resolve()
        # Validate the report path FIRST, and drop it before refusing:
        # writing the refusal report to a rejected location would be the
        # very modification the check exists to prevent.
        if report_p:
            rp = report_p.resolve()
            if rp.is_relative_to(srcr):
                report_p = None
                raise Refusal(2, "report path lies inside the input")
            if rp.is_relative_to(out.resolve()):
                report_p = None
                raise Refusal(2, "report path lies beneath the output")
        if out.exists():
            raise Refusal(2, f"output path exists: {out}")
        if out.resolve().is_relative_to(srcr):
            raise Refusal(2, "output path lies inside the input")
        ov_exe = None
        if argv_ns.official_validator:
            # Resolve ONCE (PATH lookup included) and use the resolved
            # absolute path for both invocation and hashing -- hashing the
            # unresolved argument while PATH decides what actually runs
            # would bind the certificate to the wrong binary.
            which = shutil.which(argv_ns.official_validator)
            if not which:
                raise Refusal(2, "official validator not found or not "
                                 f"executable: {argv_ns.official_validator}")
            ov_exe = str(Path(which).resolve())
            rep["official_validator_path"] = ov_exe
        subdirs = [p.name for p in src.iterdir() if p.is_dir()]
        if subdirs:
            raise Refusal(2, f"unknown subdirectories in input: {subdirs}")
        names = {p.name for p in src.iterdir() if p.is_file()}
        for req in ("x.tif", "y.tif", "z.tif", "meta.json"):
            if req not in names:
                raise Refusal(2, f"missing {req}")
        known = CORE_FILES | ADAPTER_FILES[adapter]
        unknown = sorted(n for n in names
                         if n not in known
                         and not (adapter == "scrollfiesta"
                                  and n.endswith("_facemap.i32")))
        if unknown:
            raise Refusal(2, f"unknown sidecars for adapter '{adapter}': "
                             f"{unknown} -- refusing to guess their "
                             "invalidation semantics")
        rep["input_files_sha256"] = {n: sha(src / n)
                                     for n in sorted(names)}

        out.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=f".{out.name}.",
                                         dir=out.parent) as td:
            work = Path(td)
            r = subprocess.run(
                [sys.executable,
                 str(REPO_ROOT / "bench" / "fiesta_adapter.py"),
                 "--export", str(src),
                 "--workdir", str(work / "w"),
                 "--out-dir", str(work / "staged"),
                 "--report", str(work / "adapter_report.json"),
                 "--segment-name", "transaction_input"],
                capture_output=True, text=True, cwd=REPO_ROOT)
            arep_p = work / "adapter_report.json"
            if not arep_p.exists():
                raise Refusal(1, f"transform stage produced no report; "
                                 f"stderr: {r.stderr[-400:]}")
            arep = json.loads(arep_p.read_text())
            rep["transform"] = {k: arep.get(k) for k in
                                ("terminal_disposition",
                                 "retained_fraction", "excised_px",
                                 "input_valid_px", "problems",
                                 "handed_over", "output_census",
                                 "certificate_sha256",
                                 "sidecars_preserved_and_invalidated")}
            rep["transform"]["adapter_returncode"] = r.returncode
            if not arep.get("handed_over"):
                raise Refusal(3, "transform refused handover: "
                                 f"{arep.get('problems')}")
            if r.returncode != 0:
                # The report claims success but the process died after (or
                # while) writing it. Trusting either side of a disagreement
                # is guessing; refuse as an internal failure.
                raise Refusal(1, "transform report claims handover but the "
                                 f"adapter exited {r.returncode}; refusing "
                                 "the disagreement")
            staged = work / "staged"
            # preserve the certificate inside the output (inert metadata)
            meta_dir = staged / "windcheck_transaction"
            meta_dir.mkdir()
            cert_src = arep.get("certificate")
            if cert_src and Path(cert_src).exists():
                shutil.copy2(cert_src, meta_dir / "certificate.json")
                rep["transform"]["certificate_preserved"] = \
                    "windcheck_transaction/certificate.json"
            else:
                raise Refusal(1, "transform certificate not found on disk")

            if ov_exe:
                ov = subprocess.run(
                    [ov_exe, str(staged),
                     "-o", str(work / "official.json")],
                    capture_output=True, text=True)
                try:
                    oj = json.loads((work / "official.json").read_text())
                except Exception:
                    raise Refusal(3, "official validator produced no "
                                     f"report (rc={ov.returncode})")
                rep["official_validator"] = {
                    "clean": oj.get(
                        "clean_of_transverse_self_intersection"),
                    "census": [{k: c.get(k) for k in
                                ("diagonal", "transverse")}
                               for c in oj.get("census", [])]}
                rep["official_validator"]["returncode"] = ov.returncode
                rep["official_validator"]["binary_sha256"] = sha(Path(ov_exe))
                if ov.returncode != 0:
                    raise Refusal(3, "official validator nonzero exit "
                                     f"({ov.returncode})")
                if not oj.get("clean_of_transverse_self_intersection"):
                    raise Refusal(3, "official validator reports the "
                                     "staged result NOT clean")
                shutil.copy2(work / "official.json",
                             meta_dir / "official_validator.json")

            rep["output_files_sha256"] = {
                p.name: sha(p) for p in sorted(staged.iterdir())
                if p.is_file()}
            rep["output_files_sha256"].update({
                f"windcheck_transaction/{p.name}": sha(p)
                for p in sorted(meta_dir.iterdir())})
            disp = arep.get("terminal_disposition")
            if disp == "already_clean":
                code, note = 0, "committed: input already clean"
            else:
                code, note = 10, ("committed: transformed, retained "
                                  f"{arep.get('retained_fraction'):.6f}, "
                                  "verified clean")
            # The authoritative report rides inside the output, written
            # BEFORE promotion: it exists iff the output does, and a failed
            # write here refuses the whole transaction while nothing is
            # promoted yet. It cannot list its own hash, so it is excluded
            # from output_files_sha256 by construction.
            final_rep = dict(rep, committed=True, exit_code=code, note=note)
            (meta_dir / "transaction_report.json").write_text(
                json.dumps(final_rep, indent=1))
            # atomic promotion (recheck immediately before)
            if out.exists():
                raise Refusal(2, f"output path appeared: {out}")
            os.replace(staged, out)
            rep["committed"] = True
            return finish(code, note)
    except Refusal as e:
        return finish(e.code, str(e))
    except Exception as e:  # noqa: BLE001
        return finish(1, f"internal failure: {e!r}")


def add_parser(sub):
    tr = sub.add_parser(
        "transaction",
        help="stage a tifxyz, transform only if needed, verify, and "
             "atomically commit a transverse-clean result")
    tr.add_argument("input", help="candidate tifxyz directory")
    tr.add_argument("--out", required=True,
                    help="final output path (must not exist; promoted "
                         "atomically)")
    tr.add_argument("--adapter", choices=sorted(ADAPTER_FILES),
                    default="generic",
                    help="sidecar semantics; unknown files are refused")
    tr.add_argument("--official-validator", default=None,
                    help="path to vc_tifxyz_selfcross; if given, the "
                         "staged result must also pass it")
    tr.add_argument("--report", default=None,
                    help="write the transaction report JSON here")
    tr.set_defaults(func=lambda a: run(a))
    return tr
