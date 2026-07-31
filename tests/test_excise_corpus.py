"""Driver-logic pins for bench/excise_corpus.py -- no real corpus, no
solver, no census engine. Everything runs against a synthetic base
manifest in tmp_path with the child subprocess supplied via --driver-cmd.

What these pin is the part of a corpus pass that decides WHAT HAPPENS TO
EACH SEGMENT, because that is where a corpus claim can quietly lose an
entry:

1. the checkpoint -- skip only when base hashes AND policy hash AND code
   commit all match; each mismatch independently names itself;
2. duplicate aliases -- certificated, never executed (a cut performed once
   must not be counted twice);
3. a timeout is a RESULT: a terminal record, not a driver crash;
4. a child that exits non-zero leaves a terminal record with the log tail,
   and the driver still finishes having recorded every segment;
5. the summary is idempotent across resumes;
6. selection filters apply in the documented order and sharding is a
   partition.
"""
from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, "bench")
sys.path.insert(0, "src")
import excise_corpus as ec                                     # noqa: E402


def hashes(seed: str) -> dict:
    return {ax: f"{seed}-{ax}-sha" for ax in ec.AXES} | {"mask": None}


def entry(segment: str, corpus="TestCorpus", canonical=True,
          duplicate_of=None, base_kind="original", seed=None) -> dict:
    seed = seed or segment
    return {
        "segment": segment, "corpus": corpus, "volume": "20990101000000",
        "voxel_um": 9.0,
        "original_mesh": f"data/{segment}/mesh/{segment}.tifxyz",
        "original_hashes": hashes(seed),
        "base_kind": base_kind,
        "base_mesh": f"out/base/{segment}.tifxyz",
        "base_hashes": hashes(seed),
        "repair_certificate": None, "repair_certificate_sha256": None,
        "geometry_key": f"geom-{seed}",
        "original_geometry_key": f"geom-{seed}",
        "duplicate_of": duplicate_of,
        "is_canonical": canonical}


def make_manifest(tmp_path: Path, entries=None) -> Path:
    entries = entries if entries is not None else [
        entry("segB"), entry("segA"),
        entry("segC", corpus="OtherCorpus", base_kind="displacement_repaired"),
    ]
    p = tmp_path / "corpus_bases.json"
    p.write_text(json.dumps({"schema": "corpus_bases/v1",
                             "n_segments": len(entries),
                             "entries": entries}))
    return p


def write_cert(out: Path, e: dict, commit: str, policy: str,
               disposition="transformed", **over) -> Path:
    """A certificate shaped like the one excise_shadow --certificate
    writes, with the fields the checkpoint reads."""
    cert = {
        "record_kind": ec.RECORD_KIND,
        "segment": e["segment"],
        "terminal_disposition": disposition,
        "status": "transverse_clean",
        "policy_hash": policy,
        "source_tree_digest": commit,
        "base_hashes": e["base_hashes"],
        "operational_retained_fraction": 0.9993,
        "headline_retained_fraction": 0.9994,
        "core_gate_pass": True,
        "claimed_clean": True,
        "wall_seconds": 12.5,
        "output_mesh": f"{out}/{e['segment']}_excised.tifxyz",
        **over}
    out.mkdir(parents=True, exist_ok=True)
    p = ec.cert_path(out, e["segment"])
    p.write_text(json.dumps(cert))
    return p


# ------------------------------------------------------------- selection

def test_manifest_entries_are_sorted_by_segment(tmp_path):
    man = ec.load_manifest(make_manifest(tmp_path))
    assert [e["segment"] for e in man["entries"]] == ["segA", "segB", "segC"]
    assert man["sha256"] == ec.sha(Path(man["path"]))


def test_sharding_is_a_partition(tmp_path):
    entries = [entry(f"seg{i:02d}") for i in range(7)]
    rows = ec.load_manifest(make_manifest(tmp_path, entries))["entries"]
    s0 = ec.select_segments(rows, None, None, "0/2", None)
    s1 = ec.select_segments(rows, None, None, "1/2", None)
    n0 = {r["segment"] for r in s0}
    n1 = {r["segment"] for r in s1}
    assert n0 | n1 == {r["segment"] for r in rows}      # covers everything
    assert n0 & n1 == set()                             # disjoint
    assert s0 == ec.select_segments(rows, None, None, "0/2", None)


def test_selection_filter_order(tmp_path):
    rows = ec.load_manifest(make_manifest(tmp_path))["entries"]
    assert ec.select_segments(rows, "ghost", None, None, None) == []
    assert [r["segment"] for r in
            ec.select_segments(rows, "othercorp", None, None, None)] == \
        ["segC"]
    f = tmp_path / "picks.txt"
    f.write_text("segC\n# a comment\n\nsegA\nnot_a_segment\n")
    assert [r["segment"] for r in
            ec.select_segments(rows, None, f, None, None)] == ["segA", "segC"]
    # order is corpus -> names -> shard -> limit: the name list runs BEFORE
    # the limit, so limit 1 over {segA, segC} keeps segA, not segB
    assert [r["segment"] for r in
            ec.select_segments(rows, None, f, None, 1)] == ["segA"]
    # ... and the shard runs before the limit too
    got = ec.select_segments(rows, None, None, "1/2", 1)
    assert [r["segment"] for r in got] == ["segB"]


# ------------------------------------------------------------ checkpoint

def test_checkpoint_skip_and_each_mismatch(tmp_path):
    e = entry("segA")
    out = tmp_path / "out"
    # no certificate at all -> run
    assert ec.checkpoint(e, "HEAD1", "POL1", out=out)[0] == "run"
    write_cert(out, e, "HEAD1", "POL1")
    # everything matches -> skip, and the certificate comes back
    action, cert = ec.checkpoint(e, "HEAD1", "POL1", out=out)
    assert action == "skip" and cert["source_tree_digest"] == "HEAD1"
    # --force always runs
    assert ec.checkpoint(e, "HEAD1", "POL1", force=True,
                         out=out)[0] == "rerun_forced"
    # the three mismatches, each independently
    assert ec.checkpoint(e, "HEAD2", "POL1",
                         out=out)[0] == "rerun_code_changed"
    assert ec.checkpoint(e, "HEAD1", "POL2",
                         out=out)[0] == "rerun_policy_changed"
    moved = {**e, "base_hashes": hashes("rebuilt-base")}
    assert ec.checkpoint(moved, "HEAD1", "POL1",
                         out=out)[0] == "rerun_base_changed"


def test_checkpoint_rejects_partial_and_unreadable_certificates(tmp_path):
    e = entry("segA")
    out = tmp_path / "out"
    # a certificate that records no base hashes is a MISMATCH, not a pass
    write_cert(out, e, "HEAD1", "POL1", base_hashes={})
    assert ec.checkpoint(e, "HEAD1", "POL1",
                         out=out)[0] == "rerun_base_changed"
    write_cert(out, e, "HEAD1", "POL1",
               base_hashes={"x": e["base_hashes"]["x"]})
    assert ec.checkpoint(e, "HEAD1", "POL1",
                         out=out)[0] == "rerun_base_changed"
    ec.cert_path(out, "segA").write_text("{not json")
    assert ec.checkpoint(e, "HEAD1", "POL1",
                         out=out)[0] == "rerun_unreadable_cert"


def test_checkpoint_retries_driver_written_failures(tmp_path):
    """A timeout/crash record is a record, not a result about the geometry:
    it must not freeze the segment out of the next pass."""
    e = entry("segA")
    out = tmp_path / "out"
    write_cert(out, e, "HEAD1", "POL1", disposition="error", status="timeout")
    assert ec.checkpoint(e, "HEAD1", "POL1",
                         out=out)[0] == "rerun_previous_failure"
    for disp in ("transformed", "already_clean", "triangle_empty_invalid",
                 "not_censusable", "residual_transverse"):
        write_cert(out, e, "HEAD1", "POL1", disposition=disp)
        assert ec.checkpoint(e, "HEAD1", "POL1", out=out)[0] == "skip"


# --------------------------------------------------------------- summary

def test_summary_write_is_one_line_per_segment(tmp_path):
    p = tmp_path / "corpus_summary.jsonl"
    assert ec.load_summary(p) == {}
    recs = {"segA": {"segment": "segA", "status": "error"},
            "segB": {"segment": "segB", "status": "ok"}}
    ec.write_summary(p, recs)
    recs = ec.load_summary(p)
    recs["segA"] = {"segment": "segA", "status": "ok"}
    ec.write_summary(p, recs)
    lines = [json.loads(x) for x in p.read_text().splitlines()]
    assert [r["segment"] for r in lines] == ["segA", "segB"]
    assert all(r["status"] == "ok" for r in lines)


# ------------------------------------------------- mocked end-to-end main

def fake_driver(tmp_path: Path, mode="ok") -> str:
    """A --driver-cmd stand-in. It records every invocation (so a test can
    assert a segment was NEVER executed) and writes a certificate built
    from a table baked by `bake_table`. It never imports windcheck."""
    table = tmp_path / "fake_table.json"
    calls = tmp_path / "invocations.txt"
    script = tmp_path / f"fake_excise_{mode}.py"
    script.write_text(textwrap.dedent(f"""\
        import json, sys, time
        from pathlib import Path
        argv = sys.argv
        seg = argv[argv.index("--segment") + 1]
        out = Path(argv[argv.index("--out-root") + 1])
        man = Path(argv[argv.index("--base-manifest") + 1])
        assert "--certificate" in argv, argv
        assert man.exists(), man
        with open({str(calls)!r}, "a") as fh:
            fh.write(seg + "\\n")
        table = json.loads(Path({str(table)!r}).read_text())
        print("fake excision of", seg)
        if seg not in table:
            print("no table entry: exiting 3")
            sys.exit(3)
        cert = table[seg]
        if cert.get("__sleep"):
            time.sleep(cert["__sleep"])
        out.mkdir(parents=True, exist_ok=True)
        (out / (seg + {ec.CERT_SUFFIX!r})).write_text(json.dumps(cert))
        """))
    return f"{sys.executable} {script}"


def bake_table(tmp_path: Path, entries, commit, policy, skip=(), sleep=None):
    table = {}
    for e in entries:
        if e["segment"] in skip:
            continue
        table[e["segment"]] = {
            "record_kind": ec.RECORD_KIND,
            "segment": e["segment"],
            "terminal_disposition": "transformed",
            "status": "transverse_clean",
            "policy_hash": policy, "source_tree_digest": commit,
            "base_hashes": e["base_hashes"],
            "operational_retained_fraction": 0.999,
            "headline_retained_fraction": 0.9991,
            "core_gate_pass": True, "claimed_clean": True,
            "wall_seconds": 1.0,
            "output_mesh": str(tmp_path / f"{e['segment']}_excised.tifxyz"),
            **({"__sleep": sleep} if sleep else {})}
    (tmp_path / "fake_table.json").write_text(json.dumps(table))
    return table


def invocations(tmp_path: Path) -> list[str]:
    p = tmp_path / "invocations.txt"
    return p.read_text().split() if p.exists() else []


def base_argv(man: Path, out: Path, cmd: str, *rest) -> list[str]:
    return ["--base-manifest", str(man), "--out-root", str(out),
            "--driver-cmd", cmd, "--jobs", "2", *rest]


def test_main_runs_then_checkpoint_skips(tmp_path, capsys):
    man = make_manifest(tmp_path)
    entries = ec.load_manifest(man)["entries"]
    out = tmp_path / "out"
    commit, policy = ec.code_identity(), ec.frozen_policy_hash()
    bake_table(tmp_path, entries, commit, policy)
    argv = base_argv(man, out, fake_driver(tmp_path))
    assert ec.main(argv) == 0
    recs = ec.load_summary(out / "corpus_summary.jsonl")
    assert len(recs) == 3
    assert all(r["terminal_disposition"] == "transformed"
               for r in recs.values())
    assert all(r["certificate_sha256"] for r in recs.values())
    assert sorted(invocations(tmp_path)) == ["segA", "segB", "segC"]
    # second pass: everything checkpoint-skips, nothing is re-executed
    assert ec.main(argv) == 0
    assert "[checkpoint] SKIP" in capsys.readouterr().out
    assert sorted(invocations(tmp_path)) == ["segA", "segB", "segC"]
    recs2 = ec.load_summary(out / "corpus_summary.jsonl")
    assert len(recs2) == 3
    assert all(r["terminal_disposition"] == "transformed"
               for r in recs2.values())
    # --force reruns despite matching checkpoints
    assert ec.main(argv + ["--force"]) == 0
    assert "RERUN_FORCED" in capsys.readouterr().out
    assert len(invocations(tmp_path)) == 6


def test_summary_idempotent_across_resumes(tmp_path):
    man = make_manifest(tmp_path)
    entries = ec.load_manifest(man)["entries"]
    out = tmp_path / "out"
    bake_table(tmp_path, entries, ec.code_identity(), ec.frozen_policy_hash())
    argv = base_argv(man, out, fake_driver(tmp_path))
    ec.main(argv)
    ec.main(argv)
    lines = (out / "corpus_summary.jsonl").read_text().splitlines()
    assert len(lines) == 3
    assert [json.loads(x)["segment"] for x in lines] == \
        ["segA", "segB", "segC"]


def test_duplicate_alias_is_certificated_but_never_executed(tmp_path):
    entries = [entry("segA"),
               entry("segDup", canonical=False, duplicate_of="segA",
                     seed="segA")]
    man = make_manifest(tmp_path, entries)
    out = tmp_path / "out"
    bake_table(tmp_path, ec.load_manifest(man)["entries"],
               ec.code_identity(), ec.frozen_policy_hash())
    assert ec.main(base_argv(man, out, fake_driver(tmp_path))) == 0
    # the alias was never handed to the child
    assert invocations(tmp_path) == ["segA"]
    cert = json.loads(ec.cert_path(out, "segDup").read_text())
    assert cert["record_kind"] == ec.RECORD_KIND
    assert cert["terminal_disposition"] == "duplicate_alias"
    assert cert["canonical_segment"] == "segA"
    assert cert["is_canonical"] is False
    assert cert["output_mesh"] is None
    recs = ec.load_summary(out / "corpus_summary.jsonl")
    assert set(recs) == {"segA", "segDup"}
    assert recs["segDup"]["terminal_disposition"] == "duplicate_alias"
    assert recs["segDup"]["duplicate_of"] == "segA"
    # ... and it stays an alias across a resume
    assert ec.main(base_argv(man, out, fake_driver(tmp_path))) == 0
    assert invocations(tmp_path) == ["segA"]


def test_timeout_is_a_terminal_record_not_a_crash(tmp_path):
    entries = [entry("segA"), entry("segB")]
    man = make_manifest(tmp_path, entries)
    out = tmp_path / "out"
    bake_table(tmp_path, entries, ec.code_identity(), ec.frozen_policy_hash(),
               sleep=30)
    rc = ec.main(base_argv(man, out, fake_driver(tmp_path, "slow"),
                           "--segment-timeout", "0.5"))
    assert rc == 2                       # a recorded failure, not an exception
    recs = ec.load_summary(out / "corpus_summary.jsonl")
    assert set(recs) == {"segA", "segB"}
    for seg in ("segA", "segB"):
        assert recs[seg]["status"] == "timeout"
        assert recs[seg]["terminal_disposition"] == "error"
        assert recs[seg]["wall_seconds"] < 20
        cert = json.loads(ec.cert_path(out, seg).read_text())
        assert cert["terminal_disposition"] == "error"
        assert cert["status"] == "timeout"
        assert cert["driver_terminal_record"] is True
        assert cert["timeout_s"] == 0.5
        assert cert["elapsed_seconds"] > 0
        assert (out / "logs" / f"{seg}.log").exists()


def test_child_failure_records_log_tail_and_every_segment(tmp_path):
    entries = [entry("segA"), entry("segB"), entry("segC")]
    man = make_manifest(tmp_path, entries)
    out = tmp_path / "out"
    bake_table(tmp_path, entries, ec.code_identity(), ec.frozen_policy_hash(),
               skip=("segB",))          # segB's child exits 3
    rc = ec.main(base_argv(man, out, fake_driver(tmp_path)))
    assert rc == 2
    recs = ec.load_summary(out / "corpus_summary.jsonl")
    assert set(recs) == {"segA", "segB", "segC"}       # every segment recorded
    assert recs["segA"]["terminal_disposition"] == "transformed"
    assert recs["segC"]["terminal_disposition"] == "transformed"
    bad = recs["segB"]
    assert bad["terminal_disposition"] == "error"
    assert bad["returncode"] == 3
    assert "no table entry" in bad["log_tail"]
    cert = json.loads(ec.cert_path(out, "segB").read_text())
    assert cert["terminal_disposition"] == "error"
    assert cert["driver_terminal_record"] is True
    assert "no table entry" in cert["log_tail"]
    assert cert["base_hashes"] == entries[1]["base_hashes"]
    assert (out / "logs" / "segB.log").exists()


def test_child_that_writes_no_certificate_is_an_error(tmp_path):
    """Exit 0 without a certificate is not a pass."""
    entries = [entry("segA")]
    man = make_manifest(tmp_path, entries)
    out = tmp_path / "out"
    quiet = tmp_path / "quiet.py"
    quiet.write_text("print('did nothing at all')\n")
    rc = ec.main(base_argv(man, out, f"{sys.executable} {quiet}"))
    assert rc == 2
    rec = ec.load_summary(out / "corpus_summary.jsonl")["segA"]
    assert rec["terminal_disposition"] == "error"
    assert "wrote no certificate" in rec["reason"]
    assert json.loads(ec.cert_path(out, "segA").read_text())["status"] \
        == "error"


def test_terminal_certificate_carries_manifest_provenance(tmp_path):
    e = entry("segA", base_kind="displacement_repaired")
    man_path = make_manifest(tmp_path, [e])
    manifest = ec.load_manifest(man_path)
    out = tmp_path / "out"
    p = ec.write_terminal_certificate(
        out, e, manifest, "HEAD1", "POL1", disposition="error",
        status="timeout", wall_seconds=601.2)
    cert = json.loads(p.read_text())
    assert cert["base_manifest_sha256"] == manifest["sha256"]
    assert cert["base_mesh"] == e["base_mesh"]
    assert cert["base_kind"] == "displacement_repaired"
    assert cert["geometry_key"] == e["geometry_key"]
    assert cert["policy_version"] == ec.FROZEN_POLICY_VERSION
    assert cert["policy_hash"] == "POL1"
    assert cert["source_tree_digest"] == "HEAD1"
    assert cert["wall_seconds"] == 601.2
    assert cert["core_gate_pass"] is None       # nothing was measured


def test_no_segments_selected_is_a_clean_refusal(tmp_path):
    man = make_manifest(tmp_path)
    out = tmp_path / "out"
    assert ec.main(base_argv(man, out, "false", "--corpus", "nope")) == 1
