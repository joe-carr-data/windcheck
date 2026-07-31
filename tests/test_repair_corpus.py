"""Driver-logic tests for bench/repair_corpus.py -- no real repairs.

Pins: enumeration follows the executor's mesh-glob discipline; sharding
is a deterministic partition; the checkpoint skips only on cert + input
hashes + code_commit all matching; the summary is idempotent across
resumes; the child subprocess is always mocked via --driver-cmd.
"""
from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, "bench")
import repair_corpus as rc                                     # noqa: E402

VOL = "20990101000000"


def make_corpus(tmp_path: Path, names=("segB", "segA", "segC")):
    """Synthetic corpus root: mesh dirs matching mesh/*<VOL>*.tifxyz,
    plus decoys the enumerator must skip."""
    root = tmp_path / "rootA"
    for name in names:
        m = root / name / "mesh" / f"{name}-on-{VOL}-9.0um.tifxyz"
        m.mkdir(parents=True)
        for ax in rc.AXES:
            (m / f"{ax}.tif").write_bytes(f"{name}-{ax}".encode())
    (root / "seg_without_mesh" / "mesh").mkdir(parents=True)
    (root / "seg_wrong_volume" / "mesh"
     / f"x-on-11111111-9.0um.tifxyz").mkdir(parents=True)
    (root / "stray_file.txt").write_text("not a segment dir")
    return [("TestCorpus", str(root), VOL, "out/unused"),
            ("GhostCorpus", str(tmp_path / "missing_root"), VOL, "out/x")]


def write_cert(out: Path, row: dict, commit: str) -> Path:
    """Fabricated certificate whose recorded input hashes match the
    current synthetic mesh files."""
    cert = {
        "segment_class": "prefilter_all_infeasible",
        "transactions": [{"txn": 0}],
        "final_events": {"0": 2, "1": 1},
        "instrumentation": {"wall_seconds": 1.5},
        "code_commit": commit,
        "hashes": {"input": {f"{ax}.tif": rc.sha(row["mesh"] / f"{ax}.tif")
                             for ax in rc.AXES}},
    }
    out.mkdir(parents=True, exist_ok=True)
    p = out / f"{row['segment']}_multi_certificate.json"
    p.write_text(json.dumps(cert))
    return p


# ------------------------------------------------------------- enumeration

def test_enumeration_sorted_and_mesh_gated(tmp_path):
    corp = make_corpus(tmp_path)
    rows = rc.enumerate_segments(corp)
    assert [r["segment"] for r in rows] == ["segA", "segB", "segC"]
    assert all(r["corpus"] == "TestCorpus" for r in rows)
    # the enumerated mesh is exactly what repair_multi will resolve
    for r in rows:
        assert r["mesh"].name.startswith(r["segment"])
        assert VOL in r["mesh"].name
    # deterministic across calls
    assert rows == rc.enumerate_segments(corp)


def test_sharding_partition(tmp_path):
    corp = make_corpus(tmp_path,
                       names=[f"seg{i:02d}" for i in range(7)])
    rows = rc.enumerate_segments(corp)
    s0 = rc.select_segments(rows, None, None, "0/2", None)
    s1 = rc.select_segments(rows, None, None, "1/2", None)
    n0 = {r["segment"] for r in s0}
    n1 = {r["segment"] for r in s1}
    assert n0 | n1 == {r["segment"] for r in rows}     # union = all
    assert n0 & n1 == set()                            # disjoint
    # deterministic: same shard again is identical, order included
    assert s0 == rc.select_segments(rows, None, None, "0/2", None)


def test_selection_filters(tmp_path):
    rows = rc.enumerate_segments(make_corpus(tmp_path))
    assert rc.select_segments(rows, "ghost", None, None, None) == []
    assert len(rc.select_segments(rows, "testcorp", None, None, 2)) == 2
    f = tmp_path / "picks.txt"
    f.write_text("segC\n# comment\n\nsegA\nnot_a_segment\n")
    picked = rc.select_segments(rows, None, f, None, None)
    assert [r["segment"] for r in picked] == ["segA", "segC"]


# -------------------------------------------------------------- checkpoint

def test_checkpoint_logic(tmp_path):
    rows = rc.enumerate_segments(make_corpus(tmp_path))
    row = rows[0]
    out = tmp_path / "out"
    # no certificate -> run
    assert rc.checkpoint(row, "HEAD1", out=out)[0] == "run"
    write_cert(out, row, "HEAD1")
    # everything matches -> skip, cert returned
    action, cert = rc.checkpoint(row, "HEAD1", out=out)
    assert action == "skip" and cert["code_commit"] == "HEAD1"
    # --force -> rerun even though everything matches
    assert rc.checkpoint(row, "HEAD1", force=True,
                         out=out)[0] == "rerun_forced"
    # commit mismatch -> rerun
    assert rc.checkpoint(row, "HEAD2",
                         out=out)[0] == "rerun_code_commit_mismatch"
    # input hash mismatch -> rerun
    (row["mesh"] / "x.tif").write_bytes(b"mutated mesh")
    assert rc.checkpoint(row, "HEAD1",
                         out=out)[0] == "rerun_input_hash_mismatch"
    # corrupt certificate -> rerun, not crash
    (out / f"{row['segment']}_multi_certificate.json").write_text("{oops")
    assert rc.checkpoint(row, "HEAD1",
                         out=out)[0] == "rerun_unreadable_cert"


def test_cert_fields_parses_json_and_int_keys():
    base = {"segment_class": "clean", "transactions": [1, 2],
            "instrumentation": {"wall_seconds": 3.0}, "code_commit": "c"}
    a = rc.cert_fields({**base, "final_events": {"0": 4, "1": 5}})
    b = rc.cert_fields({**base, "final_events": {0: 4, 1: 5}})
    assert a["residual"] == b["residual"] == {"d0": 4, "d1": 5}
    assert a["repaired"] == 2 and a["class"] == "clean"


# ----------------------------------------------------------------- summary

def test_summary_idempotent(tmp_path):
    p = tmp_path / "corpus_summary.jsonl"
    recs = rc.load_summary(p)
    assert recs == {}
    recs["segA"] = {"segment": "segA", "status": "error"}
    recs["segB"] = {"segment": "segB", "status": "ok"}
    rc.write_summary(p, recs)
    # rerun of segA replaces its record; still one line per segment
    recs = rc.load_summary(p)
    recs["segA"] = {"segment": "segA", "status": "ok", "wall_s": 9}
    rc.write_summary(p, recs)
    lines = [json.loads(x) for x in p.read_text().splitlines()]
    assert len(lines) == 2
    assert {r["segment"]: r["status"] for r in lines} == \
        {"segA": "ok", "segB": "ok"}
    # sorted, stable ordering
    assert [r["segment"] for r in lines] == ["segA", "segB"]


# ------------------------------------------------- mocked end-to-end main

def fake_driver(tmp_path: Path, out: Path) -> str:
    """A --driver-cmd stand-in: writes a valid certificate for the
    requested segment from a pre-baked table. Never runs repair_multi."""
    table = tmp_path / "fake_certs.json"
    script = tmp_path / "fake_repair.py"
    script.write_text(textwrap.dedent(f"""\
        import json, sys
        from pathlib import Path
        seg = sys.argv[sys.argv.index("--segment") + 1]
        table = json.loads(Path({str(table)!r}).read_text())
        if seg not in table:
            sys.exit(3)
        out = Path({str(out)!r})
        out.mkdir(parents=True, exist_ok=True)
        (out / (seg + "_multi_certificate.json")).write_text(
            json.dumps(table[seg]))
        print("fake repair of", seg)
        """))
    return f"{sys.executable} {script}"


def bake_table(tmp_path: Path, rows, commit: str, skip=()):
    table = {}
    for row in rows:
        if row["segment"] in skip:
            continue
        table[row["segment"]] = {
            "segment_class": "clean", "transactions": [],
            "final_events": {"0": 0, "1": 0},
            "instrumentation": {"wall_seconds": 0.1},
            "code_commit": commit,
            "hashes": {"input": {f"{ax}.tif":
                                 rc.sha(row["mesh"] / f"{ax}.tif")
                                 for ax in rc.AXES}}}
    (tmp_path / "fake_certs.json").write_text(json.dumps(table))


def test_main_run_then_checkpoint_skip(tmp_path, monkeypatch, capsys):
    corp = make_corpus(tmp_path)
    monkeypatch.setattr(rc, "CORPORA", corp)
    out = tmp_path / "out"
    rows = rc.enumerate_segments(corp)
    commit = rc.head_commit()
    bake_table(tmp_path, rows, commit)
    cmd = fake_driver(tmp_path, out)
    argv = ["--jobs", "2", "--driver-cmd", cmd, "--out-root", str(out),
            "--segment-timeout", "60"]
    assert rc.main(argv) == 0
    recs = rc.load_summary(out / "corpus_summary.jsonl")
    assert len(recs) == 3
    assert all(r["status"] == "ok" for r in recs.values())
    assert all((out / f"{s}_multi_certificate.json").exists() for s in recs)
    # second invocation: everything checkpoint-skips, records stay "ok"
    assert rc.main(argv) == 0
    txt = capsys.readouterr().out
    assert txt.count("[checkpoint] SKIP") == 3
    recs2 = rc.load_summary(out / "corpus_summary.jsonl")
    assert len(recs2) == 3
    assert all(r["status"] == "ok" for r in recs2.values())
    # --force reruns despite matching checkpoints
    assert rc.main(argv + ["--force"]) == 0
    assert "RERUN_FORCED" in capsys.readouterr().out


def test_main_child_failure_recorded(tmp_path, monkeypatch):
    corp = make_corpus(tmp_path)
    monkeypatch.setattr(rc, "CORPORA", corp)
    out = tmp_path / "out"
    rows = rc.enumerate_segments(corp)
    commit = rc.head_commit()
    bake_table(tmp_path, rows, commit, skip=("segB",))   # segB exits 3
    cmd = fake_driver(tmp_path, out)
    assert rc.main(["--jobs", "1", "--driver-cmd", cmd,
                    "--out-root", str(out)]) == 2
    recs = rc.load_summary(out / "corpus_summary.jsonl")
    assert recs["segB"]["status"] == "error"
    assert recs["segB"]["returncode"] == 3
    assert recs["segA"]["status"] == "ok"
    assert recs["segC"]["status"] == "ok"
    assert (out / "corpus_logs" / "segB.log").exists()


def test_main_segment_timeout(tmp_path, monkeypatch):
    corp = make_corpus(tmp_path, names=("segA",))
    monkeypatch.setattr(rc, "CORPORA", corp)
    out = tmp_path / "out"
    slow = tmp_path / "slow.py"
    slow.write_text("import time\ntime.sleep(30)\n")
    assert rc.main(["--jobs", "1", "--segment-timeout", "0.5",
                    "--driver-cmd", f"{sys.executable} {slow}",
                    "--out-root", str(out)]) == 2
    recs = rc.load_summary(out / "corpus_summary.jsonl")
    assert recs["segA"]["status"] == "timeout"
    assert recs["segA"]["wall_s"] < 10
