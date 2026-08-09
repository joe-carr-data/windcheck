"""One row per segment, and every population claim derived from those rows.

This exists because of a failure class this project has now shipped three
times: a number that is TRUE OF A SUBSET, narrated as true of the whole.

  "215 traces carry a defect"     -- 215 were PROCESSED; 190 carry one.
  "one trace fails the core gate" -- true of the 99 expansion traces;
                                     six of the 274 fail it.
  "every E-condition met"         -- true of the arithmetic; several of
                                     the gates had never run.

Each was a sentence typed from memory of a narrower record. None was a
computational error, so no test caught any of them. The fix is not more
care: it is to stop retyping aggregates. Every population number the
submission and the benchmark card make is derived HERE, from a table with
exactly one row per segment, and any two sources that disagree about a
segment stop the build instead of one of them silently winning.

What a row must survive:

  it exists exactly once, in exactly one inventory;
  every source that mentions it agrees about its inventory, disposition
  and base kind;
  its censusability agrees with whether the benchmark scored it;
  the evaluator and the certificate agree about its core gate;
  "officially clean" means transverse 0 on BOTH diagonals over a
  POSITIVE triangle count, never a clean verdict over an empty surface;
  a trace with a detected event was either transformed, or had its
  defects removed by displacement repair before the excision stage --
  and nothing else;
  a voxel scale is present or explicitly absent, never assumed.

Building it immediately surfaced a distinction nobody had written down.
`disposition` describes the INPUT BASE, the mesh actually cut, while the
pinned intrinsic spectrum was measured on the PUBLISHED ORIGINAL. For a
displacement-repaired base those are different surfaces, so six pinned
traces carry detected events on the original while their base needed no
excision at all. Both statements are true; "defect-bearing" simply has
to name its operand, and now it does.

The `claims` object at the end is the only place population figures come
from. If a claim is not derivable from the rows, it is not a claim this
project makes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

DISPOSITIONS = ("transformed", "already_clean", "not_censusable")
BASE_KINDS = ("original", "displacement_repaired")

SOURCES = {
    "release_pinned": "out/release/index.json",
    "release_expansion": "out/release/expansion_index.json",
    "benchmark_index": "out/benchmark/index.json",
    "acceptance": "out/benchmark/B7-acceptance.json",
    "intrinsic_pinned_d0": "out/spectrum_final_d0.json",
    "intrinsic_pinned_d1": "out/spectrum_final_d1.json",
    "intrinsic_expansion": "out/spectrum_expansion.json",
}


class Contradiction(SystemExit):
    pass


def _count(v, what: str) -> int:
    """A non-negative integer, and not a bool. `274.0 == 274` in Python,
    so an equality check alone accepts a float denominator."""
    if isinstance(v, bool) or not isinstance(v, int) or v < 0:
        refuse(f"{what} is {v!r}, which is not a non-negative integer")
    return v


def _fraction(v, what: str) -> float:
    if isinstance(v, bool) or not isinstance(v, (int, float)) \
            or not math.isfinite(v) or not 0.0 <= v <= 1.0:
        refuse(f"{what} is {v!r}, which is not a finite fraction in [0, 1]")
    return float(v)


def _strict_bool(v, what: str) -> bool:
    if not isinstance(v, bool):
        refuse(f"{what} is {v!r}, which is not a boolean")
    return v


def refuse(msg: str):
    raise Contradiction(f"REFUSING: {msg}")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def agree(rows: dict, seg: str, field: str, value, source: str):
    """Record a field, or refuse if a previous source said otherwise.

    Silent last-writer-wins is how two records drift apart without
    anything noticing until a human retypes one of them into prose.
    """
    row = rows[seg]
    if field in row and row[field] != value:
        refuse(f"{seg}: {source} says {field}={value!r}, "
               f"{row[field + '__from']} says {row[field]!r}")
    row[field] = value
    row[field + "__from"] = source


def build(paths: dict[str, Path]) -> dict:
    raw = {k: json.loads(p.read_text()) for k, p in paths.items()}
    rows: dict[str, dict] = {}

    # ---- 1. membership: exactly one inventory, exactly once ------------
    for key, inventory in (("release_pinned", "pinned"),
                           ("release_expansion", "expansion")):
        segs = raw[key]["segments"]
        for seg, ent in segs.items():
            if seg in rows:
                refuse(f"{seg} appears in more than one release index")
            rows[seg] = {"segment": seg}
            agree(rows, seg, "inventory", inventory, key)
            d = ent.get("disposition")
            if d not in DISPOSITIONS:
                refuse(f"{seg}: {key} disposition {d!r} is not one of "
                       f"{list(DISPOSITIONS)}")
            agree(rows, seg, "disposition", d, key)
            bk = ent.get("base_kind")
            if bk not in BASE_KINDS:
                refuse(f"{seg}: {key} base_kind {bk!r} is not one of "
                       f"{list(BASE_KINDS)}")
            agree(rows, seg, "base_kind", bk, key)
            v = ent.get("voxel_um")
            if v is not None:
                if isinstance(v, bool) or not isinstance(v, (int, float)) \
                        or not math.isfinite(v) or v <= 0:
                    refuse(f"{seg}: voxel_um {v!r} is neither a positive "
                           "finite number nor null")
            agree(rows, seg, "voxel_um", (float(v) if v is not None else None),
                  key)

    # ---- 2. the benchmark index, which must agree about all of it ------
    seen = set()
    for t in raw["benchmark_index"]["traces"]:
        seg = t["segment"]
        if seg in seen:
            refuse(f"{seg} appears twice in the benchmark index")
        seen.add(seg)
        if seg not in rows:
            refuse(f"{seg} is in the benchmark index but in no release index")
        agree(rows, seg, "inventory", t.get("inventory"), "benchmark_index")
        agree(rows, seg, "disposition", t.get("disposition"),
              "benchmark_index")
        scored = _strict_bool(t.get("scored"), f"{seg} scored")
        censusable = t.get("disposition") != "not_censusable"
        if scored != censusable:
            refuse(f"{seg}: benchmark scored={scored} but disposition says "
                   f"censusable={censusable}")
        rows[seg]["censusable"] = censusable
        rows[seg]["n_contacts"] = t.get("n_contacts")
        rows[seg]["retained_fraction"] = t.get("retained_fraction_certificate")
        if censusable:
            _count(rows[seg]["n_contacts"], f"{seg} n_contacts")
            _fraction(rows[seg]["retained_fraction"],
                      f"{seg} retained_fraction")
        else:
            for f in ("n_contacts", "retained_fraction"):
                if rows[seg][f] not in (None, 0):
                    refuse(f"{seg} is not censusable but records "
                           f"{f}={rows[seg][f]!r}")
    missing = sorted(set(rows) - seen)
    if missing:
        refuse(f"{len(missing)} released segments are absent from the "
               f"benchmark index, e.g. {missing[:3]}")

    # ---- 3. the official census and the core gate ----------------------
    acc = set()
    for r in raw["acceptance"]["results"]:
        seg = r["segment"]
        if seg in acc:
            refuse(f"{seg} appears twice in the acceptance record")
        acc.add(seg)
        if seg not in rows:
            refuse(f"{seg} is in the acceptance record but in no index")
        if not rows[seg]["censusable"]:
            refuse(f"{seg} is not censusable but carries an acceptance row")
        ev = (r.get("evaluator") or {}).get("core_gate_pass")
        ce = (r.get("certificate") or {}).get("core_gate_pass")
        for v, who in ((ev, "evaluator"), (ce, "certificate")):
            _strict_bool(v, f"{seg}: {who} core_gate_pass")
        if ev != ce:
            refuse(f"{seg}: core-gate verdict disagrees -- evaluator {ev}, "
                   f"certificate {ce}")
        rows[seg]["core_gate_pass"] = ev

        tr, tri = r.get("census_transverse") or {}, r.get("census_triangles") or {}
        for k in ("d0", "d1"):
            for src, what in ((tr, "transverse"), (tri, "triangles")):
                v = src.get(k)
                if isinstance(v, bool) or not isinstance(v, int) or v < 0:
                    refuse(f"{seg}: census {what} {k} is {v!r}")
        # A clean verdict over zero triangles is the shape this project's
        # own published mask defect took. It is not cleanliness.
        positive = all(tri[k] > 0 for k in ("d0", "d1"))
        rows[seg]["official_triangles_positive"] = positive
        rows[seg]["official_clean"] = positive and all(
            tr[k] == 0 for k in ("d0", "d1"))
    short = sorted(s for s, r in rows.items() if r["censusable"] and s not in acc)
    if short:
        refuse(f"{len(short)} censusable segments have no acceptance row, "
               f"e.g. {short[:3]}")

    # ---- 4. the intrinsic-scale analysis -------------------------------
    # Processed and defect-bearing are DIFFERENT populations. Conflating
    # them is the error this whole file exists to prevent, so they are
    # recorded as two independent fields and never derived from each
    # other.
    for r in rows.values():
        r["intrinsic_processed"] = False
        r["intrinsic_events"] = None
        r["intrinsic_events_d0"] = None
        r["intrinsic_events_d1"] = None
        r["intrinsic_operand"] = None

    # Both diagonals, for every pinned segment, exactly once. Without
    # this a missing d1 file leaves every pinned trace still marked
    # processed from d0 alone, with event totals that are quietly half a
    # measurement.
    pinned_inventory = {s for s, r in rows.items() if r["inventory"] == "pinned"}
    per_diagonal: dict[str, dict[str, int]] = {}
    for key, diag in (("intrinsic_pinned_d0", "d0"),
                      ("intrinsic_pinned_d1", "d1")):
        seen_d = set()
        for e in raw[key]:
            seg = e["segment"]
            if seg in seen_d:
                refuse(f"{seg} appears twice in {key}")
            seen_d.add(seg)
            if seg not in rows:
                refuse(f"{seg} is in {key} but in no index")
            if rows[seg]["inventory"] != "pinned":
                refuse(f"{seg} is in {key} but its inventory is "
                       f"{rows[seg]['inventory']!r}")
            n = _count(e.get("events_total"), f"{seg}: {key} events_total")
            m = _count(e.get("events_measured"), f"{seg}: {key} events_measured")
            # The submission claims every detected event was measured.
            # That claim is here or it is nowhere.
            if m != n:
                refuse(f"{seg}: {key} measured {m} of {n} events")
            rows[seg]["intrinsic_processed"] = True
            rows[seg]["intrinsic_operand"] = "published original trace"
            per_diagonal.setdefault(seg, {})[diag] = n
        if seen_d != pinned_inventory:
            missing_d = sorted(pinned_inventory - seen_d)
            extra_d = sorted(seen_d - pinned_inventory)
            refuse(f"{key} does not cover the pinned inventory exactly: "
                   f"missing {len(missing_d)} (e.g. {missing_d[:3]}), "
                   f"unexpected {len(extra_d)} (e.g. {extra_d[:3]})")

    exp = raw["intrinsic_expansion"]
    seen_expansion: set[str] = set()
    for t in exp["traces"]:
        seg = t["segment"]
        # One row per segment means ONE row. A duplicated expansion entry
        # left the processed and defect-bearing counts unchanged and so
        # was invisible -- which is precisely the kind of silence a
        # reconciliation exists to break.
        if seg in seen_expansion:
            refuse(f"{seg} appears twice in intrinsic_expansion")
        seen_expansion.add(seg)
        if seg not in rows:
            refuse(f"{seg} is in intrinsic_expansion but in no index")
        if rows[seg]["inventory"] != "expansion":
            refuse(f"{seg} is in intrinsic_expansion but its inventory is "
                   f"{rows[seg]['inventory']!r}")
        if t.get("ok") is not True:
            refuse(f"{seg}: intrinsic_expansion records ok={t.get('ok')!r}")
        diags = t.get("diagonals")
        if not isinstance(diags, dict) or set(diags) != {"d0", "d1"}:
            refuse(f"{seg}: intrinsic_expansion records diagonals "
                   f"{sorted(diags) if isinstance(diags, dict) else diags!r}, "
                   "require exactly d0 and d1")
        n = 0
        for k, dd in sorted(diags.items()):
            v = _count((dd or {}).get("events_total"),
                       f"{seg}: intrinsic_expansion {k} events_total")
            m = _count((dd or {}).get("events_measured"),
                       f"{seg}: intrinsic_expansion {k} events_measured")
            if m != v:
                refuse(f"{seg}: intrinsic_expansion {k} measured {m} of {v} "
                       "events")
            per_diagonal.setdefault(seg, {})[k] = v
            n += v
        # Every expansion base is an original, so for this half the
        # operand and the published trace are the same surface. Asserted
        # rather than assumed, because if that ever stops being true the
        # 190 count silently changes meaning.
        if rows[seg]["base_kind"] != "original":
            refuse(f"{seg}: expansion intrinsic measurement assumes the "
                   "operand is the published original, but base_kind is "
                   f"{rows[seg]['base_kind']!r}")
        rows[seg]["intrinsic_processed"] = True
        rows[seg]["intrinsic_operand"] = "published original trace"
        rows[seg]["intrinsic_events"] = n

    # Exactly the transformed expansion traces, no more and no less.
    expected_expansion = {s for s, r in rows.items()
                          if r["inventory"] == "expansion"
                          and r["disposition"] == "transformed"}
    if seen_expansion != expected_expansion:
        missing_e = sorted(expected_expansion - seen_expansion)
        extra_e = sorted(seen_expansion - expected_expansion)
        refuse("intrinsic_expansion does not cover the transformed expansion "
               f"traces exactly: missing {len(missing_e)} "
               f"(e.g. {missing_e[:3]}), unexpected {len(extra_e)} "
               f"(e.g. {extra_e[:3]})")

    for seg, per in per_diagonal.items():
        # d0 and d1 kept separately as well as summed. The sum is a sum of
        # two independent clusterings, not a count of distinct physical
        # defects, and storing only the total invites it being read as one.
        rows[seg]["intrinsic_events_d0"] = per.get("d0")
        rows[seg]["intrinsic_events_d1"] = per.get("d1")
        if rows[seg]["intrinsic_events"] is None:
            rows[seg]["intrinsic_events"] = sum(per.values())

    # Two populations that look like one. `disposition` describes the
    # INPUT BASE -- the mesh actually cut -- while the pinned intrinsic
    # spectrum was measured on the PUBLISHED ORIGINAL. For a trace whose
    # base is displacement-repaired those are different surfaces, and a
    # trace can carry detected events on the original while its base
    # needed no excision at all. Six pinned traces are exactly that, and
    # they are why "defect-bearing" must always name its operand.
    for seg, r in rows.items():
        ev, dis, bk = r["intrinsic_events"], r["disposition"], r["base_kind"]
        if ev and dis == "transformed":
            pass                                  # the ordinary case
        elif ev and dis == "already_clean" and bk == "displacement_repaired":
            r["defects_removed_by_repair_not_excision"] = True
        elif ev:
            refuse(f"{seg} carries {ev} detected events with disposition "
                   f"{dis!r} and base_kind {bk!r}: a trace with defects and "
                   "no excision must be one whose base was repaired")
        if r["intrinsic_processed"] and not r["censusable"] and ev:
            refuse(f"{seg} is below the census floor but reports {ev} "
                   "detected events")
        r.setdefault("defects_removed_by_repair_not_excision", False)
        r["defect_bearing"] = bool(ev)
        r["voxel_scale_available"] = r["voxel_um"] is not None

    for r in rows.values():
        for k in [k for k in r if k.endswith("__from")]:
            del r[k]

    # ---- 5. the claims, DERIVED ----------------------------------------
    R = list(rows.values())
    def n(pred):
        return sum(1 for r in R if pred(r))

    claims = {
        "n_segments": len(R),
        "n_pinned": n(lambda r: r["inventory"] == "pinned"),
        "n_expansion": n(lambda r: r["inventory"] == "expansion"),
        "n_censusable": n(lambda r: r["censusable"]),
        "n_excluded_not_censusable": n(lambda r: not r["censusable"]),
        "n_transformed": n(lambda r: r["disposition"] == "transformed"),
        "n_already_clean": n(lambda r: r["disposition"] == "already_clean"),
        "n_base_original": n(lambda r: r["censusable"]
                             and r["base_kind"] == "original"),
        "n_base_displacement_repaired":
            n(lambda r: r["censusable"]
              and r["base_kind"] == "displacement_repaired"),
        "n_intrinsic_processed": n(lambda r: r["intrinsic_processed"]),
        "n_defect_bearing": n(lambda r: r["defect_bearing"]),
        "n_processed_without_events": n(lambda r: r["intrinsic_processed"]
                                        and not r["defect_bearing"]),
        "n_official_clean_positive_triangles":
            n(lambda r: r.get("official_clean") is True
              and r.get("official_triangles_positive") is True),
        "n_core_gate_failures": n(lambda r: r.get("core_gate_pass") is False),
        # Named because it is the one place "defect-bearing" and
        # "transformed" legitimately disagree.
        "n_defects_removed_by_repair_not_excision":
            n(lambda r: r["defects_removed_by_repair_not_excision"]),
        "n_processed_below_census_floor":
            n(lambda r: r["intrinsic_processed"] and not r["censusable"]),
        "defect_bearing_operand": "published original trace",
        # Named because neither is what a reader assumes. The contacts are
        # counted on the meshes actually cut -- 103 of which are this
        # project's repaired bases, not published traces -- and the scale
        # count is an INVENTORY property, not a denominator for events
        # with millimetre separations.
        "contacts_operand": "benchmark input base (the mesh actually cut)",
        "n_traces_with_voxel_scale_in_inventory":
            n(lambda r: r["voxel_scale_available"]),
        "voxel_scale_note": (
            "inventory-level availability: how many traces publish a voxel "
            "size. NOT the denominator for events carrying a millimetre "
            "separation, which is an event-level count."),
        "intrinsic_event_semantics": (
            "sum of independently clustered per-diagonal events; not a "
            "cross-diagonal union, not contact rows, and not a count of "
            "unique physical defects. The 190 defect-bearing count IS a "
            "trace-level union -- any event on either diagonal."),
        "n_contacts_total": sum(r["n_contacts"] or 0 for r in R),
        "core_gate_failures": sorted(
            (r["segment"] for r in R if r.get("core_gate_pass") is False),
            key=lambda s: rows[s]["retained_fraction"]),
    }

    # The identities. Each is an equation the prose relies on; a broken
    # one means a population sentence somewhere is now false.
    C = claims
    for lhs, rhs, what in (
            (C["n_pinned"] + C["n_expansion"], C["n_segments"],
             "pinned + expansion = segments"),
            (C["n_censusable"] + C["n_excluded_not_censusable"],
             C["n_segments"], "censusable + excluded = segments"),
            (C["n_transformed"] + C["n_already_clean"], C["n_censusable"],
             "transformed + already clean = censusable"),
            (C["n_base_original"] + C["n_base_displacement_repaired"],
             C["n_censusable"], "base kinds = censusable"),
            (C["n_defect_bearing"] + C["n_processed_without_events"],
             C["n_intrinsic_processed"],
             "defect-bearing + zero-event = processed"),
            (C["n_transformed"] + C["n_defects_removed_by_repair_not_excision"],
             C["n_defect_bearing"],
             "transformed + repaired-not-excised = defect-bearing"),
            (C["n_official_clean_positive_triangles"], C["n_censusable"],
             "officially clean over positive triangles = censusable")):
        if lhs != rhs:
            refuse(f"identity broken -- {what}: {lhs} != {rhs}")

    return {
        "schema": "windcheck/v1#reconciliation",
        "what": ("one row per segment; every population claim this project "
                 "makes is derived from these rows and none is retyped"),
        "sources": {k: {"path": str(p), "sha256": sha256(p)}
                    for k, p in paths.items()},
        "claims": claims,
        "rows": sorted(R, key=lambda r: r["segment"]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    for k, v in SOURCES.items():
        ap.add_argument(f"--{k.replace('_', '-')}", default=v)
    ap.add_argument("--out", default="out/reconciliation.json")
    a = ap.parse_args()
    paths = {k: Path(getattr(a, k)) for k in SOURCES}
    for k, p in paths.items():
        if not p.is_file():
            refuse(f"source {k} is missing: {p}")
    rec = build(paths)
    Path(a.out).write_text(json.dumps(rec, indent=1) + "\n")
    c = rec["claims"]
    print(f"wrote {a.out}: {c['n_segments']} segments")
    for k in ("n_pinned", "n_expansion", "n_censusable",
              "n_excluded_not_censusable", "n_transformed", "n_already_clean",
              "n_base_original", "n_base_displacement_repaired",
              "n_intrinsic_processed", "n_defect_bearing",
              "n_processed_without_events",
              "n_official_clean_positive_triangles", "n_core_gate_failures",
              "n_traces_with_voxel_scale_in_inventory", "n_contacts_total",
              "n_defects_removed_by_repair_not_excision",
              "n_processed_below_census_floor"):
        print(f"  {k:38s} {c[k]:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
