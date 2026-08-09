"""One page describing the benchmark, generated from its own records.

No new inference: every number is read from index.json, locators.json,
the acceptance/self-check/witness-binding records, or the release
indices. If a figure is not already established, it does not appear.

That rule is enforced rather than intended, because this file previously
broke it in five ways at once:

  a missing record       `load()` returned {} and the figure printed as
                         "?", so a card generated with no acceptance test
                         at all looked like a card with an unknown number
                         rather than an ungrounded one;
  an unknown base_kind   anything not spelled "original" was counted as
                         displacement-repaired, so a typo or an absent
                         field silently inflated the repaired count;
  the core gate          the card named the LOWEST-RETENTION trace and
                         asserted it "fails the 99.9%-area core gate".
                         The verdict is recorded per segment in the
                         acceptance record; the card never read it. SIX
                         of 274 fail that gate, and the card said one;
  a partial population   the acceptance record was read for its headline
                         counts without checking that its rows cover the
                         packaged segments exactly -- a record for 3 of
                         274 traces reads identically at the top;
  truthiness             `n_removed_cells` was tested for truth to
                         separate transformed from already-clean, so an
                         absent count and a genuine zero were the same
                         thing.

All five now refuse.
"""
from __future__ import annotations
import argparse, json, math, statistics
from pathlib import Path

# Every value `base_kind` is allowed to take. A base whose kind is not on
# this list is not classifiable, and a card that guesses is worse than no
# card.
BASE_KINDS = {"original": "published original traces",
              "displacement_repaired": "certified displacement-repaired bases"}

# Likewise for disposition: which reference required a cut is recorded,
# and is not to be inferred from whether a count is non-zero.
DISPOSITIONS = ("transformed", "already_clean")


def _bool(v, what: str) -> bool:
    """A JSON boolean, and nothing else. `bool("false")` is True."""
    if not isinstance(v, bool):
        raise SystemExit(f"REFUSING: {what} is {v!r}, which is not a boolean")
    return v


def _count(v, what: str) -> int:
    """A non-negative integer. Counts of cells, contacts and components
    are not 3.5 and are not -1; accepting either means the distribution
    row below was computed from something that is not a count."""
    if isinstance(v, bool) or not isinstance(v, int) or v < 0:
        raise SystemExit(f"REFUSING: {what} is {v!r}, which is not a "
                         "non-negative integer")
    return v


def _fraction(v, what: str) -> float:
    """A finite number in [0, 1]. NaN is a float and compares false
    against every bound, so it must be excluded explicitly."""
    if isinstance(v, bool) or not isinstance(v, (int, float)) \
            or not math.isfinite(v) or not 0.0 <= v <= 1.0:
        raise SystemExit(f"REFUSING: {what} is {v!r}, which is not a finite "
                         "fraction in [0, 1]")
    return float(v)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", default="out/benchmark")
    ap.add_argument("--indices", nargs="+",
                    default=["out/release/index.json",
                             "out/release/expansion_index.json"])
    ap.add_argument("--reconciliation", default="out/reconciliation.json",
                    help="the per-segment reconciliation. Population "
                         "figures are taken FROM it, and every figure the "
                         "card derives itself must equal its claim.")
    ap.add_argument("--out", default="out/benchmark/CARD.md")
    a = ap.parse_args()
    B = Path(a.benchmark)
    idx = json.loads((B / "index.json").read_text())
    ent: dict = {}
    for p in a.indices:
        ent.update(json.loads(Path(p).read_text())["segments"])
    scored = [t for t in idx["traces"] if t.get("scored")]
    exc = [t for t in idx["traces"] if t.get("excluded_reason")]
    segs = [t["segment"] for t in scored]
    if len(set(segs)) != len(segs):
        raise SystemExit("REFUSING: index.json lists a scored segment twice")

    # ---- evaluation input and disposition, refusing the unclassifiable --
    kinds: dict[str, int] = {k: 0 for k in BASE_KINDS}
    disp: dict[str, int] = {d: 0 for d in DISPOSITIONS}
    for t in scored:
        e = ent.get(t["segment"])
        if e is None:
            raise SystemExit(f"REFUSING: {t['segment']} is scored but appears "
                             "in none of the release indices")
        k = e.get("base_kind")
        if k not in BASE_KINDS:
            raise SystemExit(
                f"REFUSING: {t['segment']} has base_kind {k!r}, which is not "
                f"one of {sorted(BASE_KINDS)}. The card will not classify an "
                "unknown base as repaired.")
        kinds[k] += 1
        d = t.get("disposition")
        if d not in DISPOSITIONS:
            raise SystemExit(f"REFUSING: {t['segment']} has disposition "
                             f"{d!r}, not one of {list(DISPOSITIONS)}")
        disp[d] += 1
        _count(t.get("n_contacts"), f"{t['segment']} n_contacts")
        _count(t.get("n_removed_cells"), f"{t['segment']} n_removed_cells")
        _fraction(t.get("retained_fraction_certificate"),
                  f"{t['segment']} retained_fraction_certificate")

    contacts = sorted(t["n_contacts"] for t in scored)
    removed = sorted(t["n_removed_cells"] for t in scored
                     if t["disposition"] == "transformed")
    ret = sorted(t["retained_fraction_certificate"] for t in scored)

    def q(v, f):
        return f"{f(v):,.6g}" if v else "-"

    def load(name):
        # Fail closed. A verification record that is absent is not a
        # verification whose result is unknown; it is a card asserting a
        # release was checked when nothing checked it.
        p = B / name
        if not p.is_file():
            raise SystemExit(f"REFUSING: {p} is missing. The card reports "
                             "what was verified before release and will not "
                             "print a placeholder in place of a check.")
        return json.loads(p.read_text())

    def need(rec, name, *keys):
        out = []
        for k in keys:
            if rec.get(k) is None:
                raise SystemExit(f"REFUSING: {name} records no {k!r}")
            out.append(rec[k])
        return out

    def verdict(rec, name):
        v = rec.get("verdict")
        if v != "PASS":
            raise SystemExit(f"REFUSING: {name} records verdict {v!r}. The "
                             "card describes a release that passed its "
                             "checks; it does not narrate a failing one.")
        return v

    def covers(n, name):
        # A denominator smaller than the packaged population is a check
        # that ran on part of it. The headline count reads the same.
        if n != len(scored):
            raise SystemExit(f"REFUSING: {name} covers {n} traces, the "
                             f"package has {len(scored)}")
        return n

    b7, sc, wb, nc = (load("B7-acceptance.json"), load("selfcheck.json"),
                      load("witness-binding.json"),
                      load("negative-controls.json"))
    loc = load("locators.json").get("coverage") or {}

    # The controls are a FROZEN population: seven mutations and one
    # positive control. Checking only numerator == denominator accepts a
    # record declaring 0 of each and a PASS verdict.
    N_NEGATIVE_CONTROLS, N_POSITIVE_CONTROLS = 7, 1

    def full(num, den, name):
        """For a record with no `verdict` field, numerator == denominator
        IS the verdict. Checking only the denominator lets the card print
        `12/274` and finish successfully, under a heading that says these
        are the checks the release passed."""
        _count(num, f"{name} numerator")
        _count(den, f"{name} denominator")
        covers(den, name)
        if num != den:
            raise SystemExit(f"REFUSING: {name} reports {num}/{den}. The "
                             "card reports what was verified before "
                             "release; a partial result is not that.")
        return num, den

    b7_agree, b7_n = full(*need(b7, "B7-acceptance.json", "n_agree",
                                "n_traces"), "B7-acceptance.json")
    verdict(b7, "B7-acceptance.json")
    for k in ("n_disagree", "n_missing_report"):
        if b7.get(k) != 0:
            raise SystemExit(f"REFUSING: B7-acceptance.json reports {k} = "
                             f"{b7.get(k)!r}, not 0")
    sc_ok, sc_n = full(*need(sc, "selfcheck.json", "n_ok", "n_traces"),
                       "selfcheck.json")
    wb_id, wb_n = full(*need(wb, "witness-binding.json", "n_identical",
                             "n_traces"), "witness-binding.json")
    loc_in, loc_ref = need(loc, "locators.json coverage",
                           "inputs_located", "references_located")
    for v, what in ((loc_in, "locators.json inputs_located"),
                    (loc_ref, "locators.json references_located")):
        covers(_count(v, what), what)
    # Counted by polarity: seven mutations must be REJECTED and one
    # positive control must be ACCEPTED. "8/8" states neither.
    nc_neg_ok, nc_neg, nc_pos_ok, nc_pos = need(
        nc, "negative-controls.json", "n_negative_rejected", "n_negative",
        "n_positive_accepted", "n_positive")
    for got, want, frozen, what in (
            (nc_neg_ok, nc_neg, N_NEGATIVE_CONTROLS,
             "negative controls rejected"),
            (nc_pos_ok, nc_pos, N_POSITIVE_CONTROLS,
             "positive controls accepted")):
        _count(got, what); _count(want, f"{what} denominator")
        if got != want or want != frozen:
            raise SystemExit(f"REFUSING: {what} is {got}/{want}; the suite "
                             f"has {frozen} of them")
    verdict(nc, "negative-controls.json")

    # ---- the reconciliation is the authority on population --------------
    # Every population number this project states is derived once, in
    # bench/reconcile.py, from one row per segment. The card recomputes
    # them from its own inputs and REQUIRES agreement rather than picking
    # one: two records that disagree about how many traces there are is
    # exactly the condition that produced three published mis-statements.
    rp = Path(a.reconciliation)
    if not rp.is_file():
        raise SystemExit(f"REFUSING: {rp} is missing. Population figures "
                         "come from the reconciliation, not from this "
                         "file's own arithmetic. Run bench/reconcile.py.")
    C = json.loads(rp.read_text())["claims"]
    for mine, theirs, what in (
            (len(scored), "n_censusable", "scorable cases"),
            (len(exc), "n_excluded_not_censusable", "exclusions"),
            (len(scored) + len(exc), "n_segments", "indexed traces"),
            (disp["transformed"], "n_transformed", "transformed references"),
            (disp["already_clean"], "n_already_clean", "already-clean references"),
            (kinds["original"], "n_base_original", "original bases"),
            (kinds["displacement_repaired"], "n_base_displacement_repaired",
             "displacement-repaired bases"),
            (sum(contacts), "n_contacts_total", "total witness contacts")):
        if mine != C[theirs]:
            raise SystemExit(f"REFUSING: the card counts {mine} {what}, the "
                             f"reconciliation says {C[theirs]}")

    # ---- the aggregates must add up to the population -------------------
    # The card previously reprinted summary figures without checking that
    # they described the rows beneath them. Every one of these is a
    # sentence the card makes; each is now derived and cross-checked.
    S = idx["summary"]
    for got, want, what in (
            (len(scored) + len(exc), S["n_traces_indexed"],
             "scored + excluded vs n_traces_indexed"),
            (len(scored), S["n_packaged"], "scored vs n_packaged"),
            (len(exc), S["n_excluded_not_censusable"],
             "excluded vs n_excluded_not_censusable"),
            (sum(contacts), S["n_contacts_total"],
             "summed per-case contacts vs n_contacts_total"),
            (disp["transformed"] + disp["already_clean"], len(scored),
             "transformed + already clean vs scored"),
            (sum(kinds.values()), len(scored), "base kinds vs scored")):
        if got != want:
            raise SystemExit(f"REFUSING: {what} -- {got} != {want}")

    # ---- the core gate, read rather than assumed ------------------------
    # The acceptance record carries `core_gate_pass` per segment, from the
    # evaluator AND from the certificate. Naming the worst-retention trace
    # and asserting it is the gate failure is a different claim from the
    # one the evidence supports.
    rows = [r["segment"] for r in b7["results"]]
    if len(set(rows)) != len(rows):
        raise SystemExit("REFUSING: the acceptance record lists a segment "
                         "twice")
    if set(rows) != set(segs):
        missing = sorted(set(segs) - set(rows))
        extra = sorted(set(rows) - set(segs))
        raise SystemExit(f"REFUSING: the acceptance record does not cover the "
                         f"packaged population exactly; missing {missing[:3]} "
                         f"({len(missing)}), unexpected {extra[:3]} "
                         f"({len(extra)})")
    gate_rows = {}
    for r in b7["results"]:
        s = r["segment"]
        ev = _bool((r.get("evaluator") or {}).get("core_gate_pass"),
                   f"{s} evaluator core_gate_pass")
        ce = _bool((r.get("certificate") or {}).get("core_gate_pass"),
                   f"{s} certificate core_gate_pass")
        if ev != ce:
            raise SystemExit(f"REFUSING: {s} core-gate verdict disagrees -- "
                             f"evaluator {ev}, certificate {ce}")
        gate_rows[s] = ev
    by_ret = {t["segment"]: t["retained_fraction_certificate"] for t in scored}
    fails = sorted((s for s, ok in gate_rows.items() if not ok),
                   key=lambda s: by_ret[s])
    # The same authority applies to WHICH traces fail, not only how many:
    # the count was right in one half of the corpus and wrong overall, so
    # the identities are what must agree.
    if len(fails) != C["n_core_gate_failures"]:
        raise SystemExit(f"REFUSING: the card counts {len(fails)} core-gate "
                         f"failures, the reconciliation says "
                         f"{C['n_core_gate_failures']}")
    if sorted(fails) != sorted(C["core_gate_failures"]):
        raise SystemExit("REFUSING: the card and the reconciliation name "
                         "different core-gate failures")

    lines = [
        "# Benchmark card — surface-topology benchmark",
        "",
        "Generated from this package's own records, plus the per-segment",
        "reconciliation shipped beside them as `reconciliation.json`. Every",
        "population figure below is DERIVED from that table's one row per",
        "segment rather than retyped, and this card refuses to build if its",
        "own arithmetic disagrees with it.",
        "",
        "The reconciliation names all seven of its sources by sha256. Four",
        "of them — the two release indices and the pinned intrinsic spectra",
        "— are not re-hosted here, so the intrinsic-scale rows cannot be",
        "re-derived from this archive alone; the benchmark-side rows can, "
        "against",
        "`index.json` and `B7-acceptance.json`.",
        "",
        "No figure appears here that is not already established in those",
        "records.",
        "",
        "## Population",
        "",
        f"- **{len(scored)} scorable cases**, {len(exc)} excluded, "
        f"{idx['summary']['n_traces_indexed']} indexed.",
        "- Evaluation input: "
        + ", ".join(f"**{kinds[k]} {BASE_KINDS[k]}**" for k in BASE_KINDS)
        + ".",
        f"- Reference: **{disp['transformed']} transformed**, "
        f"**{disp['already_clean']} already clean** (no cut required).",
        "",
        "## Distributions",
        "",
        "| | median | p90 | max |",
        "|---|---|---|---|",
        f"| witness contacts per case | {q(contacts, statistics.median)} | "
        f"{contacts[int(.9*len(contacts))]:,} | {max(contacts):,} |",
        f"| removed cells (transformed cases) | {q(removed, statistics.median)} | "
        f"{removed[int(.9*len(removed))]:,} | {max(removed):,} |",
        f"| retained fraction | {q(ret, statistics.median)} | "
        f"{ret[int(.1*len(ret))]:.6f} (p10) | {min(ret):.6f} (min) |",
        "",
        f"Total witness contacts shipped: **{idx['summary']['n_contacts_total']:,}**.",
        "",
        "## Cost outliers, named rather than averaged away",
        "",
        f"**{len(fails)} of {len(gate_rows)}** references are transverse-clean "
        "but fail the preregistered 99.9%-area core fragmentation gate. Every "
        "one is named:",
        "",
    ]
    ev_by_seg = {r["segment"]: r["evaluator"] for r in b7["results"]}
    for s in fails:
        ev = ev_by_seg[s]
        lines.append(
            f"- `{s}` — retains {by_ret[s]:.3f}, min core `R_main` "
            f"{_fraction(ev.get('min_R_main_core'), f'{s} min_R_main_core'):.3f} "
            f"over {_count(ev.get('n_core_components'), f'{s} n_core_components')} "
            "core components")
    lines += [
        "",
        "## Exclusions",
        "",
        f"{len(exc)} traces are excluded as not censusable, each with its "
        "valid-cell count and confirmation the census declined it:",
        "",
    ]
    for t in exc:
        # The two inventories record the count under different keys, and a
        # card printing "? valid cells" is worse than no card: it looks
        # like the evidence is missing when it is present.
        ev = t.get("excluded_evidence") or {}
        n = ev.get("n_valid", ev.get("n_valid_vertices"))
        if n is None:
            raise SystemExit(f"REFUSING: {t['segment']} has exclusion "
                             f"evidence with no valid-cell count: {ev}")
        # The heading says "confirmation the census declined it". The two
        # inventories record that confirmation differently, and a
        # valid-cell count alone is not it.
        inv = t.get("inventory")
        if inv == "expansion":
            if ev.get("decline_confirmed") is not True:
                raise SystemExit(f"REFUSING: {t['segment']} is excluded with "
                                 "no confirmed census decline")
            extra = "; census decline confirmed"
        elif inv == "pinned":
            if ev.get("n_retained_quads") != 0:
                raise SystemExit(f"REFUSING: {t['segment']} is excluded but "
                                 f"retains {ev.get('n_retained_quads')!r} "
                                 "quads")
            extra = "; 0 retained quads"
        else:
            raise SystemExit(f"REFUSING: {t['segment']} has inventory "
                             f"{inv!r}, so its exclusion evidence cannot be "
                             "checked")
        lines.append(f"- `{t['segment']}` — {n:,} valid cells{extra}")
    lines += [
        "",
        "## The evaluation contract",
        "",
        "- `clean` is **primary and binary**; an unclean candidate is not "
        "scored on cost at all.",
        "- A census over **zero triangles** is refused, not reported clean.",
        "- `retained_fraction` is priced on the **input**, so added or "
        "stretched geometry cannot inflate it.",
        "- **No single scalar score**, ever.",
        "- The reference is **one admissible answer, not the target**.",
        "- Frozen census parameters: `exclude=1`, `maxedge=60`, `cell=40`, "
        "`touch_tolerance=0.001`, both diagonals.",
        "",
        "## What was verified before release",
        "",
        f"- Acceptance vs certificates: **{b7_agree}/{b7_n}** — PASS",
        f"- Controls: **{nc_neg_ok}/{nc_neg} negative rejected; "
        f"{nc_pos_ok}/{nc_pos} positive accepted** — PASS",
        f"- Witness bound to a fresh census by identity: **{wb_id}/{wb_n}**",
        f"- Packer-independent verification: **{sc_ok}/{sc_n}**",
        f"- Operand locators: **{loc_in} inputs, {loc_ref} references**",
        "",
        "## Start here",
        "",
        "```sh",
        "./tools/quickstart.sh",
        "```",
    ]
    Path(a.out).write_text("\n".join(lines) + "\n")
    print(f"wrote {a.out} ({len(lines)} lines)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
