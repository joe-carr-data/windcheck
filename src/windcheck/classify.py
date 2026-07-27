"""The single source of truth for how a result is described.

Six different thresholds had accumulated across this codebase -- 0.25 and 1.4 in
the certificate, 0.15 and 1.6 in the command-line tool and the bundler, 1.8 in
the written documentation, and 0.094/0.823 as the observed gap. They disagreed,
so one segment could receive two different verdicts depending on which entry
point produced it. That is worse than any of the thresholds being wrong.

Everything that names a result now goes through this module.

## What is reported, and what is deliberately not

Three independent facts, kept separate:

    crossing_status   none | present
    period_status     agreed | disagreed | unavailable
    separation        a continuous number of revolutions, or None

**No semantic band.** Earlier versions sorted segments into "local", "one
revolution" and "wrap-scale", the last carrying an implied cause. Two problems:
adding a fifth scroll closed the gap between the upper two, so the three-way
split was an artifact of the four corpora it was derived from; and the causal
reading was never established -- three separate attempts to connect crossings to
sheet misassignment failed. The measurement supports a distance. It does not
support a diagnosis, and the vocabulary should not imply one.

Filters over separation are still useful, and are labelled literally --
"separation >= 0.15 revolutions" -- never "sheet switch".

## Why period_status exists

Separation is a distance divided by an estimated revolution period. When the two
independent estimators of that period disagree, the ratio is not interpretable,
and reporting it anyway would present a number whose denominator is in dispute.
Such segments are reported as intersecting with the scale unavailable, which is
the honest statement, rather than being silently placed in a band.
"""
from __future__ import annotations

# Literal filter boundaries. These name ranges of a measured quantity; they are
# NOT claims about what a segment is. Chosen as round numbers straddling the
# largest observed gap (0.094 -> 0.823 revolutions) rather than tuned.
SEP_NEARBY = 0.15
SEP_WIDE = 1.6

# Two geometric estimates of the same period should agree within a quarter.
# Fixed before it was known what it excluded.
PERIOD_RATIO_LO = 0.8
PERIOD_RATIO_HI = 1.25


def crossing_status(n_pairs: int) -> str:
    return "present" if n_pairs > 0 else "none"


def period_status(turning: float | None, neighbour: float | None,
                  angle_consistent: bool = True) -> str:
    """Whether the revolution period is trustworthy enough to divide by.

    `unavailable` and `disagreed` are kept distinct on purpose: the first means
    the question could not be asked -- a single winding has no neighbouring wrap
    inside itself to measure against, which is a property of the segment and not
    a fault -- and the second means it was asked and the answers conflicted.
    """
    if turning is None or turning != turning:
        return "unavailable"
    if not angle_consistent:
        return "disagreed"
    if neighbour is None or neighbour != neighbour or neighbour <= 0:
        return "unavailable"
    r = turning / neighbour
    return "agreed" if PERIOD_RATIO_LO <= r <= PERIOD_RATIO_HI else "disagreed"


def separation_band(sep: float | None) -> str:
    """A literal range label, for filtering and grouping only."""
    if sep is None:
        return "no crossing"
    if sep < SEP_NEARBY:
        return f"< {SEP_NEARBY} rev"
    if sep < SEP_WIDE:
        return f"{SEP_NEARBY} - {SEP_WIDE} rev"
    return f">= {SEP_WIDE} rev"


def verdict(crossing: str, period: str, sep: float | None,
            covering_span: float | None = None) -> str:
    """One sentence, stating only what was measured."""
    if crossing == "none":
        return "no self-intersection detected"
    if period != "agreed" or sep is None:
        return ("self-intersection present; revolution-scale classification "
                "unavailable because the two period estimates "
                + ("disagree" if period == "disagreed" else "could not both be made"))
    s = (f"self-intersection present; widest separation {sep:.2f} revolutions "
         f"along the trace's own parameter")
    if covering_span is not None:
        s += f", in a segment covering {covering_span:.2f}"
    return s


def describe(sep: float | None, covering_span: float | None) -> str:
    """A neutral note for a reader, with no causal claim.

    The one interpretive statement kept is geometric and checkable: a segment
    covering more than a full turn has two ends in the same angular sector, so
    it can meet itself there without anything having gone wrong. That is a fact
    about the parameterisation, not a diagnosis of the trace.
    """
    if sep is None:
        return "No transverse self-intersection was found in this surface."
    if (covering_span is not None and covering_span > 1.005
            and 0.75 <= sep <= 1.35):
        return ("This segment covers more than one full turn, so its two ends "
                "occupy the same angular sector and can meet there without any "
                "error. The measurement does not distinguish that case from a "
                "trace that has returned to a wrap it already traced.")
    if sep >= SEP_WIDE:
        return ("Two parts of this trace, more than "
                f"{SEP_WIDE} revolutions apart along its own parameter, occupy "
                "the same place. What caused that is not established here.")
    return ("A short-range self-contact. These are common across the corpus and "
            "are associated with elevated quad twist.")
