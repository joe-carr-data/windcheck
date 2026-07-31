"""Read the normal profiles and decide: tight packing, or doubled back?

Statistic fixed here before looking at the arrays.

Along the surface normal, papyrus is bright and the gap between wraps is dark.
For a correctly traced point the next sheet sits at the normal spacing, ~17 vx
(calibration: 17.05). So:

    NEAR band   |offset| in [4, 10] vx   -- empty for a normal wrap spacing
    FAR band    |offset| in [13, 21] vx  -- where the neighbouring sheet lives

`near_excess` = mean(NEAR) - mean(FAR), per profile, after normalising each
profile by its own central peak so overall brightness cancels.

  tight packing  a real second sheet sits a few vx away, so flagged points have
                 a BRIGHT near band -> near_excess markedly higher than controls
  doubled back   the trace visited one sheet twice, nothing is there, so the
                 near band stays dark -> near_excess indistinguishable

Controls come from the same boxes, so local packing and contrast cancel.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("out/routeA")
REACH = 26


def load():
    f = np.load(OUT / "normal_profile.flagged.npy")
    u = np.load(OUT / "normal_profile.unflagged.npy")
    return f.astype(np.float64), u.astype(np.float64)


MIN_PEAK = 20.0     # below this the ray sampled masked-out void, not papyrus


def norm(prof, offs):
    """Normalise by the central peak, after dropping rays that never hit sheet.

    Without the cut, one flagged profile with a central peak of exactly 0.0
    divided through to ~65,000 and single-handedly produced a "TIGHT PACKING"
    verdict with a mean of +65 on a quantity that should be O(1). The cut is
    applied identically to both arms.
    """
    c = np.abs(offs) <= 2
    peak = prof[:, c].mean(1)
    keep = peak >= MIN_PEAK
    return prof[keep] / peak[keep, None], keep


def stats(prof, offs):
    near = (np.abs(offs) >= 4) & (np.abs(offs) <= 10)
    far = (np.abs(offs) >= 13) & (np.abs(offs) <= 21)
    return prof[:, near].mean(1) - prof[:, far].mean(1)


def main():
    f, u = load()
    offs = np.arange(-REACH, REACH + 1, 1.0)
    print(f"flagged {f.shape}  unflagged {u.shape}")

    fn, kf = norm(f, offs)
    un, ku = norm(u, offs)
    print(f"dropped {(~kf).sum()} flagged and {(~ku).sum()} unflagged rays "
          f"with central peak < {MIN_PEAK} (sampled void)")
    sf, su = stats(fn, offs), stats(un, offs)

    # bootstrap CI on the difference of means
    rng = np.random.default_rng(0)
    boot = np.array([
        rng.choice(sf, len(sf)).mean() - rng.choice(su, len(su)).mean()
        for _ in range(4000)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    diff = sf.mean() - su.mean()

    print(f"\nnear_excess  flagged   {sf.mean():+.4f} (sd {sf.std():.4f}) "
          f"median {np.median(sf):+.4f}")
    print(f"near_excess  unflagged {su.mean():+.4f} (sd {su.std():.4f}) "
          f"median {np.median(su):+.4f}")
    print(f"median difference      {np.median(sf) - np.median(su):+.4f}"
          "   (robustness check)")
    print(f"difference             {diff:+.4f}   95% CI [{lo:+.4f}, {hi:+.4f}]")
    verdict = ("TIGHT PACKING: flagged points have a real near sheet"
               if lo > 0 else
               "no near-sheet excess: consistent with doubling back"
               if hi < 0 else
               "inconclusive: CI spans zero")
    print(f"verdict: {verdict}")

    mf, mu = fn.mean(0), un.mean(0)
    plt.figure(figsize=(11, 6))
    plt.plot(offs, mf, lw=2, label=f"flagged (n={len(fn)})", color="#d62728")
    plt.plot(offs, mu, lw=2, label=f"unflagged, same boxes (n={len(un)})",
             color="#1f77b4")
    plt.axvspan(4, 10, alpha=.10, color="red")
    plt.axvspan(-10, -4, alpha=.10, color="red")
    plt.axvspan(13, 21, alpha=.10, color="blue")
    plt.axvspan(-21, -13, alpha=.10, color="blue")
    plt.axvline(0, color="k", lw=.6)
    plt.xlabel("offset along surface normal (vx)   [7.91 um/vx]")
    plt.ylabel("CT intensity, normalised to the central peak")
    plt.title("Is there a second sheet where the trace doubles?\n"
              "red band = NEAR (4-10 vx), blue band = FAR (13-21 vx, "
              "normal wrap spacing 17.05)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "normal_profile.png", dpi=115)

    (OUT / "normal_verdict.json").write_text(json.dumps({
        "n_flagged": int(len(fn)), "n_unflagged": int(len(un)),
        "near_excess_flagged": float(sf.mean()),
        "near_excess_unflagged": float(su.mean()),
        "difference": float(diff), "ci95": [float(lo), float(hi)],
        "verdict": verdict}, indent=2))
    print("wrote", OUT / "normal_profile.png")


if __name__ == "__main__":
    main()
