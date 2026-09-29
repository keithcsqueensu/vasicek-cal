# SPDX-License-Identifier: Apache-2.0
"""S-15: score the pre-registered predictions (PREDICTION.md, G1-G5) against the run.

Reads the per-draw table draws.parquet (export_draws.py, from study_bayes_sbc --out). For each
setting, prior, arm (with the resolution rule, or the parity grid without it) and parameter: the
20-bin histogram of the rank statistic (the marginal posterior CDF at the true value), its chi-squared
p-value against uniform (19 degrees of freedom), and the tail share, the fraction of ranks below
0.025 or above 0.975. Prints the histograms and one Markdown row per prediction. Committed before
the run's results existed.

    uv run --no-project --with scipy==1.16.2 --with pyarrow==25.0.1 python studies/bayes-sbc/compare.py DRAWS_PARQUET

Readings fixed before the run (the addendum of 2026-09-28 in PREDICTION.md states the same): a
draw the estimator refuses has no rank; it is left out of that arm's histogram and tail share, and
the refusals are reported. The tail-share range 0.027-0.073 is inclusive; "below 0.027" is strict.
G5 compares tail shares under the flat prior (the only prior run at (10^4, 40) and (1,000, 100)) and
is held when (10^4, 40) lies strictly between the other two, in either order.
"""

import math
import sys

import pyarrow.parquet as pq
from scipy.stats import chi2

BINS = 20
SETTINGS = {0: "(1,000, 20)", 1: "(10^4, 40)", 2: "(1,000, 100)"}


def stats(ranks):
    ranks = [r for r in ranks if not math.isnan(r)]
    n = len(ranks)
    counts = [0] * BINS
    for r in ranks:
        counts[min(int(r * BINS), BINS - 1)] += 1
    e = n / BINS
    x2 = sum((c - e) ** 2 / e for c in counts)
    return {"n": n, "counts": counts, "p": float(chi2.sf(x2, BINS - 1)),
            "tail": sum(r < 0.025 or r > 0.975 for r in ranks) / n}


def calibrated(s):
    return s["p"] >= 0.01 and 0.027 <= s["tail"] <= 0.073


def main():
    t = pq.read_table(sys.argv[1]).to_pylist()
    res = {}
    for setting in SETTINGS:
        for prior in ("flat", "jeffreys"):
            rows = [r for r in t if r["setting"] == setting and r["prior"] == prior]
            if not rows:
                continue
            for arm, suffix in (("rule", ""), ("norule", "_norule")):
                for p in ("pd", "rho"):
                    s = stats([r[f"rank_{p}{suffix}"] for r in rows])
                    s["draws"] = len(rows)
                    res[(setting, prior, arm, p)] = s
                    print(f"{SETTINGS[setting]:>13} {prior:8} {arm:6} {p:3} n {s['n']:4} (refused "
                          f"{len(rows) - s['n']}) p {s['p']:.3g} tail {s['tail']:.3f} | {' '.join(map(str, s['counts']))}")
    print()
    print("| # | Prediction | Result | Held |")
    print("|---|---|---|---|")

    def cal_row(label, text, key):
        a, b = res[key + ("pd",)], res[key + ("rho",)]
        result = f"PD p {a['p']:.3g}, tail {a['tail']:.3f}; rho p {b['p']:.3g}, tail {b['tail']:.3f}"
        print(f"| {label} | {text} | {result} | {'yes' if calibrated(a) and calibrated(b) else 'no'} |")

    cal_row("G1", "(1,000, 20), rule, flat: p >= 0.01 and tail 0.027-0.073 for PD and rho", (0, "flat", "rule"))
    cal_row("G2", "(1,000, 20), rule, Jeffreys: the same", (0, "jeffreys", "rule"))
    a, b = res[(2, "flat", "norule", "pd")], res[(2, "flat", "norule", "rho")]
    print(f"| G3 | (1,000, 100), no rule: tail share below 0.027 for rho and PD (a hump) | PD {a['tail']:.3f}, "
          f"rho {b['tail']:.3f} | {'yes' if a['tail'] < 0.027 and b['tail'] < 0.027 else 'no'} |")
    cal_row("G4", "(1,000, 100), rule: p >= 0.01 and tail 0.027-0.073 for PD and rho", (2, "flat", "rule"))
    lo, mid, hi = (res[(s, "flat", "norule", "rho")]["tail"] for s in (0, 1, 2))
    held = min(lo, hi) < mid < max(lo, hi)
    print(f"| G5 | (10^4, 40), no rule: rho tail share between those of (1,000, 20) and (1,000, 100) | {mid:.3f} "
          f"against {lo:.3f} and {hi:.3f} | {'yes' if held else 'no'} |")


if __name__ == "__main__":
    main()
