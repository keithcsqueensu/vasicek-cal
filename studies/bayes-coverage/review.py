# SPDX-License-Identifier: Apache-2.0
"""S-9: review the estimator's out-of-band coverage verdicts (D-155's practice).

Reads summary.csv (compare.py --summary-out) and writes reviewed.csv: one row per out-of-band verdict of
the estimator (the resolution rule on; flat and Jeffreys; equal-tailed and HPD; PD and rho), with its
raw class and a diagnosis assigned by the rules below, in order. A verdict no rule explains is labelled
"unexplained" and listed. The diagnostic arm (rule off) is not reviewed: its failures are what it was
registered to show. Written after the run, as a review; the predictions were scored by compare.py.

  shared_with_profile   below the band, and the pinned profile interval is below it too on the same
                        panels: the panels' finding (short T), not the posterior's.
  prior_pulls_up        below the band with at least 90% of the misses entirely above the truth: the
                        prior's mass drags the interval up where the data are thin (the flat prior at
                        PD 0.1% or rho 0.02; HPD, the shortest set in logit, the same way near the floor).
  small_t               below the band at T = 20 with misses on both sides of the truth.
  conservative          above the band in group A (near-uninformative data; the interval spans much of
                        the box).

    uv run --no-project python studies/bayes-coverage/review.py SUMMARY_CSV OUT_CSV
"""

import csv
import math
import sys

R = 1000
HALF = 3.29 * math.sqrt(0.95 * 0.05 / R)


def diagnose(r):
    prof = float(r["profile_coverage"])
    below, above = int(r["below"]), int(r["above"])
    T = [20, 40, 100][(int(r["scenario"]) // 3) % 3]
    if r["class"] == "ABOVE":
        return "conservative" if r["group"] == "A" else None
    if prof < 0.95 - HALF:
        return "shared_with_profile"
    if above >= 0.9 * (below + above):
        return "prior_pulls_up"
    if T == 20 and below > 0 and above > 0:
        return "small_t"
    return None


def main():
    rows = [r for r in csv.DictReader(open(sys.argv[1], newline="")) if r["arm"] in ("flat", "jeffreys")
            and r["class"] != "PASS"]
    out, counts = [], {}
    for r in rows:
        d = diagnose(r) or "unexplained"
        counts[d] = counts.get(d, 0) + 1
        out.append([r["scenario"], r["group"], r["arm"], r["interval"], r["param"], r["coverage"], r["class"],
                    r["below"], r["above"], d])
        if d == "unexplained":
            print("unexplained:", out[-1])
    with open(sys.argv[2], "w", newline="") as fh:
        fh.write("# Reviewed out-of-band coverage verdicts of the grid-Bayesian estimator (S-9; rule on), one row\n"
                 "# per verdict, with the diagnosis assigned by review.py's rules (its docstring explains each).\n"
                 "# below / above: misses with the interval entirely below / above the true value.\n")
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["scenario", "group", "prior", "interval", "param", "coverage", "class", "misses_below",
                    "misses_above", "diagnosis"])
        w.writerows(out)
    print(f"{len(out)} verdicts outside the band: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))


if __name__ == "__main__":
    main()
