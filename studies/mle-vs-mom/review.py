# SPDX-License-Identifier: Apache-2.0
"""The estimator-comparison pass: review the rate MLE's out-of-band coverage verdicts (D-155's practice).

Reads summary.csv (compare.py --summary-out) and writes reviewed.csv: one row per verdict (scenario,
treatment, parameter) whose profile coverage is outside the band, with its raw class and a diagnosis.
The diagnoses are assigned by the rules below, in order, from each verdict's own statistics; a verdict
no rule explains is labelled "unexplained" and listed, never forced into a class. Written after the
run, as a review (the predictions were scored by compare.py, committed before it).

  shared_small_t          the binomial MLE's own profile coverage of the parameter is below the band
                          on the same panels: the dip belongs to the panels (short T), not to the
                          rate likelihood (as for D-164's rate recovery).
  noise_as_factor_variance  rho-hat biased up by at least 10% (relative): the rate model reads the
                          binomial noise in d/n as factor variance (S-27's mechanism). For PD, the same
                          scenario's rho-hat is biased up that much and PD-hat is too (PD-hat follows
                          the widened latent distribution through Phi(mu / sqrt(1 + v))).
  refusal_selection       parity (refuse), rho-hat biased down by at least 10%: the panels it accepts
                          are those with no zero-default period, which selects low dispersion.
  low_tail_lost           drop or substitute, rho-hat biased down by at least 10%: dropping removes the
                          low tail of the rates, substitution piles it on one value (1/(2n)); either
                          shrinks the variance the model reads as rho.
  drop_keeps_worse_periods  drop, PD-hat biased up by at least 10%: the periods kept are the worse ones.
  substituted_rate_off    substitute, PD-hat biased by at least 10% either way: every zero period is put
                          at 1/(2n), not where its rate was.
  narrow_without_bias     coverage at least 0.88 with |relative bias| under 10%: a small shortfall with
                          no systematic error in the estimate (the interval is slightly narrow).
  conservative            above the band.

    uv run --no-project python studies/mle-vs-mom/review.py SUMMARY_CSV OUT_CSV
"""

import csv
import math
import sys

BAND_Z = 3.29
TREATMENTS = ("refuse", "drop", "censor", "substitute")
REL = 0.10


def f(r, k):
    return float(r[k])


def band(share):
    m = round(share * 1000)
    half = BAND_Z * math.sqrt(0.95 * 0.05 / m)
    return 0.95 - half, 0.95 + half


def diagnose(r, t, p, cls):
    if cls == "ABOVE":
        return "conservative"
    truth = f(r, "pd") if p == "pd" else f(r, "rho")
    rel = f(r, f"{t}_{p}_bias") / truth
    rho_rel = f(r, f"{t}_rho_bias") / f(r, "rho")
    if f(r, f"bin_{p}_coverage") < band(1.0)[0]:
        return "shared_small_t"
    if p == "rho":
        if rel >= REL:
            return "noise_as_factor_variance"
        if rel <= -REL:
            return "refusal_selection" if t == "refuse" else "low_tail_lost" if t in ("drop", "substitute") else None
    else:
        if t == "drop" and rel >= REL:
            return "drop_keeps_worse_periods"
        if t == "substitute" and abs(rel) >= REL:
            return "substituted_rate_off"
        if rel >= REL and rho_rel >= REL:
            return "noise_as_factor_variance"
    if abs(rel) < REL and f(r, f"{t}_{p}_coverage") >= 0.88:
        return "narrow_without_bias"
    return None


def main():
    rows = list(csv.DictReader(open(sys.argv[1], newline="")))
    out, unexplained = [], []
    for r in rows:
        for t in TREATMENTS:
            share = f(r, f"{t}_share")
            if share == 0:
                continue
            lo, hi = band(share)
            for p in ("pd", "rho"):
                c = f(r, f"{t}_{p}_coverage")
                if lo <= c <= hi:
                    continue
                cls = "BELOW" if c < lo else "ABOVE"
                d = diagnose(r, t, p, cls)
                truth = f(r, "pd") if p == "pd" else f(r, "rho")
                rec = [r["scenario"], t, p, f"{c:.3f}", cls, f"{f(r, f'{t}_{p}_bias') / truth:+.3f}",
                       d or "unexplained"]
                out.append(rec)
                if d is None:
                    unexplained.append(rec)
    with open(sys.argv[2], "w", newline="") as fh:
        fh.write("# Reviewed out-of-band coverage verdicts of the rate MLE in the estimator-comparison pass (S-27,\n"
                 "# S-28), one row per verdict: raw class and diagnosis, assigned by review.py's rules (its\n"
                 "# docstring explains each diagnosis). Coverage over the replicates the treatment estimates;\n"
                 "# rel_bias is the mean error of the estimate over the true value.\n")
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["scenario", "treatment", "param", "coverage", "class", "rel_bias", "diagnosis"])
        w.writerows(out)
    counts = {}
    for rec in out:
        counts[rec[6]] = counts.get(rec[6], 0) + 1
    print(f"{len(out)} verdicts outside the band: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    for rec in unexplained:
        print("unexplained:", rec)


if __name__ == "__main__":
    main()
