# SPDX-License-Identifier: Apache-2.0
"""S-9: score the pre-registered predictions (PREDICTION.md, F1-F5) against the run.

Reads the per-replicate table fits.parquet (export_fits.py, from study_bayes_coverage --out) and,
for the groups of D-136, tests/golden/recovery/summary.csv (group A: the scenarios whose PD Wald
verdict is DEFERRED; B, C, D: the rest at T = 20, 40, 100). Prints one Markdown row per prediction
(predicted, result, held or not) and, with --summary-out, writes one row per scenario, parameter,
prior, arm and interval type: coverage, raw class against the band, misses below and above the
truth, median logit width, refinements and refusals. Committed before the full run's results existed.

    uv run --no-project --with pyarrow==25.0.1 python studies/bayes-coverage/compare.py FITS_PARQUET \
        [--summary-out SUMMARY_CSV]

How each statistic is computed (fixed before the run; the addendum of 2026-09-28 in PREDICTION.md
states the same readings):
  - Coverage is over all replicates. A fit the estimator refuses (kPosteriorRefused) or without a
    finite posterior (kPosteriorNumeric) does not cover, as in D-131. The band is 0.95 +/- 3.29
    sqrt(0.95 x 0.05 / R). The profile coverage is the pinned one (not computed = not covering).
  - "Refines" is kPosteriorRefined set; "refuses" is kPosteriorRefused set.
  - F2 is scored for each prior: refinement pooled over the B-D replicates (>= 95%), refusal per B-D
    scenario (at most 1 of 1,000, i.e. <= 0.1%). Held if both priors hold.
  - F5 is scored for each prior, the median logit width over the replicates with an interval; held if
    both priors hold.
"""

import argparse
import csv
import math
from pathlib import Path

import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
BAND_Z = 3.29  # TOL_RECOVERY_COVERAGE_BAND_Z (tests/tolerances.toml)
REFINED, REFUSED, NUMERIC = 1 << 0, 1 << 1, 1 << 3  # engine::kPosterior* (engine/posterior.hpp)
NOT_COMPUTED = 1 << 2  # engine::kIntervalNotComputed (engine/profile.hpp)
ARMS = ("flat", "jeffreys", "flat_norule", "jeffreys_norule")
COARSE = (8, 32, 34, 35, 58, 59, 61, 62)


def logit(v):
    return math.log(v) - math.log1p(-v)


def median(x):
    x = sorted(x)
    n = len(x)
    if n == 0:
        return float("nan")
    return x[n // 2] if n % 2 else 0.5 * (x[n // 2 - 1] + x[n // 2])


def groups():
    out = {}
    with open(REPO / "tests/golden/recovery/summary.csv", newline="") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            sid = int(r["scenario"])
            T = [20, 40, 100][(sid // 3) % 3]
            out[sid] = "A" if r["pd_verdict"] == "DEFERRED" else {20: "B", 40: "C", 100: "D"}[T]
    return out


class Scenario:
    def __init__(self, rows):
        r0 = rows[0]
        self.id, self.T, self.n = r0["scenario"], r0["periods"], r0["obligors"]
        self.truth = {"pd": r0["pd_true"], "rho": r0["rho_true"]}
        self.rows = rows
        self.R = len(rows)
        half = BAND_Z * math.sqrt(0.95 * 0.05 / self.R)
        self.band = (0.95 - half, 0.95 + half)

    def ok(self, r, arm):
        return not r[f"{arm}_flags"] & (REFUSED | NUMERIC)

    def ends(self, r, arm, kind, p):
        return r[f"{arm}_{kind}_lo_{p}"], r[f"{arm}_{kind}_hi_{p}"]

    def coverage(self, arm, kind, p):
        """(coverage, misses below the truth, misses above) over all replicates."""
        cov = below = above = 0
        t = self.truth[p]
        for r in self.rows:
            if not self.ok(r, arm):
                continue
            lo, hi = self.ends(r, arm, kind, p)
            if lo <= t <= hi:
                cov += 1
            elif hi < t:
                below += 1
            elif lo > t:
                above += 1
        return cov / self.R, below, above

    def profile_coverage(self, p):
        t = self.truth[p]
        return sum(not r[f"prof_flags_{p}"] & NOT_COMPUTED and r[f"prof_lo_{p}"] <= t <= r[f"prof_hi_{p}"]
                   for r in self.rows) / self.R

    def width(self, arm, kind, p):
        w = []
        for r in self.rows:
            if self.ok(r, arm):
                lo, hi = self.ends(r, arm, kind, p)
                w.append(logit(hi) - logit(lo))
        return median(w)

    def share(self, arm, flag):
        return sum(bool(r[f"{arm}_flags"] & flag) for r in self.rows) / self.R

    def cls(self, c):
        return "BELOW" if c < self.band[0] else "ABOVE" if c > self.band[1] else "PASS"


def load(path):
    by = {}
    for r in pq.read_table(path).to_pylist():
        by.setdefault(r["scenario"], []).append(r)
    return {sid: Scenario(sorted(rows, key=lambda r: r["replicate"])) for sid, rows in sorted(by.items())}


def row(label, prediction, result, held):
    print(f"| {label} | {prediction} | {result} | {'yes' if held else 'no'} |")


def score(sc, grp):
    bd = [s for s in sc.values() if grp[s.id] in ("B", "C", "D")]
    print("| # | Prediction | Result | Held |")
    print("|---|---|---|---|")

    parts, held = [], True
    for arm in ("flat_norule", "jeffreys_norule"):
        coarse = [sc[i] for i in COARSE if i in sc]
        above = [s for s in coarse if s.cls(s.coverage(arm, "et", "pd")[0]) == "ABOVE"]
        parts.append(f"{arm.split('_')[0]} {len(above)} of {len(coarse)}")
        held = held and len(coarse) == 8 and len(above) >= 6
    row("F1", "rule off: equal-tailed PD above the band in >= 6 of the 8 coarse scenarios, both priors",
        "; ".join(parts), held)

    parts, held = [], True
    for arm in ("flat", "jeffreys"):
        refined = sum(s.share(arm, REFINED) * s.R for s in bd) / sum(s.R for s in bd)
        worst = max(bd, key=lambda s: s.share(arm, REFUSED))
        ok = refined >= 0.95 and all(s.share(arm, REFUSED) <= 0.001 for s in bd)
        parts.append(f"{arm}: refined {refined:.1%}, most refusals {worst.share(arm, REFUSED):.1%} (scenario {worst.id})")
        held = held and len(bd) == 50 and ok
    row("F2", "B-D: the rule refines in >= 95% of replicates and refuses <= 0.1% in every scenario",
        "; ".join(parts), held)

    diffs = []
    for s in bd:
        for p in ("pd", "rho"):
            diffs.append(abs(s.coverage("jeffreys", "et", p)[0] - s.profile_coverage(p)))
    k = sum(d <= 0.015 for d in diffs)
    row("F3", "B-D: Jeffreys equal-tailed coverage within 0.015 of the profile in >= 80 of 100 verdicts",
        f"{k} of {len(diffs)}; largest difference {max(diffs):.3f}", len(diffs) == 100 and k >= 80)

    b = [s for s in bd if grp[s.id] == "B"]
    k = sum(s.coverage("flat", "et", "rho")[0] >= s.coverage("jeffreys", "et", "rho")[0] for s in b)
    row("F4", "group B: flat equal-tailed rho coverage >= Jeffreys' in >= 9 of 12", f"{k} of {len(b)}",
        len(b) == 12 and k >= 9)

    parts, held = [], True
    for arm in ("flat", "jeffreys"):
        narrower = sum(s.width(arm, "hpd", "rho") < s.width(arm, "et", "rho") for s in bd)
        close = sum(abs(s.coverage(arm, "hpd", "rho")[0] - s.coverage(arm, "et", "rho")[0]) <= 0.02 for s in bd)
        parts.append(f"{arm}: HPD narrower {narrower} of {len(bd)}, coverage within 0.02 {close} of {len(bd)}")
        held = held and len(bd) == 50 and narrower >= 45 and close >= 40
    row("F5", "B-D: HPD median rho width below equal-tailed in >= 45 of 50, coverages within 0.02 in >= 40 of 50",
        "; ".join(parts), held)


def write_summary(sc, grp, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["scenario", "group", "param", "arm", "interval", "coverage", "class", "below", "above",
                    "median_logit_width", "refined_share", "refused_share", "profile_coverage"])
        for s in sc.values():
            for p in ("pd", "rho"):
                for arm in ARMS:
                    for kind in ("et", "hpd"):
                        c, lo, hi = s.coverage(arm, kind, p)
                        w.writerow([s.id, grp[s.id], p, arm, kind, f"{c:.4f}", s.cls(c), lo, hi,
                                    f"{s.width(arm, kind, p):.6g}", f"{s.share(arm, REFINED):.4f}",
                                    f"{s.share(arm, REFUSED):.4f}", f"{s.profile_coverage(p):.4f}"])
    print(f"\nwrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fits")
    ap.add_argument("--summary-out")
    a = ap.parse_args()
    sc = load(a.fits)
    grp = groups()
    numeric = {arm: sum(bool(r[f"{arm}_flags"] & NUMERIC) for s in sc.values() for r in s.rows) for arm in ARMS}
    print(f"{len(sc)} scenarios, {sum(s.R for s in sc.values())} panels; no finite posterior (counted as not "
          f"covering): {numeric}\n")
    score(sc, grp)
    if a.summary_out:
        write_summary(sc, grp, a.summary_out)


if __name__ == "__main__":
    main()
