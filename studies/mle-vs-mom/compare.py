# SPDX-License-Identifier: Apache-2.0
"""The estimator-comparison pass (S-8, S-27, S-28): score predictions E1-E13 against the run.

Reads the per-replicate table fits.parquet (export_fits.py, from study_estimator_pass --out) and,
for the groups of D-136, tests/golden/recovery/summary.csv (group A: the scenarios whose PD Wald
verdict is DEFERRED; B, C, D: the rest at T = 20, 40, 100). Prints one Markdown row per
prediction (predicted, result, held or not) and, with --summary-out, writes one row per scenario of
the reported statistics. Committed before the full run's results existed.

    uv run --no-project --with scipy==1.16.2 --with pyarrow==25.0.1 python studies/mle-vs-mom/compare.py \
        FITS_PARQUET [--summary-out SUMMARY_CSV]

How each statistic is computed (fixed before the run; the addendum of 2026-09-28 in PREDICTION.md
states the same readings):
  - An estimator's bias, RMSE and mean relative error of PD-hat are over the replicates it estimates:
    the binomial MLE all of them; MoM those without kMomRefused; a rate-MLE treatment those with
    status 0. The refusal share is 1 minus the estimated share.
  - A rate-MLE profile coverage is over the replicates the treatment estimates; an interval not
    computed does not cover (D-131). The band is 0.95 +/- 3.29 sqrt(0.95 x 0.05 / m), m = that count.
  - S-27's pairwise statistics (the mean difference of rho-hat, the coverage difference) are over the
    replicates both estimators estimate; the SE unit is the root mean square of the binomial fit's
    Hessian SE of rho over those replicates where it is finite.
  - The model references are computed at full precision with model_reference.py's own functions
    (P0 = p_zero; the S-27 bias in SE units = binomial_probit_variance / (s^2 sqrt(2/T))), not read
    from the rounded tables: a printed 0.05 can be either side of E8's threshold. E3's large-n
    RMSE ratio is the registered table's value (a Monte Carlo number; the +/-25% band dwarfs its
    rounding).
  - E9 covers the 81 scenarios with n <= 10^4, as registered. E10-E13 are conditioned on P0 or on the
    reference bias, not on n, so they range over all 135 scenarios that meet the condition.
  - E12 is held when it holds for each of the three treatments separately.
"""

import argparse
import csv
import math
import sys
from pathlib import Path

import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model_reference as ref  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
BAND_Z = 3.29  # TOL_RECOVERY_COVERAGE_BAND_Z (tests/tolerances.toml)
MOM_REFUSED = 1 << 0  # engine::kMomRefused (engine/moments.hpp)
TREATMENTS = ("refuse", "drop", "censor", "substitute")
PDS, RHOS, TS = (0.001, 0.01, 0.05), (0.02, 0.12, 0.24), (20, 40, 100)
INDISTINGUISHABLE_SE = 0.1
INDISTINGUISHABLE_COVERAGE = 0.0227

# The registered large-n RMSE ratio MoM/MLE for rho (PREDICTION.md, the S-8 table), by (PD, rho, T).
S8_REFERENCE = {
    (0.001, 0.02): (1.12, 1.18, 1.22), (0.001, 0.12): (1.44, 1.66, 1.96), (0.001, 0.24): (1.91, 2.20, 2.68),
    (0.01, 0.02): (1.08, 1.11, 1.13), (0.01, 0.12): (1.34, 1.49, 1.63), (0.01, 0.24): (1.63, 1.86, 2.12),
    (0.05, 0.02): (1.04, 1.06, 1.07), (0.05, 0.12): (1.23, 1.30, 1.33), (0.05, 0.24): (1.42, 1.55, 1.60),
}


def band(m):
    half = BAND_Z * math.sqrt(0.95 * 0.05 / m)
    return 0.95 - half, 0.95 + half


def mean(x):
    return sum(x) / len(x) if x else float("nan")


def rmse(errs):
    return math.sqrt(sum(e * e for e in errs) / len(errs)) if errs else float("nan")


def groups():
    """D-136's groups for the recovery scenarios 0-80."""
    out = {}
    with open(REPO / "tests/golden/recovery/summary.csv", newline="") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            sid = int(r["scenario"])
            T = [20, 40, 100][(sid // 3) % 3]
            out[sid] = "A" if r["pd_verdict"] == "DEFERRED" else {20: "B", 40: "C", 100: "D"}[T]
    return out


class Scenario:
    """One scenario's replicates and the statistics the predictions read."""

    def __init__(self, rows):
        r0 = rows[0]
        self.id, self.pd, self.rho = r0["scenario"], r0["pd_true"], r0["rho_true"]
        self.T, self.n = r0["periods"], r0["obligors"]
        self.R = len(rows)
        s2 = self.rho / (1 - self.rho)
        self.p0 = ref.p_zero(self.pd, self.rho, self.n)
        self.refused_expected = 1 - (1 - self.p0) ** self.T
        self.ref_bias_se = ref.binomial_probit_variance(self.pd, self.rho, self.n) / (s2 * math.sqrt(2 / self.T))
        self.est = {}  # estimator -> list of (row index) it estimates
        self.est["bin"] = list(range(self.R))
        self.est["mom"] = [i for i, r in enumerate(rows) if not r["mom_flags"] & MOM_REFUSED]
        for t in TREATMENTS:
            self.est[t] = [i for i, r in enumerate(rows) if r[f"{t}_status"] == 0]
        self.rows = rows

    def values(self, e, p, idx=None):
        key = {"bin": "bin", "mom": "mom"}.get(e, e)
        idx = self.est[e] if idx is None else idx
        return [self.rows[i][f"{key}_{p}"] for i in idx]

    def truth(self, p):
        return self.pd if p == "pd" else self.rho

    def bias(self, e, p, idx=None):
        return mean([v - self.truth(p) for v in self.values(e, p, idx)])

    def rmse(self, e, p):
        return rmse([v - self.truth(p) for v in self.values(e, p)])

    def rel_err_pd(self, e):
        return mean([(v - self.pd) / self.pd for v in self.values(e, "pd")])

    def share(self, e):
        return len(self.est[e]) / self.R

    def coverage(self, e, p, idx=None):
        idx = self.est[e] if idx is None else idx
        key = "bin" if e == "bin" else e
        return sum(self.rows[i][f"{key}_cover_{p}"] == 1 for i in idx) / len(idx) if idx else float("nan")

    def below_band(self, e, p):
        m = len(self.est[e])
        return m > 0 and self.coverage(e, p) < band(m)[0]

    def pairwise(self, e):
        """S-27: (mean rho difference in binomial SE units, rho coverage difference), on common replicates."""
        idx = self.est[e]
        if not idx:
            return float("nan"), float("nan")
        se = [self.rows[i]["bin_se_rho"] for i in idx if math.isfinite(self.rows[i]["bin_se_rho"])]
        unit = math.sqrt(mean([s * s for s in se]))
        diff = mean(self.values(e, "rho", idx)) - mean(self.values("bin", "rho", idx))
        return diff / unit, self.coverage(e, "rho", idx) - self.coverage("bin", "rho", idx)

    def indistinguishable(self, e="censor"):
        d, c = self.pairwise(e)
        return abs(d) <= INDISTINGUISHABLE_SE and abs(c) <= INDISTINGUISHABLE_COVERAGE


def load(path):
    t = pq.read_table(path).to_pylist()
    by = {}
    for r in t:
        by.setdefault(r["scenario"], []).append(r)
    return {sid: Scenario(sorted(rows, key=lambda r: r["replicate"])) for sid, rows in sorted(by.items())}


def row(label, prediction, result, held):
    print(f"| {label} | {prediction} | {result} | {'yes' if held else 'no'} |")


def score(sc, grp):
    bd = [s for s in sc.values() if grp.get(s.id) in ("B", "C", "D")]
    ratio_rho = {s.id: s.rmse("mom", "rho") / s.rmse("bin", "rho") for s in bd}
    print("| # | Prediction | Result | Held |")
    print("|---|---|---|---|")

    k = sum(ratio_rho[s.id] > 1 for s in bd)
    row("E1", "RMSE(MoM)/RMSE(MLE) for rho > 1 in >= 45 of the 50 B-D scenarios",
        f"{k} of {len(bd)}; ratios {min(ratio_rho.values()):.2f}-{max(ratio_rho.values()):.2f}", len(bd) == 50 and k >= 45)

    cells, inc = 0, 0
    for pd in PDS:
        for T in TS:
            for n in (100, 1000, 10000):
                trio = [s for s in bd if s.pd == pd and s.T == T and s.n == n]
                if len(trio) != 3:
                    continue
                trio.sort(key=lambda s: s.rho)
                cells += 1
                a, b, c = (ratio_rho[s.id] for s in trio)
                inc += a < b < c
    row("E2", "the rho ratio increases 0.02 -> 0.12 -> 0.24 in >= 80% of (PD, T, n) cells all in B-D",
        f"{inc} of {cells} ({inc / cells:.0%})" if cells else "no such cell", cells > 0 and inc >= 0.8 * cells)

    big = [s for s in bd if s.n == 10000]
    within = []
    for s in big:
        r_ref = S8_REFERENCE[(s.pd, s.rho)][TS.index(s.T)]
        within.append(abs(ratio_rho[s.id] / r_ref - 1) <= 0.25)
    row("E3", "at n = 10^4 (B-D), the rho ratio within +/-25% of the large-n reference in >= 75%",
        f"{sum(within)} of {len(big)} ({sum(within) / len(big):.0%})", sum(within) >= 0.75 * len(big))

    k = sum(s.bias("mom", "rho") < s.bias("bin", "rho") for s in bd)
    row("E4", "MoM's rho-hat mean error more negative than the MLE's in >= 40 of the 50 B-D scenarios",
        f"{k} of {len(bd)}", len(bd) == 50 and k >= 40)

    ratio_pd = [s.rmse("mom", "pd") / s.rmse("bin", "pd") for s in bd]
    k = sum(0.95 <= x <= 1.15 for x in ratio_pd)
    row("E5", "RMSE ratio for PD in 0.95-1.15 in >= 45 of the 50 B-D scenarios",
        f"{k} of {len(bd)}; ratios {min(ratio_pd):.3f}-{max(ratio_pd):.3f}", len(bd) == 50 and k >= 45)

    top = [s for s in sc.values() if s.n == 1_000_000]
    k = sum(s.indistinguishable() for s in top)
    row("E6", "n = 10^6: censored rate MLE indistinguishable from the binomial MLE in >= 24 of 27",
        f"{k} of {len(top)}", len(top) == 27 and k >= 24)

    e7 = [s for s in sc.values() if s.ref_bias_se >= 0.5 and s.share("censor") >= 0.5]
    up = [s for s in e7 if mean(s.values("censor", "rho")) > mean(s.values("bin", "rho"))]
    row("E7", "mean rho-hat (rate) > mean rho-hat (binomial) in every scenario with reference bias >= 0.5 SE",
        f"{len(up)} of {len(e7)}" + ("" if len(up) == len(e7) else
                                     f"; not: {', '.join(str(s.id) for s in e7 if s not in up)}"),
        len(e7) > 0 and len(up) == len(e7))

    e8 = [s for s in sc.values() if s.n >= 1000 and (s.ref_bias_se < 0.05 or s.ref_bias_se > 0.3)]
    agree = [s for s in e8 if s.indistinguishable() == (s.ref_bias_se < 0.1)]
    row("E8", "'indistinguishable' agrees with 'reference below 0.1' in >= 80% of unambiguous cells, n >= 10^3",
        f"{len(agree)} of {len(e8)} ({len(agree) / len(e8):.0%})", len(agree) >= 0.8 * len(e8))

    small = [s for s in sc.values() if s.n <= 10000]
    ok = []
    for s in small:
        q = s.refused_expected
        ok.append(abs((1 - s.share("refuse")) - q) <= 3.29 * math.sqrt(q * (1 - q) / s.R))
    row("E9", "parity refusal share within +/-3.29 sqrt(q(1-q)/R) of 1-(1-P0)^T in >= 79 of 81 (a bug check)",
        f"{sum(ok)} of {len(small)}" + ("" if all(ok) else
                                        f"; outside: {', '.join(str(s.id) for s, g in zip(small, ok) if not g)}"),
        len(small) == 81 and sum(ok) >= 79)

    e10 = [s for s in sc.values() if s.p0 >= 0.2 and s.share("drop") >= 0.5]
    k = [s for s in e10 if s.rel_err_pd("drop") > 0.10]
    row("E10", "P0 >= 0.2: drop's mean relative error of PD-hat > +10% in every scenario it estimates in >= half",
        f"{len(k)} of {len(e10)}" + (f"; range {min(s.rel_err_pd('drop') for s in e10):+.0%} to "
                                     f"{max(s.rel_err_pd('drop') for s in e10):+.0%}" if e10 else ""),
        len(e10) > 0 and len(k) == len(e10))

    e11 = [s for s in sc.values() if s.p0 >= 0.05 and s.share("censor") >= 0.5 and s.share("drop") >= 0.5]
    k = sum(abs(s.rel_err_pd("censor")) < abs(s.rel_err_pd("drop")) for s in e11)
    row("E11", "P0 >= 0.05: |mean rel. error of PD-hat| censored < drop in >= 90% of scenarios both estimate",
        f"{k} of {len(e11)} ({k / len(e11):.0%})" if e11 else "no scenario", len(e11) > 0 and k >= 0.9 * len(e11))

    parts, held = [], True
    for t in ("drop", "censor", "substitute"):
        e12 = [s for s in sc.values() if s.ref_bias_se >= 1 and s.share(t) >= 0.5]
        k = sum(s.below_band(t, "rho") for s in e12)
        parts.append(f"{t} {k} of {len(e12)}" + (f" ({k / len(e12):.0%})" if e12 else ""))
        held = held and len(e12) > 0 and k >= 0.8 * len(e12)
    row("E12", "reference bias >= 1 SE: rho profile coverage below the band in >= 80%, for each of drop, "
        "censored, substituted", "; ".join(parts), held)

    e13 = [s for s in sc.values() if s.p0 >= 0.2 and s.share("censor") >= 0.5 and s.share("substitute") >= 0.5]
    k = sum(abs(s.rel_err_pd("censor")) < abs(s.rel_err_pd("substitute")) for s in e13)
    row("E13", "P0 >= 0.2: |mean rel. error of PD-hat| censored < substituted in >= 80% of scenarios both estimate",
        f"{k} of {len(e13)} ({k / len(e13):.0%})" if e13 else "no scenario", len(e13) > 0 and k >= 0.8 * len(e13))


def write_summary(sc, grp, path):
    cols = ["scenario", "group", "pd", "rho", "T", "n", "p0", "refused_expected", "ref_bias_se"]
    for e in ("bin", "mom") + TREATMENTS:
        cols += [f"{e}_share", f"{e}_pd_bias", f"{e}_pd_rmse", f"{e}_pd_rel_err", f"{e}_rho_bias", f"{e}_rho_rmse"]
        if e != "mom":
            cols += [f"{e}_pd_coverage", f"{e}_rho_coverage"]
    for e in TREATMENTS:
        cols += [f"{e}_rho_diff_se", f"{e}_rho_coverage_diff"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(cols)
        for s in sc.values():
            out = [s.id, grp.get(s.id, "-"), s.pd, s.rho, s.T, s.n, f"{s.p0:.6g}", f"{s.refused_expected:.6g}",
                   f"{s.ref_bias_se:.6g}"]
            for e in ("bin", "mom") + TREATMENTS:
                if not s.est[e]:
                    out += [f"{s.share(e):.4f}"] + ["nan"] * (5 if e == "mom" else 7)
                    continue
                out += [f"{s.share(e):.4f}", f"{s.bias(e, 'pd'):.6g}", f"{s.rmse(e, 'pd'):.6g}",
                        f"{s.rel_err_pd(e):.6g}", f"{s.bias(e, 'rho'):.6g}", f"{s.rmse(e, 'rho'):.6g}"]
                if e != "mom":
                    out += [f"{s.coverage(e, 'pd'):.4f}", f"{s.coverage(e, 'rho'):.4f}"]
            for e in TREATMENTS:
                d, c = s.pairwise(e)
                out += [f"{d:.4f}", f"{c:.4f}"]
            w.writerow(out)
    print(f"\nwrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fits")
    ap.add_argument("--summary-out")
    a = ap.parse_args()
    sc = load(a.fits)
    grp = groups()
    counts = {g: sum(v == g for v in grp.values()) for g in "ABCD"}
    if counts != {"A": 31, "B": 12, "C": 17, "D": 21}:
        sys.exit(f"groups {counts} differ from D-136's A 31, B 12, C 17, D 21")
    errors = {t: sum(r[f"{t}_status"] == 2 for s in sc.values() for r in s.rows) for t in TREATMENTS}
    print(f"{len(sc)} scenarios, {sum(s.R for s in sc.values())} panels; rate-MLE engine errors (counted as not "
          f"estimated): {errors}\n")
    score(sc, grp)
    if a.summary_out:
        write_summary(sc, grp, a.summary_out)


if __name__ == "__main__":
    main()
