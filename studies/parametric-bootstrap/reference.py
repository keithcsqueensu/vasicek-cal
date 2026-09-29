# SPDX-License-Identifier: Apache-2.0
"""S-10 and S-4 (the parametric-bootstrap run): reference numbers from pinned results and the model only.

No parametric-bootstrap panel is simulated and no likelihood-ratio statistic at the truth is computed;
the inputs are the pinned recovery summary (tests/golden/recovery/summary.csv: profile and iid bootstrap
coverage, the bias and spread of logit rho-hat) and S-9's per-replicate rows (Jeffreys equal-tailed
coverage on the same panels, D-175). Three parts:

  bartlett  If the profile LR statistic W at the truth is k chi2_1, the profile interval covers with
            probability F(3.8415 / k), F the chi2_1 CDF. Inverting each scenario's pinned coverage gives
            the factor k it implies, the factor an oracle Bartlett correction (S-4a) would apply.
  shift     A percentile interval from a parametric bootstrap at theta-hat is centred near theta-hat
            + b, b the estimator's bias, so it misses the truth as if the bias were 2b. In logit units
            with sd s: coverage ~ Phi(1.96 + 2b/s) - Phi(-1.96 + 2b/s), before skew and SE variability.
  gap       Per group and T: profile and Jeffreys coverage on the same replicates (S-9's rows), the
            targets S-10 and S-4b are compared with.

    uv run --no-project --with scipy==1.16.2 --with pyarrow==25.0.1 python studies/parametric-bootstrap/reference.py
"""

import csv
import math
from pathlib import Path

import pyarrow.parquet as pq
from scipy.stats import chi2, norm

REPO = Path(__file__).resolve().parents[2]
C95 = chi2.ppf(0.95, 1)


def logit(v):
    return math.log(v) - math.log1p(-v)


def summary():
    with open(REPO / "tests/golden/recovery/summary.csv", newline="") as f:
        return [r for r in csv.DictReader(line for line in f if not line.startswith("#"))]


def group(r):
    sid = int(r["scenario"])
    return "A" if r["pd_verdict"] == "DEFERRED" else "BCD"[(sid // 3) % 3]


def bartlett(rows):
    print("## Implied Bartlett factor k = 3.8415 / F^-1(coverage), groups B-D (pinned profile coverage, R = 1,000)")
    print("| group | T | PD: coverage, k (median, range) | rho: coverage, k (median, range) |")
    print("|---|---|---|---|")
    for g, T in (("B", 20), ("C", 40), ("D", 100)):
        sel = [r for r in rows if group(r) == g]
        cells = []
        for p in ("pd", "rho"):
            cov = [float(r[f"{p}_profile_coverage"]) for r in sel]
            ks = sorted(C95 / chi2.ppf(c, 1) for c in cov)
            pooled = sum(cov) / len(cov)
            cells.append(f"{pooled:.4f}, k {C95 / chi2.ppf(pooled, 1):.3f} ({ks[0]:.3f}-{ks[-1]:.3f})")
        print(f"| {g} | {T} | {cells[0]} | {cells[1]} |")
    print()


def shift(rows):
    print("## Percentile interval under a bias of b (logit sd units): coverage ~ Phi(1.96 + 2b) - Phi(-1.96 + 2b)")
    print("| scenario | T | rho bias (sd) | predicted percentile coverage | pinned profile | pinned iid percentile |")
    print("|---|---|---|---|---|---|")
    for r in rows:
        if group(r) != "B":
            continue
        sid = int(r["scenario"])
        true_rho = [0.02, 0.12, 0.24][(sid // 9) % 3]
        b = (logit(float.fromhex(r["rho_mean_hex"])) - logit(true_rho)) / float.fromhex(r["rho_sd_u_hex"])
        pred = norm.cdf(1.96 + 2 * b) - norm.cdf(-1.96 + 2 * b)
        print(f"| {sid} | 20 | {b:+.2f} | {pred:.3f} | {r['rho_profile_coverage']} | {r['rho_boot_coverage']} |")
    print()


def gap(rows):
    grp = {int(r["scenario"]): group(r) for r in rows}
    acc = {}
    for r in pq.read_table(REPO / "studies/bayes-coverage/fits.parquet").to_pylist():
        g = grp[r["scenario"]]
        for p in ("pd", "rho"):
            t = r[f"{p}_true"]
            prof = not r[f"prof_flags_{p}"] & 4 and r[f"prof_lo_{p}"] <= t <= r[f"prof_hi_{p}"]
            jef = r[f"jeffreys_et_lo_{p}"] <= t <= r[f"jeffreys_et_hi_{p}"]
            a = acc.setdefault((g, p), [0, 0, 0])
            a[0] += prof
            a[1] += jef
            a[2] += 1
    print("## Profile and Jeffreys equal-tailed coverage, pooled per group (S-9's rows, the same panels)")
    print("| group | PD: profile, Jeffreys | rho: profile, Jeffreys |")
    print("|---|---|---|")
    for g in "ABCD":
        cells = [f"{acc[(g, p)][0] / acc[(g, p)][2]:.4f}, {acc[(g, p)][1] / acc[(g, p)][2]:.4f}" for p in ("pd", "rho")]
        print(f"| {g} | {cells[0]} | {cells[1]} |")
    print()


if __name__ == "__main__":
    rows = summary()
    bartlett(rows)
    shift(rows)
    gap(rows)
