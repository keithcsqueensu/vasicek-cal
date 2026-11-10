# SPDX-License-Identifier: Apache-2.0
"""Replicate vasicek-cal's multi-grade binomial-mixture likelihood with plain scipy (M6, D-151, D-411).

Model (docs/methodology/multigrade_mle.md). A period has K grades; grade k has n_k obligors and d_k defaults.
Every obligor loads on the same factor Z ~ N(0, 1) with the same rho, and given Z grade k defaults with
probability
    p_k(Z) = Phi((Phi^-1(PD_k) - sqrt(rho) Z) / sqrt(1 - rho)),
independently across grades, so the period log-likelihood is
    l(PD_1, ..., PD_K, rho) = log E_Z[ prod_k binom.pmf(d_k; n_k, p_k(Z)) ].

This script is the likelihood half of the replication (PR 1; the estimator half arrives with the estimator in
PR 2). With ordinary tools only (binom.logpmf in log space, scipy.optimize.minimize_scalar for the integrand's mode,
scipy.integrate.quad on the window where the integrand is within e^-50 of its maximum), it computes l on every
row of tests/golden/multigrade/likelihood.csv and checks it against the committed value plus the binomial
coefficients: the values the engine is tested against (to within 6 eps, unit_multigrade). The tolerance is
TOL_SCIPY_MULTIGRADE_LL_REL in tests/tolerances.toml.

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/multigrade_mle.py
"""

import csv
import math
import sys
import tomllib
from pathlib import Path

from scipy import integrate, optimize, special, stats

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = ROOT / "tests" / "golden" / "multigrade" / "likelihood.csv"
DROP = 50.0  # integrate where the integrand is within e^-50 of its maximum


def period_loglik(pds, rho, grades):
    """log E_Z[prod_k binom.pmf(d_k; n_k, p_k(Z))], integrating around the mode of the integrand."""
    a, b = math.sqrt(rho), math.sqrt(1.0 - rho)
    cs = [stats.norm.ppf(p) for p in pds]
    live = [(c, n, d) for c, (n, d) in zip(cs, grades) if n > 0]

    def h(z):  # log of prod_k binom.pmf(d_k; n_k, p_k(z)) * phi(z), with log Phi from norm.logcdf in both tails
        v = stats.norm.logpdf(z)
        for c, n, d in live:
            x = (c - a * z) / b
            v += special.gammaln(n + 1) - special.gammaln(d + 1) - special.gammaln(n - d + 1)
            if d:
                v += d * stats.norm.logcdf(x)
            if n - d:
                v += (n - d) * stats.norm.logcdf(-x)
        return v

    # Start at the n-weighted mean of the grades' own large-n modes, then let Brent's method find the mode.
    z0 = sum(n * (c - b * stats.norm.ppf((d + 0.5) / (n + 1))) / a for c, n, d in live) / sum(n for _, n, _ in live)
    mode = optimize.minimize_scalar(lambda z: -h(z), bracket=(z0 - 1.0, z0 + 1.0)).x
    top = h(mode)

    def edge(direction):  # walk outwards, doubling the step, until h has dropped by DROP
        z, step = mode, 0.25 * direction
        while h(z) > top - DROP:
            z, step = z + step, 2.0 * step
        return z

    lo, hi = edge(-1.0), edge(+1.0)
    area, _ = integrate.quad(lambda z: math.exp(h(z) - top), lo, hi, points=[mode], epsabs=0.0, epsrel=1e-12,
                             limit=1000)
    return top + math.log(area)


def lchoose(n, d):
    return special.gammaln(n + 1) - special.gammaln(d + 1) - special.gammaln(n - d + 1)


def nan_inf(x):
    """x, or inf for a NaN: a running max() would otherwise drop a NaN and a check would pass (D-218)."""
    return math.inf if math.isnan(x) else x


def main():
    with open(ROOT / "tests" / "tolerances.toml", "rb") as f:
        tol = {k: v["value"] for k, v in tomllib.load(f).items()}
    with open(GOLDEN, newline="") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    worst, where = 0.0, ""
    for r in rows:
        rho = float.fromhex(r["rho_hex"])
        grades = [tuple(int(v) for v in g.split(":")) for g in r["grades"].split(";")]
        pds = [float.fromhex(p) for p in r["pd_hex"].split(";")]
        want = float.fromhex(r["log_integral_hex"]) + sum(lchoose(n, d) for n, d in grades if n > 0)
        got = period_loglik(pds, rho, grades)
        rel = nan_inf(abs(got - want) / max(1.0, abs(want)))
        if not rel <= worst:
            worst, where = rel, f"{r['kind']}, rho {rho:.4g}, K = {len(grades)}"
    passed = worst <= tol["TOL_SCIPY_MULTIGRADE_LL_REL"]
    print(f"{len(rows)} periods: worst relative difference {worst:.3g} ({where}) {'ok' if passed else 'FAIL'}")
    print("all checks passed" if passed else "SOME CHECKS FAILED")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
