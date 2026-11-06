# SPDX-License-Identifier: Apache-2.0
"""Replicate the engine's Pluto-Tasche most-prudent upper bounds with plain scipy (S-6, D-214, D-215).

The method (docs/methodology/pluto_tasche.md). A panel of T periods with n obligors each and D defaults in total.
At confidence gamma the bound is the PD p with P(S <= D | p) = 1 - gamma:
    independent:  S ~ Bin(N, p), N = n T: the one-sided Clopper-Pearson bound, here from the beta quantile
                  p = B^-1(gamma; D + 1, N - D);
    one factor:   S = sum_t S_t, the S_t independent with the Vasicek mixed binomial PMF at (p, rho_a):
                  P(S_t = d) = integral C(n, d) p(z)^d (1 - p(z))^(n - d) phi(z) dz,
                  p(z) = Phi((Phi^-1(p) - sqrt(rho_a) z) / sqrt(1 - rho_a)); the T-fold convolution, truncated
                  at D, solved for p by brentq in logit p.
The one-period PMF here is a trapezoid sum over z in [-12, 12] at step 1/200, in log space (log_ndtr, gammaln): no
code is shared with core's quadrature. The integrand is analytic and Gaussian-tailed, so the rule converges
geometrically; studies/pluto-tasche/reference.py checks the same sum against scipy.integrate.quad at 1.6e-13.

It checks every row of tests/golden/pluto_tasche/bounds.csv (written by tests/studies/pluto_tasche_goldens.cpp) in
logit p, to TOL_SCIPY_PT_LOGIT from tests/tolerances.toml.

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/pluto_tasche.py
"""

import csv
import math
import sys
from pathlib import Path

import numpy as np
import tomllib
from scipy.optimize import brentq
from scipy.special import gammaln, log_ndtr, logsumexp, ndtri
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = ROOT / "tests" / "golden" / "pluto_tasche" / "bounds.csv"

Z = np.linspace(-12.0, 12.0, 4801)
LOG_W = np.log(np.full(Z.size, Z[1] - Z[0])) - 0.5 * Z * Z - 0.5 * math.log(2.0 * math.pi)
LOG_W[0] -= math.log(2.0)
LOG_W[-1] -= math.log(2.0)
XTOL = 1e-13  # brentq's tolerance in logit p


def tolerance(name):
    with open(ROOT / "tests" / "tolerances.toml", "rb") as f:
        return tomllib.load(f)[name]["value"]


def expit(u):
    return 1.0 / (1.0 + math.exp(-u))


def logit(p):
    return math.log(p) - math.log1p(-p)


def period_pmf(n, pd, rho, dmax):
    """P(S_t = d), d = 0..dmax (zero above n)."""
    x = (ndtri(pd) - math.sqrt(rho) * Z) / math.sqrt(1.0 - rho)
    lp, lq = log_ndtr(x), log_ndtr(-x)
    d = np.arange(min(dmax, n) + 1, dtype=float)
    lc = gammaln(n + 1.0) - gammaln(d + 1.0) - gammaln(n - d + 1.0)
    terms = lc[None, :] + d[None, :] * lp[:, None] + (n - d)[None, :] * lq[:, None] + LOG_W[:, None]
    out = np.zeros(dmax + 1)
    out[: d.size] = np.exp(logsumexp(terms, axis=0))
    return out


def total_cdf(n, t, pd, rho, d_total):
    """P(S <= D): the T-fold convolution of the one-period PMF, truncated at D."""
    f = period_pmf(n, pd, rho, d_total)
    m = f.size
    out = np.zeros(m)
    out[0] = 1.0
    while t:
        if t & 1:
            out = np.convolve(out, f)[:m]
        t >>= 1
        if t:
            f = np.convolve(f, f)[:m]
    return out.sum()


def bound_factor(n, t, d_total, rho_a, gamma):
    def g(u):
        return math.log(max(total_cdf(n, t, expit(u), rho_a, d_total), 1e-300)) - math.log1p(-gamma)

    return brentq(g, -30.0, 15.0, xtol=XTOL, rtol=4 * np.finfo(float).eps, maxiter=500)


def main():
    tol = tolerance("TOL_SCIPY_PT_LOGIT")
    with open(GOLDEN, newline="") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    worst = {"independent": 0.0, "factor": 0.0}
    bad = 0
    for r in rows:
        n, t, d = int(r["n"]), int(r["periods"]), int(r["d_total"])
        rho, gamma = float.fromhex(r["rho_a_hex"]), float.fromhex(r["gamma_hex"])
        want = float.fromhex(r["logit_bound_hex"])
        if r["method"] == "independent":
            got = logit(float(beta.ppf(gamma, d + 1, n * t - d)))
        else:
            got = bound_factor(n, t, d, rho, gamma)
        diff = abs(got - want)
        worst[r["method"]] = max(worst[r["method"]], diff)
        if not diff <= tol:
            bad += 1
            print(f"MISMATCH {r['method']} n={n} T={t} D={d} rho_a={rho} gamma={gamma}: "
                  f"engine {want!r} scipy {got!r} diff {diff:.3g}")
    print(f"{len(rows)} rows; worst |logit engine - logit scipy|: independent {worst['independent']:.3g}, "
          f"factor {worst['factor']:.3g} (tolerance {tol:g})")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
