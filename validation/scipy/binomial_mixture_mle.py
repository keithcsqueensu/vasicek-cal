# SPDX-License-Identifier: Apache-2.0
"""Replicate vasicek-cal's binomial-mixture MLE with plain scipy (M1.9, D-126).

Model (docs/methodology/binomial_mixture_mle.md). Period t has n_t obligors and d_t defaults.
Given the systematic factor Z ~ N(0, 1), each obligor defaults with probability
    p(Z) = Phi((Phi^-1(PD) - sqrt(rho) Z) / sqrt(1 - rho)),
so d_t is Binomial(n_t, p(Z)) given Z, and the panel log-likelihood is
    ll(PD, rho) = sum_t log E_Z[ binom.pmf(d_t; n_t, p(Z)) ].

This script uses ordinary tools only (binom.logpmf in log space, scipy.integrate.quad around the
integrand's mode, scipy.optimize.minimize) and checks three things against the engine's
committed results on the six M1.7 panels:
  1. the panel log-likelihood at each of the engine's surface points;
  2. the log-likelihood at the engine's estimate against that at scipy's own optimum over the
     same box (the engine's estimate should be as good, to optimiser accuracy);
  3. the estimates, in units of the engine's standard errors.
Tolerances come from tests/tolerances.toml (the TOL_SCIPY_* entries).

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/binomial_mixture_mle.py
"""

import csv
import math
import sys
import tomllib
from collections import Counter
from pathlib import Path

from scipy import integrate, optimize, stats

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = ROOT / "tests" / "golden"
PD_BOX, RHO_BOX = (1e-4, 0.2), (1e-3, 0.5)  # the engine's estimation box (D-115)
DROP = 50.0  # integrate where the integrand is within e^-50 of its maximum


def period_loglik(pd, rho, n, d):
    """log E_Z[binom.pmf(d; n, p(Z))], integrating around the mode of the integrand."""
    c, a, b = stats.norm.ppf(pd), math.sqrt(rho), math.sqrt(1.0 - rho)

    def h(z):  # log of binom.pmf(d; n, p(z)) * phi(z)
        return stats.binom.logpmf(d, n, stats.norm.cdf((c - a * z) / b)) + stats.norm.logpdf(z)

    # Start where p(z) = (d + 1/2) / (n + 1), then let Brent's method find the mode.
    z0 = (c - b * stats.norm.ppf((d + 0.5) / (n + 1))) / a
    mode = optimize.minimize_scalar(lambda z: -h(z), bracket=(z0 - 1.0, z0 + 1.0)).x
    top = h(mode)

    def edge(direction):  # walk outwards, doubling the step, until h has dropped by DROP
        z, step = mode, 0.25 * direction
        while h(z) > top - DROP:
            z, step = z + step, 2.0 * step
        return z

    lo, hi = edge(-1.0), edge(+1.0)
    area, _ = integrate.quad(lambda z: math.exp(h(z) - top), lo, hi, points=[mode], epsabs=0.0,
                             epsrel=1e-12, limit=500)
    return top + math.log(area)


def panel_loglik(pd, rho, panel):
    """Sum over periods; periods with the same (n, d) are integrated once."""
    return sum(k * period_loglik(pd, rho, n, d) for (n, d), k in Counter(panel).items())


def logit(v):
    return math.log(v / (1.0 - v))


def inv_logit(u):
    return 1.0 / (1.0 + math.exp(-u))


def scipy_fit(panel):
    """Maximise the panel log-likelihood over the box, in logit coordinates like the engine."""
    start = (logit(min(max(sum(d for _, d in panel) / sum(n for n, _ in panel), 2e-4), 0.19)), logit(0.1))
    bounds = [(logit(PD_BOX[0]), logit(PD_BOX[1])), (logit(RHO_BOX[0]), logit(RHO_BOX[1]))]
    res = optimize.minimize(lambda u: -panel_loglik(inv_logit(u[0]), inv_logit(u[1]), panel), start,
                            method="Nelder-Mead", bounds=bounds, options={"xatol": 1e-7, "fatol": 1e-10,
                                                                          "maxiter": 2000})
    return inv_logit(res.x[0]), inv_logit(res.x[1]), -res.fun


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def main():
    with open(ROOT / "tests" / "tolerances.toml", "rb") as f:
        tol = {k: v["value"] for k, v in tomllib.load(f).items()}
    panels = {}
    for r in read_csv(GOLDEN / "validation" / "panels.csv"):
        panels.setdefault(int(r["scenario"]), []).append((int(r["n"]), int(r["d"])))
    ok = True

    # 1. The panel log-likelihood at every committed surface point of the engine.
    worst, where = 0.0, ""
    for r in read_csv(GOLDEN / "validation" / "surface.csv"):
        pd, rho, core = (float.fromhex(r[k]) for k in ("pd_hex", "rho_hex", "loglik_hex"))
        rel = abs(panel_loglik(pd, rho, panels[int(r["scenario"])]) / core - 1.0)
        if rel > worst:
            worst, where = rel, f"scenario {r['scenario']}, PD {pd:.4g}, rho {rho:.4g}"
    passed = worst <= tol["TOL_SCIPY_SURFACE_LL_REL"]
    ok &= passed
    print(f"1. surface: worst relative difference {worst:.3g} ({where}) {'ok' if passed else 'FAIL'}")

    # 2 and 3. scipy's own optimum against the engine's estimate.
    for r in read_csv(GOLDEN / "xref" / "core_estimates.csv"):
        s = int(r["scenario"])
        pd_c, rho_c, se_pd, se_rho = (float.fromhex(r[k]) for k in ("pd_hex", "rho_hex", "se_pd_hex", "se_rho_hex"))
        pd_s, rho_s, ll_s = scipy_fit(panels[s])
        gap = ll_s - panel_loglik(pd_c, rho_c, panels[s])  # >= 0 up to optimiser accuracy
        z_pd, z_rho = abs(pd_s - pd_c) / se_pd, abs(rho_s - rho_c) / se_rho
        passed = (abs(gap) <= tol["TOL_SCIPY_OPTIMUM_LL_ABS"] and max(z_pd, z_rho) <= tol["TOL_SCIPY_ESTIMATE_SE"])
        ok &= passed
        print(f"scenario {s}: engine PD {pd_c:.6g} rho {rho_c:.6g} | scipy PD {pd_s:.6g} rho {rho_s:.6g} | "
              f"2. ll(scipy optimum) - ll(engine estimate) {gap:.3g} | 3. difference {z_pd:.3g} / {z_rho:.3g} SE "
              f"{'ok' if passed else 'FAIL'}")
    print("all checks passed" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
