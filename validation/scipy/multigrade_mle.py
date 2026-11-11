# SPDX-License-Identifier: Apache-2.0
"""Replicate vasicek-cal's multi-grade binomial-mixture likelihood with plain scipy (M6, D-151, D-411).

Model (docs/methodology/multigrade_mle.md). A period has K grades; grade k has n_k obligors and d_k defaults.
Every obligor loads on the same factor Z ~ N(0, 1) with the same rho, and given Z grade k defaults with
probability
    p_k(Z) = Phi((Phi^-1(PD_k) - sqrt(rho) Z) / sqrt(1 - rho)),
independently across grades, so the period log-likelihood is
    l(PD_1, ..., PD_K, rho) = log E_Z[ prod_k binom.pmf(d_k; n_k, p_k(Z)) ].

Two halves. The likelihood half: with ordinary tools only (binom.logpmf in log space, scipy.optimize.minimize_scalar for the integrand's mode,
scipy.integrate.quad on the window where the integrand is within e^-50 of its maximum), it computes l on every
row of tests/golden/multigrade/likelihood.csv and checks it against the committed value plus the binomial
coefficients: the values the engine is tested against (to within 6 eps, unit_multigrade). The tolerance is
TOL_SCIPY_MULTIGRADE_LL_REL in tests/tolerances.toml.

The estimate half (D-413): on the two panels of tests/golden/multigrade/estimates.csv it maximises the panel
log-likelihood, the sum of the periods' l, with scipy.optimize.minimize (L-BFGS-B in u_k = logit PD_k and
v = logit rho, inside the engine's box PD in [1e-5, 0.5], rho in [1e-3, 0.5], from each grade's pooled rate
(sum d + 1/2) / (sum n + 1) and rho = 0.1, with central-difference gradients). It checks that its estimate is
within TOL_SCIPY_MULTIGRADE_FIT_LOGIT of the engine's in every coordinate, that its log-likelihood at the engine's
estimate is no lower than at its own optimum minus TOL_SCIPY_MULTIGRADE_FIT_LL_ABS, and that it agrees with the
engine's committed log-likelihood there within TOL_SCIPY_MULTIGRADE_LL_REL.

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
ESTIMATES = ROOT / "tests" / "golden" / "multigrade" / "estimates.csv"
BOX = ((1e-5, 0.5), (1e-3, 0.5))  # the engine's default box: PD, rho
GRAD_STEP = 1e-5  # central differences in logit units: quad's 1e-12 noise over 1e-5 is a 1e-7 relative gradient error
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
    passed = estimates(tol) and passed
    print("all checks passed" if passed else "SOME CHECKS FAILED")
    return 0 if passed else 1


def logit(p):
    return math.log(p) - math.log1p(-p)


def expit(u):
    return special.expit(u)


def estimates(tol):
    with open(ESTIMATES, newline="") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    ok = True
    for r in rows:
        K = int(r["K"])
        panel = [[tuple(int(v) for v in g.split(":")) for g in period.split(";")] for period in r["counts"].split("|")]
        core = [logit(float.fromhex(p)) for p in r["pd_hex"].split(";")] + [logit(float.fromhex(r["rho_hex"]))]

        def loglik(x):
            pds, rho = [float(expit(u)) for u in x[:K]], float(expit(x[K]))
            return sum(period_loglik(pds, rho, grades) for grades in panel)

        def neg_and_grad(x):
            g = []
            for i in range(K + 1):
                up, dn = list(x), list(x)
                up[i] += GRAD_STEP
                dn[i] -= GRAD_STEP
                g.append(-(loglik(up) - loglik(dn)) / (2.0 * GRAD_STEP))
            return -loglik(x), g

        start = []
        for k in range(K):
            n = sum(period[k][0] for period in panel)
            d = sum(period[k][1] for period in panel)
            start.append(min(logit(BOX[0][1]), max(logit(BOX[0][0]), logit((d + 0.5) / (n + 1.0)))))
        start.append(logit(0.1))
        bounds = [tuple(logit(b) for b in BOX[0])] * K + [tuple(logit(b) for b in BOX[1])]
        res = optimize.minimize(neg_and_grad, start, jac=True, method="L-BFGS-B", bounds=bounds,
                                options={"ftol": 1e-15, "gtol": 1e-9, "maxiter": 500})
        gap = max(abs(a - b) for a, b in zip(res.x, core))
        at_core = loglik(core)
        beaten = -res.fun - at_core
        want = float.fromhex(r["loglik_hex"])
        rel = nan_inf(abs(at_core - want) / max(1.0, abs(want)))
        good = (gap <= tol["TOL_SCIPY_MULTIGRADE_FIT_LOGIT"] and beaten <= tol["TOL_SCIPY_MULTIGRADE_FIT_LL_ABS"]
                and rel <= tol["TOL_SCIPY_MULTIGRADE_LL_REL"])
        ok = ok and good
        print(f"estimate, {r['name']} (K = {K}, T = {len(panel)}): worst logit gap {gap:.3g}, scipy optimum above the "
              f"engine's by {beaten:.3g}, log-likelihood at the engine's estimate {rel:.3g} relative "
              f"({res.nit} iterations) {'ok' if good else 'FAIL'}")
    return ok


if __name__ == "__main__":
    sys.exit(main())
