# SPDX-License-Identifier: Apache-2.0
"""Replicate the AR(1)-factor binomial log-likelihood with plain scipy (M6; docs/methodology/ar1_binomial_mle.md).

Model: given Z_t, d_t ~ Binomial(n_t, p(Z_t)) with p(z) = Phi((Phi^-1(PD) - sqrt(rho) z) / sqrt(1 - rho)); the factor is
AR(1) with unit variance, Z_1 ~ N(0, 1), Z_t | Z_{t-1} ~ N(phi Z_{t-1}, 1 - phi^2). The log-likelihood is
    l = log integral of prod_t C(n_t, d_t) p(z_t)^d_t (1 - p(z_t))^(n_t - d_t) dP(z_1, ..., z_T).

For every row of tests/golden/ar1/likelihood.csv (engine and ref values at fixed parameters on hand-made panels,
written by tests/studies/ar1_goldens.cpp), this script computes l on its own, in ways that share nothing with the
engine's node placement or with ref's grid. Two formulations of the factor's law are used, so that an error derived the
same way in the engine and in one of them (the transition variance 1 - phi^2, say) shows up against the other:
  - the Markov formulation, through the transition density N(phi z_{t-1}, 1 - phi^2), which the engine also uses:
    1. a forward filter on one uniform grid of composite-Simpson nodes over [-12, 12], with the transition matrix built
       from scipy.stats.norm.logpdf and the recursion in log space (scipy.special.logsumexp); every row;
    2. for T = 2 and 3 and n <= 1000, direct integration by nested scipy.integrate.quad: given Z_2, the other periods'
       factors are independent; scaled by exp(-l) from method 1 so the integral is of order one;
  - the joint formulation, which does not go through transitions at all:
    3. for T = 2 and 3 (every n), adaptive cubature (scipy.integrate.cubature, Gauss-Kronrod 21) of the explicit
       joint density: (Z_1, ..., Z_T) multivariate normal with covariance phi^|i - j| (a Toeplitz matrix,
       scipy.stats.multivariate_normal), times the binomial terms; over the joint mode +- 20 Laplace scales, both
       found from the joint density itself; scaled the same way.
It fails if any differs from the engine's value by more than its registered tolerance (tests/tolerances.toml),
relative to max(1, |l|): TOL_SCIPY_AR1_LOGLIK_REL for methods 1 and 2, TOL_SCIPY_AR1_JOINT_REL for method 3.

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/ar1_binomial_mle.py
"""

import csv
import math
import sys
from pathlib import Path

import numpy as np
import tomllib
from scipy import integrate, optimize, special, stats
from scipy.linalg import toeplitz

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = ROOT / "tests" / "golden" / "ar1" / "likelihood.csv"
Z_MAX = 12.0  # the hardest golden row (PD 0.001, rho 0.05, 90 of 1000) has its mass near z = -8


def log_binomial_term(n, d, pd, rho, z):
    """log C(n, d) + d log p(z) + (n - d) log(1 - p(z)), with log Phi from scipy.special.log_ndtr.

    log C(n, d) is the log of the exact integer math.comb(n, d): a difference of gammaln values loses about 2e-11
    to cancellation at n = 10^4, which would dominate the comparison.
    """
    x = (stats.norm.ppf(pd) - math.sqrt(rho) * z) / math.sqrt(1.0 - rho)
    lc = math.log(math.comb(n, d))
    out = np.full_like(z, lc, dtype=float)
    if d > 0:
        out = out + d * special.log_ndtr(x)
    if n - d > 0:
        out = out + (n - d) * special.log_ndtr(-x)
    return out


def grid_filter(ns, ds, pd, rho, phi):
    """Method 1 (Markov formulation): forward filter on a uniform Simpson grid over [-Z_MAX, Z_MAX].

    The step resolves the sharpest period.
    """
    step = min(0.02, 0.5 / math.sqrt(max(ns)))
    m = math.ceil(2.0 * Z_MAX / step)
    m += m % 2  # Simpson needs an even number of intervals
    z = np.linspace(-Z_MAX, Z_MAX, m + 1)
    h = z[1] - z[0]
    w = np.full(m + 1, 2.0 * h / 3.0)  # composite Simpson: h/3 (1, 4, 2, 4, ..., 2, 4, 1)
    w[1::2] = 4.0 * h / 3.0
    w[0] = w[-1] = h / 3.0
    log_w = np.log(w)
    sd = math.sqrt(1.0 - phi * phi)
    total = 0.0
    log_mass = None  # log filtered probability mass at each node, normalised
    for t, (n, d) in enumerate(zip(ns, ds)):
        lf = log_binomial_term(n, d, pd, rho, z)
        if t == 0:
            log_pred = stats.norm.logpdf(z)
        else:
            # log_pred[i] = log sum_j mass_j K(z_j -> z_i), in chunks of target nodes to bound memory.
            log_pred = np.empty_like(z)
            for lo in range(0, z.size, 512):
                zi = z[lo : lo + 512, None]
                log_k = stats.norm.logpdf(zi, loc=phi * z[None, :], scale=sd)
                log_pred[lo : lo + 512] = special.logsumexp(log_k + log_mass[None, :], axis=1)
        a = log_w + lf + log_pred
        lse = special.logsumexp(a)
        total += lse
        log_mass = a - lse
    return total


def quad_nested(ns, ds, pd, rho, phi, scale):
    """Method 2 (Markov formulation): nested 1-D integration for T = 2 or 3, scaled by exp(-scale).

    Given Z_2, the factor values Z_1 and Z_3 are independent, each N(phi Z_2, 1 - phi^2) (the chain is reversible), so
        L = integral f_2(z_2) phi(z_2) prod_{s != 2} [integral f_s(z_s) N(z_s; phi z_2, 1 - phi^2) dz_s] dz_2,
    a 1-D integral of products of 1-D integrals, each done by scipy.integrate.quad.
    """
    sd = math.sqrt(1.0 - phi * phi)
    opts = {"epsabs": 0.0, "epsrel": 1e-12, "limit": 200}

    def log_f(t, z):
        return float(log_binomial_term(ns[t], ds[t], pd, rho, np.array([z]))[0])

    grid = np.linspace(-Z_MAX, Z_MAX, 3601)
    # each period's own likelihood mode; the outer integrand's mass lies near one of them or near 0
    modes = sorted({0.0, *(float(grid[np.argmax(log_binomial_term(n, d, pd, rho, grid))]) for n, d in zip(ns, ds))})

    def neighbour(t, z2):
        # integral f_t(z) N(z; phi z2, 1 - phi^2) dz, scaled by exp(-c) with c the integrand's largest log on a grid
        g = log_binomial_term(ns[t], ds[t], pd, rho, grid) + stats.norm.logpdf(grid, loc=phi * z2, scale=sd)
        k = int(np.argmax(g))
        c = float(g[k])
        v, _ = integrate.quad(lambda z: math.exp(log_f(t, z) + stats.norm.logpdf(z, loc=phi * z2, scale=sd) - c),
                              -Z_MAX, Z_MAX, points=[float(grid[k])], **opts)
        return c + math.log(v)

    def integrand(z2):
        v = log_f(1, z2) + stats.norm.logpdf(z2)
        for t in range(len(ns)):
            if t != 1:
                v += neighbour(t, z2)
        return math.exp(v - scale)

    val, _ = integrate.quad(integrand, -Z_MAX, Z_MAX, points=modes, **opts)
    return scale + math.log(val)


def cubature_joint(ns, ds, pd, rho, phi, scale):
    """Method 3 (joint formulation): adaptive cubature of the explicit joint density (Toeplitz covariance phi^|i-j|),
    scaled by exp(-scale).

    The joint log-density is concave (a normal times probit binomial terms), so it has one mode; it is found by BFGS on
    the joint log-density itself, its Hessian there by central differences gives Laplace scales s_i, and the cubature
    runs over mode +- 20 s_i with the mode as a breakpoint. Over the whole of [-Z_MAX, Z_MAX]^T the rule's first
    evaluation can miss a narrow ridge (n = 10^4, phi = 0.95) and report convergence on nothing.
    """
    T = len(ns)
    mvn = stats.multivariate_normal(mean=np.zeros(T), cov=toeplitz(phi ** np.arange(T)))

    def log_joint(x):  # x has shape (points, T)
        v = mvn.logpdf(x)
        for t in range(T):
            v = v + log_binomial_term(ns[t], ds[t], pd, rho, x[:, t])
        return v

    def neg(x):
        return -float(log_joint(np.atleast_2d(x))[0])

    mode = optimize.minimize(neg, np.zeros(T), method="BFGS", options={"gtol": 1e-9}).x
    h = 1e-4
    eye = np.eye(T) * h
    hess = np.array([[(neg(mode + eye[i] + eye[j]) - neg(mode + eye[i] - eye[j]) - neg(mode - eye[i] + eye[j])
                       + neg(mode - eye[i] - eye[j])) / (4 * h * h) for j in range(T)] for i in range(T)])
    s = np.sqrt(np.diag(np.linalg.inv(hess)))
    lo, hi = mode - 20.0 * s, mode + 20.0 * s
    res = integrate.cubature(lambda x: np.exp(log_joint(x) - scale), lo, hi, rule="gk21", rtol=1e-11, atol=0.0,
                             max_subdivisions=100000, points=[mode])
    if res.status != "converged" or not res.estimate > 0:
        return math.nan
    return scale + math.log(res.estimate)


def main():
    with open(GOLDEN) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    with open(ROOT / "tests" / "tolerances.toml", "rb") as f:
        reg = {k: v["value"] for k, v in tomllib.load(f).items()}
    failures = 0
    worst = {"grid filter": 0.0, "nested quad": 0.0, "joint cubature": 0.0}
    count = dict.fromkeys(worst, 0)
    tol = {"grid filter": reg["TOL_SCIPY_AR1_LOGLIK_REL"], "nested quad": reg["TOL_SCIPY_AR1_LOGLIK_REL"],
           "joint cubature": reg["TOL_SCIPY_AR1_JOINT_REL"]}
    for r in rows:
        ns = [int(v) for v in r["n"].split()]
        ds = [int(v) for v in r["d"].split()]
        pd, rho, phi = (float.fromhex(r[k]) for k in ("pd_hex", "rho_hex", "phi_hex"))
        engine = float.fromhex(r["loglik_hex"])
        scale = max(1.0, abs(engine))
        g = grid_filter(ns, ds, pd, rho, phi)
        values = {"grid filter": g}
        if len(ns) <= 3 and max(ns) <= 1000:
            values["nested quad"] = quad_nested(ns, ds, pd, rho, phi, g)
        if len(ns) <= 3:
            values["joint cubature"] = cubature_joint(ns, ds, pd, rho, phi, g)
        for method, v in values.items():
            e = abs(v - engine) / scale
            count[method] += 1
            worst[method] = max(worst[method], e) if not math.isnan(e) else math.inf
            if not e <= tol[method]:
                failures += 1
                print(f"FAIL {method}: panel {r['panel']} phi {phi}: scipy {v!r} engine {engine!r} ({e:.2e})")
    print(f"{len(rows)} rows; worst relative difference from the engine:")
    for method, w in worst.items():
        print(f"  {method}: {count[method]} rows, {w:.2e} (tolerance {tol[method]:g})")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
