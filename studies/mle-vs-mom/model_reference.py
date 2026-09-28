# SPDX-License-Identifier: Apache-2.0
"""The estimator-comparison pass (S-8, S-27, S-28): reference numbers from the model alone.

No recovery panel is read, and no estimator of the pass is run on simulated panels. Three parts:

  s8   MoM against the MLE in the large-n limit. As n -> inf the period rates p_t = Phi(x_t),
       x_t ~ N(m, s^2), m = probit(PD)/sqrt(1 - rho), s^2 = rho/(1 - rho), are observed exactly. The MLE is
       the mean and 1/T variance of x_t; MoM (joint-default-probability form) is PD-hat = mean p_t and
       rho-hat solving Phi2(c, c; rho) = mean p_t^2, c = probit(PD-hat), the large-n limit of the pooled
       d(d - 1)/(n(n - 1)); rho-hat is set to the box floor 1e-3 when mean p_t^2 <= PD-hat^2.
       Prints RMSE(MoM)/RMSE(MLE) for rho and PD, and both estimators' bias of rho.
  s27  The Vasicek-rate MLE reads probit(d/n) as the latent x, so binomial noise adds
       B = E_z[p(1 - p) / (n phi(probit p)^2)] to the variance of x (delta method; unreliable where n p(z)
       is small in much of the factor's range, which is where zero-default periods are common). In SE units
       rho-hat_rate's bias is about B / (s^2 sqrt(2/T)). n B does not depend on n, so the bias falls as 1/n.
  s28  The probability that a period has no defaults, P0 = E_z[(1 - p(z))^n], and the share of T-period
       panels with at least one, 1 - (1 - P0)^T: exactly the share the parity rate estimator refuses.

    uv run --no-project --with scipy python studies/mle-vs-mom/model_reference.py [s8|s27|s28]
"""

import math
import random
import sys

import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy import integrate
from scipy.optimize import brentq
from scipy.stats import norm

PDS = [0.001, 0.01, 0.05]
RHOS = [0.02, 0.12, 0.24]
TS = [20, 40, 100]
NODES, WEIGHTS = hermegauss(96)
WEIGHTS = WEIGHTS / WEIGHTS.sum()


def p_of_z(pd, rho, z):
    return norm.cdf((norm.ppf(pd) - math.sqrt(rho) * z) / math.sqrt(1 - rho))


def phi2_equal(c, r):
    """Phi2(c, c; r) = E_z[Phi((c - sqrt(r) z)/sqrt(1 - r))^2]."""
    return float(np.sum(WEIGHTS * norm.cdf((c - math.sqrt(r) * NODES) / math.sqrt(1 - r)) ** 2))


def mom_large_n(ps):
    pd = sum(ps) / len(ps)
    pd2 = sum(p * p for p in ps) / len(ps)
    c = norm.ppf(pd)
    if pd2 <= pd * pd:
        return pd, 1e-3
    if phi2_equal(c, 0.999) < pd2:
        return pd, 0.5
    return pd, min(0.5, max(1e-3, brentq(lambda r: phi2_equal(c, r) - pd2, 1e-9, 0.999, xtol=1e-12)))


def mle_large_n(xs):
    mu = sum(xs) / len(xs)
    v = sum((x - mu) ** 2 for x in xs) / len(xs)
    return norm.cdf(mu / math.sqrt(1 + v)), v / (1 + v)


def s8(draws=2000):
    print("| PD | ρ | T | RMSE ratio MoM/MLE, ρ | RMSE ratio, PD | bias of ρ̂: MoM | bias of ρ̂: MLE |")
    print("|---|---|---|---|---|---|---|")
    for pd in PDS:
        for rho in RHOS:
            for T in TS:
                rng = random.Random(3)
                m, s = norm.ppf(pd) / math.sqrt(1 - rho), math.sqrt(rho / (1 - rho))
                e = {k: [] for k in ("mom_pd", "mom_rho", "mle_pd", "mle_rho")}
                for _ in range(draws):
                    xs = [rng.gauss(m, s) for _ in range(T)]
                    a, b = mom_large_n([norm.cdf(x) for x in xs]), mle_large_n(xs)
                    e["mom_pd"].append(a[0] - pd)
                    e["mom_rho"].append(a[1] - rho)
                    e["mle_pd"].append(b[0] - pd)
                    e["mle_rho"].append(b[1] - rho)
                rmse = {k: math.sqrt(sum(v * v for v in x) / len(x)) for k, x in e.items()}
                bias = {k: sum(x) / len(x) for k, x in e.items()}
                print(f"| {pd} | {rho} | {T} | {rmse['mom_rho'] / rmse['mle_rho']:.2f} | "
                      f"{rmse['mom_pd'] / rmse['mle_pd']:.2f} | {bias['mom_rho']:+.4f} | {bias['mle_rho']:+.4f} |")


def binomial_probit_variance(pd, rho, n):
    def f(z):
        p = p_of_z(pd, rho, z)
        return p * (1 - p) / (n * norm.pdf(norm.ppf(p)) ** 2) * norm.pdf(z)
    return integrate.quad(f, -8, 8, limit=400)[0]


def s27():
    ns = [100, 1000, 10000, 100000, 1000000]
    print("| PD | ρ | " + " | ".join(f"n = {n:,}" for n in ns) + " |")
    print("|---|---|" + "---|" * len(ns))
    for pd in PDS:
        for rho in RHOS:
            s2 = rho / (1 - rho)
            cells = []
            for n in ns:
                b = binomial_probit_variance(pd, rho, n)
                cells.append(" / ".join(f"{b / (s2 * math.sqrt(2 / T)):.2f}" for T in TS))
            print(f"| {pd} | {rho} | " + " | ".join(cells) + " |")
    print("(each cell: rho-hat_rate's bias in SE units at T = 20 / 40 / 100)")


def p_zero(pd, rho, n):
    def f(z):
        return math.exp(n * math.log1p(-p_of_z(pd, rho, z))) * norm.pdf(z)
    return integrate.quad(f, -12, 12, limit=400, points=[-3, 0, 3])[0]


def s28():
    print("| PD | ρ | n | P(d = 0) | refused, T = 20 | T = 40 | T = 100 |")
    print("|---|---|---|---|---|---|---|")
    for pd in PDS:
        for rho in RHOS:
            for n in (100, 1000, 10000):
                p0 = p_zero(pd, rho, n)
                refused = " | ".join(f"{1 - (1 - p0) ** T:.4f}" for T in TS)
                print(f"| {pd} | {rho} | {n:,} | {p0:.4f} | {refused} |")


if __name__ == "__main__":
    parts = sys.argv[1:] or ["s8", "s27", "s28"]
    for part in parts:
        {"s8": s8, "s27": s27, "s28": s28}[part]()
        print()
