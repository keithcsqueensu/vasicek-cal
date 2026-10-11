# SPDX-License-Identifier: Apache-2.0
"""S-36 (rating philosophy and fitted correlation): reference values from the model only.

No panel is simulated or fitted here. Everything follows from the DGP's definition
(PREDICTION.md, "The DGP") by quadrature over the factor and the obligor population.

Model. Obligor i has a TTC probit c_i = Phi^-1(TTC PD_i). Its PIT probit in a period with factor z is
x_i(z) = (c_i - sqrt(rho) z) / sqrt(1 - rho) (engine convention: a higher z is better). It is rated on
the hybrid score s_i(z) = (1 - lam) c_i + lam x_i(z); lam = 0 is TTC, lam = 1 is PIT. Grade g holds
the obligors whose score lies in [b_{g-1}, b_g), the master scale's probit boundaries (open-ended at
both ends of the scale).

Narrow-grade limit (closed form). Inside a grade s is about constant, so c = (s + lam sqrt(rho) z /
sqrt(1 - rho)) / k with k = 1 - lam + lam / sqrt(1 - rho), and the grade's PIT probit is linear in z
with slope -sqrt(rho) (1 - lam) / (k sqrt(1 - rho)). The grade's probit default rate is then Gaussian
with variance v = rho (1 - lam)^2 / ((1 - rho) k^2), and its Vasicek-implied correlation is
rho_eff = v / (1 + v): rho at lam = 0, 0 at lam = 1.

Finite grades (the reference used for predictions). For each z on a Gauss-Hermite grid, the grade's
default rate is R_g(z) = E[Phi(x_i(z)) | s_i(z) in grade g], integrating c over the obligor
population, and its share of obligors is w_g(z). The n -> infinity pseudo-true values of the
Vasicek-rate MLE are the moments of X_g = Phi^-1(R_g(Z)) (a Gaussian MLE matches the first two
moments): v_g = Var X_g, rho_g = v_g / (1 + v_g), PD_g = Phi(E X_g / sqrt(1 + v_g)). The binomial MLE
at finite n targets about the same values (S-27), which PREDICTION.md's tolerances allow for.

    uv run --no-project --with scipy==1.16.2 python studies/pit-ttc-rho/reference.py [--out FILE]
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.special import log_ndtr, logsumexp, ndtr, ndtri, ndtri_exp, roots_hermitenorm

K = 17
PD_K = np.array([0.0003 * 1.5 ** (k - 1) for k in range(1, K + 1)])  # master scale (D-151)
SCALES = {"L": (5, 0.24), "M": (8, 0.12), "H": (11, 0.02)}            # centre, rho (D-151)
LAMBDAS = [0.0, 0.25, 0.5, 0.75, 1.0]
# Probit boundaries between grades: geometric midpoints of adjacent master-scale PDs.
B = ndtri(np.sqrt(PD_K[:-1] * PD_K[1:]))
LO = ndtri(PD_K / math.sqrt(1.5))   # each obligor's TTC PD is log-uniform on [PD_k/sqrt1.5, PD_k*sqrt1.5)
HI = ndtri(PD_K * math.sqrt(1.5))


def shares(centre: int) -> np.ndarray:
    w = np.exp(-((np.arange(1, K + 1) - centre) ** 2) / 32.0)
    return w / w.sum()


def population(centre: int, m: int = 400):
    """Obligors' TTC probits c and weights: within each TTC grade, log PD uniform (midpoint rule)."""
    sh = shares(centre)
    cs, ws = [], []
    u = (np.arange(m) + 0.5) / m
    for k in range(K):
        logpd = math.log(PD_K[k] / math.sqrt(1.5)) + u * math.log(1.5)
        cs.append(ndtri(np.exp(logpd)))
        ws.append(np.full(m, sh[k] / m))
    return np.concatenate(cs), np.concatenate(ws)


def grade_of(s: np.ndarray) -> np.ndarray:
    return np.searchsorted(B, s, side="right")  # 0..16


def reference(centre: int, rho: float, lam: float, nz: int = 64, n_obl: int = 10_000):
    """Per grade: pseudo-true (rho, PD) over the periods in which the grade is populated.

    A period in which the grade's expected count n_obl * w_g(z) is below 1 is treated as empty: it carries
    no default count, so it drops out of a fit. The moments are taken over the populated periods (each
    period weighted equally), and the probability of an empty period is reported.
    """
    z, wz = roots_hermitenorm(nz)
    wz = wz / wz.sum()
    c, wc = population(centre)
    lwc = np.log(wc)
    X = np.full((K, nz), np.nan)
    W = np.zeros((K, nz))
    for j, zj in enumerate(z):
        x = (c - math.sqrt(rho) * zj) / math.sqrt(1.0 - rho)
        s = (1.0 - lam) * c + lam * x
        g = grade_of(s)
        lp = log_ndtr(x)
        for k in range(K):
            m = g == k
            if not m.any():
                continue
            W[k, j] = wc[m].sum()
            X[k, j] = ndtri_exp(logsumexp(lwc[m] + lp[m]) - math.log(W[k, j]))
    k_ = 1.0 - lam + lam / math.sqrt(1.0 - rho)
    v_narrow = rho * (1.0 - lam) ** 2 / ((1.0 - rho) * k_ * k_)
    rows = []
    for g in range(K):
        pop = W[g] * n_obl >= 1.0
        p_empty = float(wz[~pop].sum())
        if p_empty > 0.5:
            rows.append((g + 1, float("nan"), float("nan"), v_narrow / (1 + v_narrow), float(W[g] @ wz), p_empty))
            continue
        q = wz[pop] / wz[pop].sum()
        mu = float(X[g, pop] @ q)
        v = float(((X[g, pop] - mu) ** 2) @ q)
        rows.append((g + 1, v / (1.0 + v), float(ndtr(mu / math.sqrt(1.0 + v))), v_narrow / (1.0 + v_narrow),
                     float(W[g] @ wz), p_empty))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).with_name("reference.csv")))
    a = ap.parse_args()
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["scale", "rho", "lambda", "grade", "master_pd", "rho_pseudo", "pd_pseudo",
                    "rho_narrow", "share_mean", "p_empty"])
        for name, (centre, rho) in SCALES.items():
            for lam in LAMBDAS:
                for g, rp, pp, rn, sm, pe in reference(centre, rho, lam):
                    w.writerow([name, rho, lam, g, f"{PD_K[g-1]:.6g}", f"{rp:.6g}", f"{pp:.6g}", f"{rn:.6g}",
                                f"{sm:.6g}", f"{pe:.6g}"])


if __name__ == "__main__":
    main()
