# SPDX-License-Identifier: Apache-2.0
"""S-36 (rating philosophy and fitted correlation): reference values from the model only.

No panel is simulated or fitted here. Everything follows from the DGP's definition
(PREDICTION.md, "The DGP") by quadrature over the factor and the obligor population.

Model. Obligor i has a TTC probit c_i = Phi^-1(TTC PD_i). Its PIT probit in a period with factor z is
x_i(z) = (c_i - sqrt(rho) z) / sqrt(1 - rho) (engine convention: a higher z is better). Its rating target
is the hybrid y_i(z) = (1 - lam) c_i + lam x_i(z) (lam = 0 TTC, lam = 1 PIT), and the rating moves only
part of the way towards it each period (rating momentum, speed theta): the score uses the smoothed factor
M_t = (1 - theta) M_{t-1} + theta z_t in place of z_t, s_i = (1 - lam) c_i + lam (c_i - sqrt(rho) M_t) /
sqrt(1 - rho). theta = 1 is no momentum. M is stationary: M_{t-1} ~ N(0, theta / (2 - theta)),
independent of z_t. Grade g holds the obligors whose score lies in [b_{g-1}, b_g), the master scale's
probit boundaries (open-ended at both ends of the scale).

Narrow-grade limit (closed form). With a = 1 / sqrt(1 - rho), k = 1 - lam + lam a and beta = lam a / k,
inside a grade the PIT probit is a s_g / k + a sqrt(rho) (beta M_t - z_t). So the grade's probit default
rate is Gaussian with variance v = rho / (1 - rho) * (1 - 2 beta theta + beta^2 theta / (2 - theta)),
rho_eff = v / (1 + v), and lag-j autocovariance (1 - theta)^j beta theta (beta / (2 - theta) - 1) rho /
(1 - rho), negative: a rating that lags the cycle leaves better obligors in a grade after a bad year.
At theta = 1, v = rho (1 - lam)^2 / ((1 - rho) k^2): rho at lam = 0, 0 at lam = 1.

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


def narrow(rho: float, lam: float, theta: float):
    a = 1.0 / math.sqrt(1.0 - rho)
    k = 1.0 - lam + lam * a
    beta = lam * a / k
    q = 1.0 - 2.0 * beta * theta + beta * beta * theta / (2.0 - theta)
    v = rho / (1.0 - rho) * q
    g1 = (1.0 - theta) * beta * theta * (beta / (2.0 - theta) - 1.0) * rho / (1.0 - rho)
    return v / (1.0 + v), (g1 / v if v > 0 else 0.0)


def reference(centre: int, rho: float, lam: float, theta: float, nz: int = 48, n_obl: int = 10_000):
    """Per grade: pseudo-true (rho, PD) over the periods in which the grade is populated.

    Quadrature over (z_t, M_{t-1}), independent normals. A period in which the grade's expected count
    n_obl * w_g is below 1 is treated as empty: it carries no count and drops out of a fit. The moments are
    over the populated periods (each weighted equally); the probability of an empty period is reported.
    """
    z, wz = roots_hermitenorm(nz)
    wz = wz / wz.sum()
    if theta >= 1.0 or lam == 0.0:
        mprev, wm = np.array([0.0]), np.array([1.0])
    else:
        mprev, wm = z * math.sqrt(theta / (2.0 - theta)), wz
    c, wc = population(centre)
    lwc = np.log(wc)
    pts_w, X, W = [], [], []
    for zj, wj in zip(z, wz):
        x = (c - math.sqrt(rho) * zj) / math.sqrt(1.0 - rho)
        lp = log_ndtr(x)
        for mj, wmj in zip(mprev, wm):
            M = (1.0 - theta) * mj + theta * zj
            s = (1.0 - lam) * c + lam * (c - math.sqrt(rho) * M) / math.sqrt(1.0 - rho)
            g = grade_of(s)
            xs, ws = np.full(K, np.nan), np.zeros(K)
            for k in range(K):
                m = g == k
                if m.any():
                    ws[k] = wc[m].sum()
                    xs[k] = ndtri_exp(logsumexp(lwc[m] + lp[m]) - math.log(ws[k]))
            pts_w.append(wj * wmj); X.append(xs); W.append(ws)
    pw, X, W = np.array(pts_w), np.array(X).T, np.array(W).T
    rn, r1 = narrow(rho, lam, theta)
    rows = []
    for g in range(K):
        pop = W[g] * n_obl >= 1.0
        p_empty = float(pw[~pop].sum())
        if p_empty > 0.5:
            rows.append((g + 1, float("nan"), float("nan"), rn, r1, float(W[g] @ pw), p_empty))
            continue
        q = pw[pop] / pw[pop].sum()
        mu = float(X[g, pop] @ q)
        v = float(((X[g, pop] - mu) ** 2) @ q)
        rows.append((g + 1, v / (1.0 + v), float(ndtr(mu / math.sqrt(1.0 + v))), rn, r1, float(W[g] @ pw), p_empty))
    return rows


RULES = [(0.0, 1.0)] + [(lam, th) for lam in (0.25, 0.5, 0.75, 1.0) for th in (1.0, 0.5, 0.25)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).with_name("reference.csv")))
    a = ap.parse_args()
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["scale", "rho", "lambda", "theta", "grade", "master_pd", "rho_pseudo", "pd_pseudo",
                    "rho_narrow", "acf1_narrow", "share_mean", "p_empty"])
        for name, (centre, rho) in SCALES.items():
            for lam, th in RULES:
                for g, rp, pp, rn, r1, sm, pe in reference(centre, rho, lam, th):
                    w.writerow([name, rho, lam, th, g, f"{PD_K[g-1]:.6g}", f"{rp:.6g}", f"{pp:.6g}", f"{rn:.6g}",
                                f"{r1:.6g}", f"{sm:.6g}", f"{pe:.6g}"])


if __name__ == "__main__":
    main()
