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

Arm I (idiosyncratic momentum) has its own closed form (narrow_idio) and finite-grade quadrature
(reference_idio), written to reference_idio.csv; --mc-check also checks both against a Monte Carlo of the
grade construction (reference_idio_mc.csv). See the "Arm I" comment below.

    uv run --no-project --with scipy==1.16.2 python studies/pit-ttc-rho/reference.py [--out FILE] [--mc-check]
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


# --- Arm I: idiosyncratic momentum (owner, 2026-10-12) -------------------------------------------------------------
# A share kappa of the idiosyncratic variance is an observable driver eta_it (AR(1), persistence psi, stationary
# N(0, 1)): x_i = (c_i - sqrt(rho) z - sqrt((1 - rho) kappa) eta_i) / sqrt((1 - rho)(1 - kappa)). The rating tracks
# eta with its own lag: N_it = (1 - thI) N_i,t-1 + thI eta_it, as M_t tracks Z_t at thS. Inside a narrow grade the
# PIT probit is a s_g / k + a sqrt(rho)(beta M - z) + a sqrt((1 - rho) kappa)(beta N - eta), a = 1/sqrt((1-rho)(1-kappa)),
# k = 1 - lam + lam a, beta = lam a / k. The idiosyncratic part varies across the grade's obligors and diversifies:
# the grade's default rate is Phi of the systematic part over sqrt(D), D = 1 + a^2 (1 - rho) kappa V_I, with
# V_I = Var(beta N - eta) = beta^2 Var N - 2 beta Cov(N, eta) + 1, Var N = thI^2/(1-(1-thI)^2) (1+(1-thI)psi)/(1-(1-thI)psi),
# Cov(N, eta) = thI / (1 - (1 - thI) psi). So v = a^2 rho (1 - 2 beta thS + beta^2 thS/(2 - thS)) / D and
# rho_eff = v / (1 + v). Checked against a Monte Carlo of the grade construction (--mc-check).
KAPPA, PSI = 0.5, 0.8
RULES_I = [(0.0, 1.0, 1.0)] + [(lam, ts, ti) for lam in (0.5, 1.0) for ts in (1.0, 0.25) for ti in (1.0, 0.5, 0.25)]


def narrow_idio(rho: float, lam: float, thS: float, thI: float, kappa: float = KAPPA, psi: float = PSI) -> float:
    a = 1.0 / math.sqrt((1.0 - rho) * (1.0 - kappa))
    k = 1.0 - lam + lam * a
    beta = lam * a / k
    q = 1.0 - 2.0 * beta * thS + beta * beta * thS / (2.0 - thS)
    var_n = thI * thI / (1.0 - (1.0 - thI) ** 2) * (1.0 + (1.0 - thI) * psi) / (1.0 - (1.0 - thI) * psi)
    cov = thI / (1.0 - (1.0 - thI) * psi)
    d = 1.0 + a * a * (1.0 - rho) * kappa * (beta * beta * var_n - 2.0 * beta * cov + 1.0)
    v = a * a * rho * q / d
    return v / (1.0 + v)


def reference_idio(centre: int, rho: float, lam: float, thS: float, thI: float, nz: int = 24, m: int = 100,
                   nn: int = 641, n_obl: int = 10_000):
    """Arm I's finite-grade pseudo-true (rho, PD) per grade, as reference() does for the main arm.

    Quadrature over (z_t, M_{t-1}); the obligor population over (c, N) with c as in population() and N on nn
    equally spaced nodes over +-8 standard deviations of its stationary law, each weighted by the normal
    density at the node, normalised to sum to 1 (a rectangle rule). The driver eta is integrated out exactly given N:
    eta | N ~ N(mu N, sig2), mu = Cov / Var N, sig2 = 1 - Cov^2 / Var N, so an obligor's conditional PD given
    (c, N, z) is Phi((c - sqrt(rho) z - sqrt((1 - rho) kappa) mu N) / sqrt((1 - rho)(1 - kappa) + (1 - rho) kappa sig2)).
    """
    z, wz = roots_hermitenorm(nz)
    wz = wz / wz.sum()
    if thS >= 1.0 or lam == 0.0:
        mprev, wm = np.array([0.0]), np.array([1.0])
    else:
        mprev, wm = z * math.sqrt(thS / (2.0 - thS)), wz
    var_n = thI * thI / (1.0 - (1.0 - thI) ** 2) * (1.0 + (1.0 - thI) * PSI) / (1.0 - (1.0 - thI) * PSI)
    cov = thI / (1.0 - (1.0 - thI) * PSI)
    mu, sig2 = cov / var_n, max(1.0 - cov * cov / var_n, 0.0)
    c1, wc1 = population(centre, m)
    u = np.linspace(-8.0, 8.0, nn)             # uniform in N, so the tails (which drive low-PD grades) are kept
    wn = np.exp(-0.5 * u * u)
    wn /= wn.sum()
    n1 = math.sqrt(var_n) * u
    c = np.repeat(c1, nn)
    n_ = np.tile(n1, c1.size)
    wc = np.repeat(wc1, nn) * np.tile(wn, c1.size)
    sq, ik = math.sqrt((1.0 - rho) * (1.0 - KAPPA)), math.sqrt((1.0 - rho) * KAPPA)
    den = math.sqrt(sq * sq + ik * ik * sig2)
    pts_w, X, W = [], [], []
    for zj, wj in zip(z, wz):
        p = ndtr((c - math.sqrt(rho) * zj - ik * mu * n_) / den) * wc
        for mj, wmj in zip(mprev, wm):
            M = (1.0 - thS) * mj + thS * zj
            s = (1.0 - lam) * c + lam * (c - math.sqrt(rho) * M - ik * n_) / sq
            g = grade_of(s)
            ws = np.bincount(g, weights=wc, minlength=K)
            ds = np.bincount(g, weights=p, minlength=K)
            with np.errstate(divide="ignore", invalid="ignore"):
                xs = np.where(ws > 0, ndtri(np.maximum(ds / np.where(ws > 0, ws, 1.0), 1e-300)), np.nan)
            pts_w.append(wj * wmj); X.append(xs); W.append(ws)
    pw, X, W = np.array(pts_w), np.array(X).T, np.array(W).T
    rn = narrow_idio(rho, lam, thS, thI)
    rows = []
    for g in range(K):
        pop = W[g] * n_obl >= 1.0
        p_empty = float(pw[~pop].sum())
        if p_empty > 0.5:
            rows.append((g + 1, float("nan"), float("nan"), rn, float(W[g] @ pw), p_empty))
            continue
        q = pw[pop] / pw[pop].sum()
        mu_x = float(X[g, pop] @ q)
        v = float(((X[g, pop] - mu_x) ** 2) @ q)
        rows.append((g + 1, v / (1.0 + v), float(ndtr(mu_x / math.sqrt(1.0 + v))), rn, float(W[g] @ pw), p_empty))
    return rows


def mc_idio(centre: int, rho: float, lam: float, thS: float, thI: float, nobl: int = 150_000, nz: int = 20, seed: int = 1):
    """Pseudo-true rho per grade from a fixed sample of obligors (c, eta, N), by quadrature over (z, M_{t-1})."""
    rng = np.random.default_rng(seed)
    sh = shares(centre)
    k = rng.choice(K, size=nobl, p=sh)
    c = ndtri(PD_K[k] * np.exp((rng.random(nobl) - 0.5) * math.log(1.5)))
    var_n = thI * thI / (1.0 - (1.0 - thI) ** 2) * (1.0 + (1.0 - thI) * PSI) / (1.0 - (1.0 - thI) * PSI)
    cov = thI / (1.0 - (1.0 - thI) * PSI)
    eta = rng.standard_normal(nobl)
    n_ = cov * eta + math.sqrt(max(var_n - cov * cov, 0.0)) * rng.standard_normal(nobl)
    z, wz = roots_hermitenorm(nz)
    wz = wz / wz.sum()
    mp, wm = (z * math.sqrt(thS / (2.0 - thS)), wz) if thS < 1.0 else (np.array([0.0]), np.array([1.0]))
    sq, ik = math.sqrt((1.0 - rho) * (1.0 - KAPPA)), math.sqrt((1.0 - rho) * KAPPA)
    xs = {g: [] for g in range(K)}
    ws = {g: [] for g in range(K)}
    for zj, wj in zip(z, wz):
        p = ndtr((c - math.sqrt(rho) * zj - ik * eta) / sq)
        for mj, wmj in zip(mp, wm):
            m = (1.0 - thS) * mj + thS * zj
            s = (1.0 - lam) * c + lam * (c - math.sqrt(rho) * m - ik * n_) / sq
            g = np.searchsorted(B, s, side="right")
            cnt = np.bincount(g, minlength=K)
            dsum = np.bincount(g, weights=p, minlength=K)
            for gg in range(K):
                if cnt[gg] >= 5e-4 * nobl:
                    xs[gg].append(ndtri(dsum[gg] / cnt[gg]))
                    ws[gg].append(wj * wmj)
    out = {}
    for gg in range(K):
        if len(ws[gg]) < 0.5 * len(z) * len(mp):
            continue
        w = np.array(ws[gg]) / sum(ws[gg])
        x = np.array(xs[gg])
        mu = w @ x
        v = w @ (x - mu) ** 2
        out[gg + 1] = v / (1.0 + v)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).with_name("reference.csv")))
    ap.add_argument("--mc-check", action="store_true", help="also check arm I's closed form by Monte Carlo (minutes)")
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
    with open(Path(a.out).with_name("reference_idio.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["scale", "rho", "lambda", "theta_s", "theta_i", "kappa", "psi", "grade", "master_pd",
                    "rho_pseudo", "pd_pseudo", "rho_narrow", "share_mean", "p_empty"])
        for name, (centre, rho) in SCALES.items():
            for lam, ts, ti in RULES_I:
                for g, rp, pp, rn, sm, pe in reference_idio(centre, rho, lam, ts, ti):
                    w.writerow([name, rho, lam, ts, ti, KAPPA, PSI, g, f"{PD_K[g-1]:.6g}", f"{rp:.6g}", f"{pp:.6g}",
                                f"{rn:.6g}", f"{sm:.6g}", f"{pe:.6g}"])
    if a.mc_check:
        inf = {"L": range(9, 15), "M": range(8, 17), "H": range(9, 17)}
        with open(Path(a.out).with_name("reference_idio_mc.csv"), "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["scale", "lambda", "theta_s", "theta_i", "rho_narrow", "mc_min", "mc_max", "grades"])
            for name, (centre, rho) in SCALES.items():
                for lam, ts, ti in RULES_I:
                    r = mc_idio(centre, rho, lam, ts, ti)
                    vals = [r[g] for g in inf[name] if g in r]
                    lo, hi = (min(vals), max(vals)) if vals else (float("nan"), float("nan"))
                    w.writerow([name, lam, ts, ti, f"{narrow_idio(rho, lam, ts, ti):.6g}", f"{lo:.6g}", f"{hi:.6g}", len(vals)])


if __name__ == "__main__":
    main()
