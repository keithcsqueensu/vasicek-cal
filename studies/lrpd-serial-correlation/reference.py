# SPDX-License-Identifier: Apache-2.0
"""S-35 (long-run PD under an AR(1) factor): reference numbers from the model and pinned results only.

No AR(1) panel is simulated or fitted here. The inputs are the model and the pinned phi = 0 coverage
(tests/golden/recovery/summary.csv for the binomial MLE, tests/golden/vasicek_rate/recovery_summary.csv
for the Vasicek-rate MLE). Under Z_t = phi Z_{t-1} + sqrt(1 - phi^2) e_t every period keeps its phi = 0
distribution, so each estimator's target is unchanged (the expected log-likelihood of a period depends
only on its marginal): what changes is the variance of a mean over T dependent periods.

  binomial  The iid binomial-mixture MLE in (logit PD, logit rho). Per period, the score s(d) of the
            marginal log-likelihood; I = E[s s'] (the information the Wald SE assumes). Under AR(1) the
            lag-k covariance of the scores is C_k = sum_{j >= 1} a_j a_j' phi^(jk), a_j the normalised
            Hermite coefficients of h(z) = E[s(d) | Z = z]: the binomial noise is independent across
            periods, so only h carries the dependence. Sandwich V = I^-1 (I + 2 sum_{k<T} (1 - k/T) C_k)
            I^-1; the inflation f_a = V_aa / (I^-1)_aa is how much the true variance of the estimate
            exceeds the variance its SE assumes.
  average   The simple average of d_t / n: mean PD exactly; per-period variance Var p(Z) + E[p(1-p)]/n,
            lag-k covariance sum_j b_j^2 phi^(jk), b_j the Hermite coefficients of p(z).
  rate      The Vasicek-rate MLE on the n -> infinity rates (probit rates x_t Gaussian AR(1), variance
            sigma^2 = rho / (1 - rho)): mu-hat inflates by f_mu = 1 + 2 sum (1 - k/T) phi^k, sigma^2-hat by
            f_s = 1 + 2 sum (1 - k/T) phi^(2k), and E[sigma^2-hat] = sigma^2 (1 - f_mu / T) exactly.

Coverage mapping: a scenario's pinned phi = 0 coverage c0 gives its effective critical value
z0 = Phi^-1((1 + c0)/2); with inflation f and an added bias of delta standard deviations, the predicted
coverage is Phi(z0/sqrt f - delta) - Phi(-z0/sqrt f - delta). delta is used for rho only:
delta = ((f_mu - 1)/T) / sqrt(2 f_s / T) in log sigma^2, the extra downward bias from the dependent mean.

    uv run --no-project --with scipy==1.16.2 python studies/lrpd-serial-correlation/reference.py [--out FILE]
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln, ndtr, ndtri

REPO = Path(__file__).resolve().parents[2]
PHIS = (0.3, 0.6)
PDS = (0.001, 0.01, 0.05)
RHOS = (0.02, 0.12, 0.24)
PERIODS = (20, 40, 100)
OBLIGORS = (100, 1000, 10000)
HERMITE = 60  # Hermite terms; the unexplained share of Var h(Z) is reported
Z = np.linspace(-9.0, 9.0, 3601)  # trapezoid grid for the factor: spectrally accurate for smooth integrands
WZ = np.exp(-0.5 * Z * Z) / math.sqrt(2 * math.pi) * (Z[1] - Z[0])
WZ[0] *= 0.5
WZ[-1] *= 0.5


def scenario(sid: int) -> tuple[float, float, int, int]:
    i_n, i_t, i_rho, i_pd = sid % 3, (sid // 3) % 3, (sid // 9) % 3, sid // 27
    return PDS[i_pd], RHOS[i_rho], PERIODS[i_t], OBLIGORS[i_n]


def cond_pd(pd: float, rho: float) -> np.ndarray:
    return ndtr((ndtri(pd) - math.sqrt(rho) * Z) / math.sqrt(1 - rho))


def hermite(j_max: int) -> np.ndarray:
    """Normalised probabilists' Hermite polynomials He_j / sqrt(j!) on Z, j = 0..j_max."""
    h = np.empty((j_max + 1, Z.size))
    h[0] = 1.0
    h[1] = Z
    for j in range(1, j_max):
        h[j + 1] = (Z * h[j] * math.sqrt(1.0 / (j + 1))) - h[j - 1] * math.sqrt(j / (j + 1))
    return h


HE = hermite(HERMITE)


def binom_given_z(n: int, p: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Binom(d; n, p(z)), shape (len(d), len(z))."""
    lp = np.log(np.clip(p, 1e-300, 1.0))
    lq = np.log1p(-np.clip(p, 0.0, 1.0 - 1e-16))
    lc = gammaln(n + 1) - gammaln(d + 1) - gammaln(n - d + 1)
    return np.exp(lc[:, None] + d[:, None] * lp[None, :] + (n - d)[:, None] * lq[None, :])


def support(n: int, pd: float, rho: float) -> np.ndarray:
    """The counts d whose marginal probability exceeds 1e-13 (the rest carry no measurable weight)."""
    d = np.arange(n + 1, dtype=float)
    marg = (binom_given_z(n, cond_pd(pd, rho), d) * WZ).sum(axis=1)
    keep = np.nonzero(marg > 1e-13)[0]
    return d[keep[0]:keep[-1] + 1]


def binomial_inflation(pd: float, rho: float, n: int, phis, periods):
    """f_PD, f_rho of the iid binomial MLE (logit axes) per (phi, T), and the Hermite remainder."""
    d = support(n, pd, rho)
    u = np.array([math.log(pd / (1 - pd)), math.log(rho / (1 - rho))])
    step = 1e-4

    def logf(uu):
        p_ = 1 / (1 + math.exp(-uu[0]))
        r_ = 1 / (1 + math.exp(-uu[1]))
        return np.log((binom_given_z(n, cond_pd(p_, r_), d) * WZ).sum(axis=1))

    s = np.empty((2, d.size))
    for a in range(2):
        e = np.zeros(2)
        e[a] = step
        s[a] = (logf(u + e) - logf(u - e)) / (2 * step)
    k_dz = binom_given_z(n, cond_pd(pd, rho), d)  # f(d | z)
    f_d = (k_dz * WZ).sum(axis=1)
    info = (s[:, None, :] * s[None, :, :] * f_d).sum(axis=2)
    h = s @ k_dz  # E[s | z], shape (2, len(Z))
    a = (h[:, None, :] * HE[None, :, :] * WZ).sum(axis=2)  # (2, J+1)
    var_h = ((h - (h * WZ).sum(axis=1, keepdims=True)) ** 2 * WZ).sum(axis=1)
    remainder = 1 - (a[:, 1:] ** 2).sum(axis=1) / var_h
    inv = np.linalg.inv(info)
    out = {}
    for phi in phis:
        for t in periods:
            b = info.copy()
            for k in range(1, t):
                ck = sum(np.outer(a[:, j], a[:, j]) * phi ** (j * k) for j in range(1, HERMITE + 1))
                b += 2 * (1 - k / t) * ck
            v = inv @ b @ inv
            out[(phi, t)] = (v[0, 0] / inv[0, 0], v[1, 1] / inv[1, 1])
    return out, remainder


def average_inflation(pd: float, rho: float, n: int, phi: float, t: int) -> float:
    p = cond_pd(pd, rho)
    b = (p * HE * WZ).sum(axis=1)
    var = ((p - pd) ** 2 * WZ).sum() + (p * (1 - p) * WZ).sum() / n
    lr = var + 2 * sum((1 - k / t) * sum(b[j] ** 2 * phi ** (j * k) for j in range(1, HERMITE + 1)) for k in range(1, t))
    return lr / var


def f_ar(phi: float, t: int, power: int = 1) -> float:
    return 1 + 2 * sum((1 - k / t) * phi ** (power * k) for k in range(1, t))


def coverage(c0: float, f: float, delta: float = 0.0) -> float:
    c0 = min(max(c0, 1e-6), 1 - 1e-6)
    z0 = float(ndtri((1 + c0) / 2))
    return float(ndtr(z0 / math.sqrt(f) - delta) - ndtr(-z0 / math.sqrt(f) - delta))


def read_summary(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "studies/lrpd-serial-correlation/reference.csv"))
    args = ap.parse_args()
    pinned = {int(r["scenario"]): r for r in read_summary(REPO / "tests/golden/recovery/summary.csv")}
    rate_rows = read_summary(REPO / "tests/golden/vasicek_rate/recovery_summary.csv")
    rows = []
    cache = {}
    for sid in range(81):
        pd, rho, t, n = scenario(sid)
        if (pd, rho, n) not in cache:
            cache[(pd, rho, n)] = binomial_inflation(pd, rho, n, PHIS, PERIODS)
        infl, rem = cache[(pd, rho, n)]
        r = pinned[sid]
        for phi in PHIS:
            f_pd, f_rho = infl[(phi, t)]
            f_mu, f_s = f_ar(phi, t), f_ar(phi, t, 2)
            delta = ((f_mu - 1) / t) / math.sqrt(2 * f_s / t)
            rows.append({
                "scenario": sid, "pd": pd, "rho": rho, "periods": t, "obligors": n, "phi": phi,
                "f_pd": f_pd, "f_rho": f_rho, "f_avg": average_inflation(pd, rho, n, phi, t),
                "hermite_remainder_pd": rem[0], "hermite_remainder_rho": rem[1],
                "rho_extra_bias_rel": -(f_mu - 1) / t * (1 - rho), "rho_delta": delta,
                "wald_pd_c0": float(r["pd_coverage"]), "wald_pd_pred": coverage(float(r["pd_coverage"]), f_pd),
                "wald_rho_c0": float(r["rho_coverage"]),
                "wald_rho_pred": coverage(float(r["rho_coverage"]), f_rho, delta),
                "profile_pd_c0": float(r["pd_profile_coverage"]),
                "profile_pd_pred": coverage(float(r["pd_profile_coverage"]), f_pd),
                "profile_rho_c0": float(r["rho_profile_coverage"]),
                "profile_rho_pred": coverage(float(r["rho_profile_coverage"]), f_rho, delta),
            })
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in row.items()})
    rate = []
    for r in rate_rows:
        pd, rho, t = float(r["pd"]), float(r["rho"]), int(r["periods"])
        reps = int(r["replicates"])
        mu, s2 = float(ndtri(pd)) * math.sqrt(1 + rho / (1 - rho)), rho / (1 - rho)
        arg = mu / math.sqrt(1 + s2)
        g_mu = math.exp(-0.5 * arg * arg) / math.sqrt(2 * math.pi) / math.sqrt(1 + s2)
        g_s = math.exp(-0.5 * arg * arg) / math.sqrt(2 * math.pi) * (-0.5 * mu) * (1 + s2) ** -1.5
        v_mu, v_s = g_mu * g_mu * s2, g_s * g_s * 2 * s2 * s2  # per-period delta-method terms
        for phi in PHIS:
            f_mu, f_s = f_ar(phi, t), f_ar(phi, t, 2)
            delta = ((f_mu - 1) / t) / math.sqrt(2 * f_s / t)
            f_pd = (v_mu * f_mu + v_s * f_s) / (v_mu + v_s)
            c_pd, c_rho = int(r["pd_covered"]) / reps, int(r["rho_covered"]) / reps
            rate.append({"pd": pd, "rho": rho, "periods": t, "phi": phi, "f_mu": f_mu, "f_s": f_s, "f_pd": f_pd,
                         "sigma2_ratio": 1 - f_mu / t, "sigma2_ratio_phi0": 1 - 1 / t,
                         "profile_pd_c0": c_pd, "profile_pd_pred": coverage(c_pd, f_pd),
                         "profile_rho_c0": c_rho, "profile_rho_pred": coverage(c_rho, f_s, delta)})
    rate_out = Path(args.out).with_name(Path(args.out).stem + "_rate.csv")
    with open(rate_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rate[0]), lineterminator="\n")
        w.writeheader()
        for row in rate:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in row.items()})
    print("rate MLE (n -> infinity): T, phi, f_mu, f_s, E[sigma2-hat]/sigma2 (phi = 0: 1 - 1/T)")
    for t in PERIODS:
        for phi in PHIS:
            print(f"  T {t:3d} phi {phi}: f_mu {f_ar(phi, t):.3f}, f_s {f_ar(phi, t, 2):.3f}, "
                  f"ratio {1 - f_ar(phi, t) / t:.4f} (phi 0: {1 - 1 / t:.4f})")
    print(f"wrote {args.out} and {rate_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
