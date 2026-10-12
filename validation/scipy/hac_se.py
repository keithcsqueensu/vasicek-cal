# SPDX-License-Identifier: Apache-2.0
"""Replicate the engine's HAC (Newey-West) sandwich SE with plain scipy (S-48, D-209).

The method (docs/methodology/hac_se.md). At the binomial-mixture MLE's estimate (PD, rho), in logit
coordinates u = (logit PD, logit rho):
    s_t   period t's score d l_t / d u, l_t = log E_Z[binom.pmf(d_t; n_t, p(Z))];
    c_t   s_t minus the mean of the s_t over the T periods;
    S     sum_t c_t c_t' + sum_{l=1}^{L} (1 - l/(L+1)) sum_{t=l}^{T-1} (c_t c_{t-l}' + c_{t-l} c_t');
    H     the panel log-likelihood's Hessian in u;
    V     H^-1 S H^-1; the HAC SE of u_a is sqrt(V_aa).
The lag rule is L = floor(4 (T/100)^(2/9)).

This script computes each piece on its own: the score as the posterior mean of d log binom.pmf / d u
(scipy.integrate.quad over the factor, around the integrand's mode, as validation/scipy/
binomial_mixture_mle.py integrates the likelihood), and H by central differences of the panel score. It
checks, on the six M1.7 validation panels (tests/golden/validation/panels.csv) at the engine's estimates
and lags (tests/golden/hac/validation.csv):
  1. the lag rule gives the golden's middle lag at T = 20;
  2. (-H)^-1 against the engine's iid SEs and correlation, which come from the engine's own analytic
     Hessian: so H is the same matrix;
  3. the HAC SEs and correlation at L = 1, the rule and twice the rule;
  4. the fits with no HAC SE are exactly those flagged kFlagGridEdge, kFlagFlatSurface or kFlagNearBound
     (bits from tests/engine_flags.json), and for each such near-bound fit this script finds the estimate
     within 2 of its own iid SEs of the box's bound in logit units.
Tolerances come from tests/tolerances.toml (TOL_SCIPY_HAC_SE_REL).

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/hac_se.py
"""

import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import tomllib
from scipy import integrate, optimize, stats

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = ROOT / "tests" / "golden"
BOX = ((1e-4, 0.2), (1e-3, 0.5))  # the recovery grid's box (tests/recovery/recovery.hpp)
NEAR_BOUND_SE = 2.0  # engine/calibrate.hpp kNearBoundSe
DROP = 50.0  # integrate where the integrand is within e^-50 of its maximum
STEP = 1e-5  # central-difference step for H, logit units


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def inv_logit(u):
    return 1.0 / (1.0 + math.exp(-u))


def period_score(u, n, d):
    """d l_t / d u for one period: E[d log binom.pmf(d; n, p(Z)) / d u | d], the posterior mean over Z."""
    pd, rho = inv_logit(u[0]), inv_logit(u[1])
    c, a, b = stats.norm.ppf(pd), math.sqrt(rho), math.sqrt(1.0 - rho)

    def x(z):
        return (c - a * z) / b

    def h(z):  # log of binom.pmf(d; n, p(z)) * phi(z)
        return stats.binom.logpmf(d, n, stats.norm.cdf(x(z))) + stats.norm.logpdf(z)

    def dlog_dx(z):  # d log binom.pmf / dp * dp/dx, in log space where p is tiny
        xz = x(z)
        log_p, log_q = stats.norm.logcdf(xz), stats.norm.logcdf(-xz)
        log_phi = stats.norm.logpdf(xz)
        return d * math.exp(log_phi - log_p) - (n - d) * math.exp(log_phi - log_q)

    # dx/dPD and dx/drho, then d/du = d/dv * v (1 - v)
    def dx_du(z):
        dpd = 1.0 / (stats.norm.pdf(c) * b)
        drho = -z / (2.0 * a * b) + x(z) / (2.0 * (1.0 - rho))
        return dpd * pd * (1.0 - pd), drho * rho * (1.0 - rho)

    z0 = (c - b * stats.norm.ppf((d + 0.5) / (n + 1))) / a
    mode = optimize.minimize_scalar(lambda z: -h(z), bracket=(z0 - 1.0, z0 + 1.0)).x
    top = h(mode)

    def edge(direction):
        z, step = mode, 0.25 * direction
        while h(z) > top - DROP:
            z, step = z + step, 2.0 * step
        return z

    lo, hi = edge(-1.0), edge(+1.0)

    def quad(f, epsabs):
        return integrate.quad(f, lo, hi, points=[mode], epsabs=epsabs, epsrel=1e-12, limit=500)[0]

    # The score integrands change sign and integrate to about 0 near the maximum, where a relative
    # tolerance cannot be met: they get an absolute one, 1e-12 of the score's scale.
    w = quad(lambda z: math.exp(h(z) - top), 0.0)
    g0 = quad(lambda z: math.exp(h(z) - top) * dlog_dx(z) * dx_du(z)[0], 1e-12 * w) / w
    g1 = quad(lambda z: math.exp(h(z) - top) * dlog_dx(z) * dx_du(z)[1], 1e-12 * w) / w
    return np.array([g0, g1])


def scores(u, panel):
    """Per-period scores in period order; each distinct (n, d) integrated once."""
    once = {obs: period_score(u, *obs) for obs in Counter(panel)}
    return np.array([once[obs] for obs in panel])


def hessian(u, panel):
    """Central differences of the panel score, symmetrised."""
    H = np.empty((2, 2))
    for b in range(2):
        e = np.zeros(2)
        e[b] = STEP
        H[:, b] = (scores(u + e, panel).sum(axis=0) - scores(u - e, panel).sum(axis=0)) / (2.0 * STEP)
    return 0.5 * (H + H.T)


def newey_west(s, lag):
    c = s - s.mean(axis=0)
    S = c.T @ c
    for l in range(1, min(lag, len(c) - 1) + 1):
        G = c[l:].T @ c[:-l]
        S += (1.0 - l / (lag + 1.0)) * (G + G.T)
    return S


def main():
    with open(ROOT / "tests" / "tolerances.toml", "rb") as f:
        tol = {k: v["value"] for k, v in tomllib.load(f).items()}["TOL_SCIPY_HAC_SE_REL"]
    panels = {}
    for r in read_csv(GOLDEN / "validation" / "panels.csv"):
        panels.setdefault(int(r["scenario"]), []).append((int(r["n"]), int(r["d"])))
    rows = read_csv(GOLDEN / "hac" / "validation.csv")
    flags = json.loads((ROOT / "tests" / "engine_flags.json").read_text())["estimate"]
    skip = sum(flags[k]["value"] for k in ("kFlagGridEdge", "kFlagFlatSurface", "kFlagNearBound"))
    near = flags["kFlagNearBound"]["value"]
    ok = True

    T = {len(p) for p in panels.values()}
    rule = {math.floor(4.0 * (t / 100.0) ** (2.0 / 9.0)) for t in T}
    lags = sorted({int(r["lag"]) for r in rows})
    passed = len(T) == 1 and len(rule) == 1 and lags == [1, *rule, 2 * min(rule)]
    ok &= passed
    print(f"1. lag rule: T = {sorted(T)}, rule {sorted(rule)}, golden lags {lags} {'ok' if passed else 'FAIL'}")

    worst_iid = worst_hac = 0.0
    computed = not_computed = 0
    for s in sorted(panels):
        mine = [r for r in rows if int(r["scenario"]) == s]
        r0 = mine[0]
        pd, rho = float.fromhex(r0["pd_hex"]), float.fromhex(r0["rho_hex"])
        u = np.array([math.log(pd / (1.0 - pd)), math.log(rho / (1.0 - rho))])
        sc = scores(u, panels[s])
        H = hessian(u, panels[s])
        A = np.linalg.inv(H)
        iid = np.linalg.inv(-H)
        iid_se = np.sqrt(np.diag(iid))
        iid_corr = iid[0, 1] / (iid_se[0] * iid_se[1])
        want = [float.fromhex(r0[k]) for k in ("iid_se_u_pd_hex", "iid_se_u_rho_hex", "iid_corr_hex")]
        d_iid = max(abs(iid_se[0] / want[0] - 1), abs(iid_se[1] / want[1] - 1), abs(iid_corr - want[2]))
        worst_iid = max(worst_iid, d_iid)
        line = [f"scenario {s}: PD {pd:.5g} rho {rho:.5g} | 2. iid SE difference {d_iid:.2g}"]
        to_bound = [min(u[a] - math.log(lo / (1 - lo)), math.log(hi / (1 - hi)) - u[a]) / iid_se[a]
                    for a, (lo, hi) in enumerate(BOX)]
        for r in mine:
            V = A @ newey_west(sc, int(r["lag"])) @ A
            se = np.sqrt(np.diag(V))
            corr = V[0, 1] / (se[0] * se[1])
            want = [float.fromhex(r[k]) for k in ("hac_se_u_pd_hex", "hac_se_u_rho_hex", "hac_corr_hex")]
            f = int(r["flags"])
            if any(math.isnan(w) for w in want):  # 4. not computed: must be a flagged fit
                not_computed += 1
                why = f & skip != 0 and all(math.isnan(w) for w in want)
                if f & near:
                    why &= min(to_bound) < NEAR_BOUND_SE
                ok &= why
                line.append(f"L {r['lag']}: not computed (flags {f}; {min(to_bound):.3g} iid SEs from the box's "
                            f"bound) {'ok' if why else 'FAIL'}")
                continue
            computed += 1
            if f & skip:
                ok = False
                line.append(f"L {r['lag']}: computed on a flagged fit (flags {f}) FAIL")
            d = max(abs(se[0] / want[0] - 1), abs(se[1] / want[1] - 1), abs(corr - want[2]))
            if not math.isfinite(d):
                d = math.inf
            worst_hac = max(worst_hac, d)
            line.append(f"L {r['lag']}: SE {se[0]:.6g}, {se[1]:.6g} (iid {iid_se[0]:.6g}, {iid_se[1]:.6g}) diff {d:.2g}")
        print(" | ".join(line))
    print(f"4. HAC SEs computed {computed}, not computed {not_computed} (flagged fits)")
    for name, worst in (("2. iid SEs from scipy's H", worst_iid), ("3. HAC SEs", worst_hac)):
        passed = worst <= tol
        ok &= passed
        print(f"{name}: worst difference {worst:.3g} (tolerance {tol:g}) {'ok' if passed else 'FAIL'}")
    print("all checks passed" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
