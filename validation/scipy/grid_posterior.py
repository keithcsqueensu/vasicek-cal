# SPDX-License-Identifier: Apache-2.0
"""Replicate the grid-Bayesian estimator with numpy/scipy (M3; engine/posterior.hpp, D-161).

Independent of the engine's quadrature, each binomial-mixture probability
    P(d | n, PD, rho) = C(n, d) E[p(Z)^d (1 - p(Z))^(n - d)],   p(z) = Phi((Phi^-1(PD) - sqrt(rho) z) / sqrt(1 - rho)),
is computed here by scipy's adaptive quad (for the Jeffreys table check, where a central difference
amplifies every error) or by a fixed 200-node Gauss-Hermite rule vectorised in numpy (for the
posteriors), where the engine uses its adaptive parity rule. With them, this script:
  - recomputes the Jeffreys table of tests/golden/posterior/jeffreys.csv (log sqrt det I_u, I_u the Fisher
    information of one period in logit u, the score by central differences with step 1e-4, over the
    counts with P(d) >= 1e-14 x max) and compares it to TOL_SCIPY_JEFFREYS_LOG_ABS;
  - re-runs the estimator on every panel of tests/golden/posterior/posterior.csv: the posterior on the
    parity grid (61 x 41 logit over the box) with the flat or Jeffreys prior and cell-uniform masses;
    checks that the resolution rule (4 points per SD) fails there exactly when the engine refined and
    holds on the engine's final grid (read from the file: its construction uses the MLE's Hessian SE);
    and recomputes the equal-tailed and HPD intervals on that grid, to TOL_SCIPY_POSTERIOR_END_U.

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/grid_posterior.py
"""

import csv
import math
import sys
from pathlib import Path

import numpy as np
import tomllib
from numpy.polynomial.hermite_e import hermegauss
from scipy import integrate, special, stats

ROOT = Path(__file__).resolve().parents[2]
GOLD = ROOT / "tests" / "golden" / "posterior"
PD_BOX, RHO_BOX = (1e-4, 0.2), (1e-3, 0.5)
NODES, WEIGHTS = hermegauss(200)
LOGW = np.log(WEIGHTS / math.sqrt(2 * math.pi))
STEP, FLOOR = 1e-4, 1e-14
POINTS_PER_SD, REFINE_PER_SD, SPAN_SD, MAX_REF, OUTSIDE, LEVEL = 4.0, 5.0, 8.0, 3, 1e-6, 0.95
REFINED, REFUSED = 1, 2


def logit(v):
    return np.log(v) - np.log1p(-v)


def inv_logit(u):
    return 1.0 / (1.0 + np.exp(-u))


def log_prob(n, d, pd, rho):
    """log P(d | n, PD, rho) for arrays pd, rho (same shape) and a 1-D array of counts d."""
    pd, rho = np.asarray(pd, float)[..., None, None], np.asarray(rho, float)[..., None, None]
    d = np.asarray(d, float)[:, None]
    x = (stats.norm.ppf(pd) - np.sqrt(rho) * NODES) / np.sqrt(1 - rho)
    lp, lq = special.log_ndtr(x), special.log_ndtr(-x)
    lc = special.gammaln(n + 1) - special.gammaln(d + 1) - special.gammaln(n - d + 1)
    return lc[..., 0] + special.logsumexp(d * lp + (n - d) * lq + LOGW, axis=-1)


def log_prob_quad(n, d, pd, rho):
    """log P(d | n, PD, rho) by adaptive quad, split at the integrand's peak."""
    c, s, s1 = stats.norm.ppf(pd), math.sqrt(rho), math.sqrt(1 - rho)
    zs = np.linspace(-12, 12, 481)
    x = (c - s * zs) / s1
    h = d * special.log_ndtr(x) + (n - d) * special.log_ndtr(-x) - zs * zs / 2
    top, zm = float(h.max()), float(zs[h.argmax()])

    def f(z):
        xx = (c - s * z) / s1
        return math.exp(d * special.log_ndtr(xx) + (n - d) * special.log_ndtr(-xx) - z * z / 2 - top)

    v = sum(integrate.quad(f, a, b, limit=400, epsabs=0.0, epsrel=1e-13)[0] for a, b in ((-40.0, zm), (zm, 40.0)))
    lc = special.gammaln(n + 1) - special.gammaln(d + 1) - special.gammaln(n - d + 1)
    return lc + top + math.log(v / math.sqrt(2 * math.pi))


def jeffreys_quad(n, u0, u1):
    """log sqrt det I_u at one point, every probability by adaptive quad."""
    lps = [log_prob_quad(n, d, inv_logit(u0), inv_logit(u1)) for d in range(n + 1)]
    floor = max(lps) + math.log(FLOOR)
    i00 = i01 = i11 = 0.0
    for d, lp in enumerate(lps):
        if lp < floor:
            continue
        g0 = (log_prob_quad(n, d, inv_logit(u0 + STEP), inv_logit(u1)) - log_prob_quad(n, d, inv_logit(u0 - STEP), inv_logit(u1))) / (2 * STEP)
        g1 = (log_prob_quad(n, d, inv_logit(u0), inv_logit(u1 + STEP)) - log_prob_quad(n, d, inv_logit(u0), inv_logit(u1 - STEP))) / (2 * STEP)
        p = math.exp(lp)
        i00, i01, i11 = i00 + p * g0 * g0, i01 + p * g0 * g1, i11 + p * g1 * g1
    return 0.5 * math.log(i00 * i11 - i01 * i01)


def jeffreys(n, u0, u1):
    """log sqrt det I_u at arrays of points (u0, u1)."""
    d = np.arange(n + 1)
    lp = log_prob(n, d, inv_logit(u0), inv_logit(u1))
    keep = lp >= lp.max(axis=-1, keepdims=True) + math.log(FLOOR)
    g0 = (log_prob(n, d, inv_logit(u0 + STEP), inv_logit(u1)) - log_prob(n, d, inv_logit(u0 - STEP), inv_logit(u1))) / (2 * STEP)
    g1 = (log_prob(n, d, inv_logit(u0), inv_logit(u1 + STEP)) - log_prob(n, d, inv_logit(u0), inv_logit(u1 - STEP))) / (2 * STEP)
    p = np.where(keep, np.exp(lp), 0.0)
    i00, i01, i11 = (p * g0 * g0).sum(-1), (p * g0 * g1).sum(-1), (p * g1 * g1).sum(-1)
    return 0.5 * np.log(i00 * i11 - i01 * i01)


class Axis:
    def __init__(self, lo, hi, n):
        self.u = np.linspace(logit(lo), logit(hi), n)
        self.u[-1] = logit(hi)
        self.step = (logit(hi) - logit(lo)) / (n - 1)
        self.width = np.full(n, self.step)
        self.width[[0, -1]] = self.step / 2
        self.lo = np.concatenate([[self.u[0]], self.u[1:] - self.step / 2])
        self.hi = np.concatenate([self.u[:-1] + self.step / 2, [self.u[-1]]])


def masses(counts, n, a0, a1, prior, jt):
    U0, U1 = np.meshgrid(a0.u, a1.u, indexing="ij")
    ds, k = np.unique(counts, return_counts=True)
    ll = (log_prob(n, ds, inv_logit(U0), inv_logit(U1)) * k).sum(-1)
    if prior == "flat":
        lpri = np.log(inv_logit(U0) * (1 - inv_logit(U0))) + np.log(inv_logit(U1) * (1 - inv_logit(U1)))
    else:
        lpri = jt(U0, U1)
    lp = ll + lpri + np.log(np.outer(a0.width, a1.width))
    m = np.exp(lp - lp.max())
    return m / m.sum()


def moments(ax, mass):
    c, w = (ax.lo + ax.hi) / 2, ax.hi - ax.lo
    m1 = (mass * c).sum()
    return m1, math.sqrt(max((mass * (c * c + w * w / 12)).sum() - m1 * m1, 0.0))


def quantile(ax, mass, p):
    cum = np.concatenate([[0.0], np.cumsum(mass)])
    i = int(np.searchsorted(cum[1:], p))
    i = min(i, len(mass) - 1)
    return ax.lo[i] + (p - cum[i]) / mass[i] * (ax.hi[i] - ax.lo[i])


def cdf_at(ax, mass, u):
    frac = np.clip((u - ax.lo) / np.where(ax.hi > ax.lo, ax.hi - ax.lo, 1.0), 0.0, 1.0)
    return float((mass * frac).sum())


def hpd(ax, mass):
    """The shortest interval holding LEVEL; among candidates within 1e-9 relative of the shortest, the one
    centred closest to the median (then the lowest), as the engine chooses."""
    cand = []
    for i in range(len(mass)):
        for e in (ax.lo[i], ax.hi[i]):
            c = cdf_at(ax, mass, e)
            if c + LEVEL <= 1.0:
                cand.append((e, quantile(ax, mass, c + LEVEL)))
            if c >= LEVEL:
                cand.append((quantile(ax, mass, c - LEVEL), e))
    shortest = min(b - a for a, b in cand)
    median = quantile(ax, mass, 0.5)
    ties = [(abs((a + b) / 2 - median), a, b) for a, b in cand if b - a <= shortest * (1 + 1e-9)]
    _, a, b = min(ties)
    return a, b


def main():
    with open(ROOT / "tests" / "tolerances.toml", "rb") as f:
        tol = {k: v["value"] for k, v in tomllib.load(f).items()}
    failures = []
    # --- the Jeffreys table ---
    with open(GOLD / "jeffreys.csv", newline="") as f:
        jrows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    h = float.fromhex
    u0 = logit(np.array([h(r["pd_hex"]) for r in jrows]))
    u1 = logit(np.array([h(r["rho_hex"]) for r in jrows]))
    got = np.array([jeffreys_quad(int(jrows[0]["n"]), a, b) for a, b in zip(u0, u1)])
    want = np.array([h(r["log_prior_hex"]) for r in jrows])
    worst_j = float(np.max(np.abs(got - want)))
    if worst_j > tol["TOL_SCIPY_JEFFREYS_LOG_ABS"]:
        failures.append(f"Jeffreys table differs by {worst_j:.3g}")
    # The parity grid's table for the posterior's Jeffreys prior (bilinear in log, as the engine).
    pa0, pa1 = Axis(*PD_BOX, 61), Axis(*RHO_BOX, 41)
    P0, P1 = np.meshgrid(pa0.u, pa1.u, indexing="ij")
    table = jeffreys(100, P0, P1)

    def jt(U0, U1):
        f0 = np.clip((U0 - pa0.u[0]) / pa0.step, 0, 60)
        f1 = np.clip((U1 - pa1.u[0]) / pa1.step, 0, 40)
        i0, i1 = np.minimum(f0.astype(int), 59), np.minimum(f1.astype(int), 39)
        a, b = f0 - i0, f1 - i1
        return ((1 - a) * ((1 - b) * table[i0, i1] + b * table[i0, i1 + 1]) +
                a * ((1 - b) * table[i0 + 1, i1] + b * table[i0 + 1, i1 + 1]))

    # --- the posteriors ---
    with open(GOLD / "posterior.csv", newline="") as f:
        prows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    worst_u = 0.0
    for r in prows:
        at = f"PD {r['pd']} rho {r['rho']} rep {r['replicate']} {r['prior']}"
        counts = np.array([int(x) for x in r["counts"].split()])
        # The parity grid, then the local grid the engine reports (checked: the rule must fail on the
        # parity grid exactly when the engine refined, and hold on the reported grid).
        m = masses(counts, 100, pa0, pa1, r["prior"], jt)
        sd_p = [moments(pa0, m.sum(1))[1], moments(pa1, m.sum(0))[1]]
        needs = not (sd_p[0] >= POINTS_PER_SD * pa0.step and sd_p[1] >= POINTS_PER_SD * pa1.step)
        if needs != bool(int(r["flags"]) & REFINED):
            failures.append(f"{at}: the rule's decision on the parity grid differs")
            continue
        a0 = Axis(h(r["pd_lo_bound_hex"]), h(r["pd_hi_bound_hex"]), int(r["pd_points"]))
        a1 = Axis(h(r["rho_lo_bound_hex"]), h(r["rho_hi_bound_hex"]), int(r["rho_points"]))
        m = masses(counts, 100, a0, a1, r["prior"], jt)
        marg = [m.sum(1), m.sum(0)]
        for axis_i, (ax, mg, name) in enumerate(((a0, marg[0], "pd"), (a1, marg[1], "rho"))):
            sd = moments(ax, mg)[1]
            if not sd >= POINTS_PER_SD * ax.step:
                failures.append(f"{at}: the rule fails on the engine's final grid, axis {axis_i}")
            ends = {"et_lo": quantile(ax, mg, (1 - LEVEL) / 2), "et_hi": quantile(ax, mg, 1 - (1 - LEVEL) / 2)}
            ends["hpd_lo"], ends["hpd_hi"] = hpd(ax, mg)
            for key, u in ends.items():
                du = abs(u - logit(h(r[f"{name}_{key}_hex"])))
                worst_u = max(worst_u, du)
                if du > tol["TOL_SCIPY_POSTERIOR_END_U"]:
                    failures.append(f"{at}: {name} {key} differs by {du:.3g} in logit")
    print(f"Jeffreys table: {len(jrows)} points, worst |log difference| {worst_j:.3g}")
    print(f"posteriors: {len(prows)} fits, worst interval end difference {worst_u:.3g} in logit")
    for line in failures:
        print("FAIL " + line, file=sys.stderr)
    print("scipy agrees with the engine" if not failures else f"{len(failures)} disagreement(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
