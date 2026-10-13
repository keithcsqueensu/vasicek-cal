# SPDX-License-Identifier: Apache-2.0
"""S-6 (Pluto-Tasche most-prudent upper bounds): exact reference numbers, from the model only.

No panel is simulated or fitted here. Every number is exact up to quadrature and root-solving error.

The model. In each period, Z_t ~ N(0, 1) independently, p(z) = Phi((Phi^-1(PD) - sqrt(rho) z) / sqrt(1 - rho)),
and d_t | Z_t ~ Bin(n, p(Z_t)) (docs/methodology/dgp.md). With n constant, the total S = sum_t d_t has the T-fold
convolution of the one-period mixed binomial PMF.

The bounds, both functions of the total D only (n is constant over the T periods):
    (a) independent: the p with P(Bin(N, p) <= D) = 1 - gamma, N = T n (one-sided Clopper-Pearson),
        here from the beta quantile, p = B^-1(gamma; D + 1, N - D); at D = 0 it is 1 - (1 - gamma)^(1/N);
    (b) one-factor: the p with P(S <= D | p, rho_a) = 1 - gamma under the model above, with rho_a assumed:
        the generating rho (oracle) or 0.12 (fixed).
P(S <= D | p, rho) falls as p rises (p(z) rises with p for every z), so each bound rises with D.

Exact coverage. bound(D) >= PD exactly when P_a(S <= D | PD) >= 1 - gamma, with P_a the bound's assumed law.
So the bound covers exactly when D >= D_min, the smallest D with P_a(S <= D | PD) >= 1 - gamma, and
    coverage = P_true(S >= D_min) = 1 - F_true(D_min - 1).
For the oracle, P_a = P_true, so coverage = 1 - F_true(D_min - 1) >= gamma by construction (F_true(D_min - 1) <
1 - gamma). No bound is solved for the coverage; the bounds at D_min - 1 and D_min are solved only to report how far
each lies from the truth in logit (a simulated bound within the solver's tolerance of the truth could flip).

Exact median ratio. The bound is increasing in D, so its median over panels is bound(D_med), D_med the smallest D
with F_true(D) >= 1/2. A simulated median over R replicates lies, with probability about 0.999, between bound(D_lo)
and bound(D_hi), the D at F_true = 1/2 -+ 3.29 sqrt(0.25 / R) (D-131's z, applied to the median's order
statistic); those two ratios are its band.

Quadrature (independent of core/, as ref/ is): the one-period PMF is a trapezoid sum over z in [-12, 12] with step
1/200, in log space (log Phi from log_ndtr, log C(n, d) from gammaln). The integrand is analytic and Gaussian-tailed,
so the trapezoid rule converges geometrically; check_quadrature compares it with scipy.integrate.quad.

The ordering check (reference_ordering.csv): for every pair of independent <= rho_a 0.02 <= 0.12 <= 0.24, every gamma
and every total D up to each (n, T)'s largest 0.9999 quantile of the true total, it lists any D where the order fails,
and each pair's closest approach. It solves about 22,000 bounds, in a process pool (--workers; about 17 min on 23).

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python studies/pluto-tasche/reference.py [--out FILE] [--workers N]
"""

from __future__ import annotations

import argparse
import csv
import math
from functools import cache
from pathlib import Path

import numpy as np
from scipy import integrate
from scipy.optimize import brentq
from scipy.special import gammaln, log_ndtr, logsumexp, ndtri
from scipy.stats import beta, binom

HERE = Path(__file__).resolve().parent

# The recovery grid (tests/recovery/recovery.hpp): id = ((i_pd * 3 + i_rho) * 3 + i_T) * 3 + i_n.
PDS = (0.001, 0.01, 0.05)
RHOS = (0.02, 0.12, 0.24)
PERIODS = (20, 40, 100)
OBLIGORS = (100, 1000, 10000)
GAMMAS = (0.5, 0.75, 0.9, 0.95, 0.99)
RHO_FIXED = 0.12
R = 1000
Z_BAND = 3.29  # D-131's two-sided 0.1% band

Z = np.linspace(-12.0, 12.0, 4801)
H = Z[1] - Z[0]
LOG_W = np.log(np.full(Z.size, H)) - 0.5 * Z * Z - 0.5 * math.log(2.0 * math.pi)
LOG_W[0] -= math.log(2.0)
LOG_W[-1] -= math.log(2.0)

SOLVE_XTOL = 1e-13  # brentq's tolerance in logit p


def scenario(sid: int) -> dict:
    i_n, i_t, i_rho, i_pd = sid % 3, (sid // 3) % 3, (sid // 9) % 3, sid // 27
    return {"scenario": sid, "pd": PDS[i_pd], "rho": RHOS[i_rho], "periods": PERIODS[i_t], "obligors": OBLIGORS[i_n]}


def low_default(s: dict) -> bool:
    """The study's cells: expected defaults per period n PD <= 1 (with a little slack for rounding)."""
    return s["obligors"] * s["pd"] <= 1.0 + 1e-9


SCENARIOS = [s for s in map(scenario, range(81)) if low_default(s)]


def log_cond(pd: float, rho: float, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """log p(z) and log(1 - p(z)), accurate in both tails."""
    x = (ndtri(pd) - math.sqrt(rho) * z) / math.sqrt(1.0 - rho)
    return log_ndtr(x), log_ndtr(-x)


@cache
def period_pmf(n: int, pd: float, rho: float, dmax: int) -> np.ndarray:
    """P(d_t = d) for d = 0..dmax under the one-period mixed binomial."""
    lp, lq = log_cond(pd, rho, Z)
    d = np.arange(dmax + 1, dtype=float)
    lc = gammaln(n + 1.0) - gammaln(d + 1.0) - gammaln(n - d + 1.0)
    terms = lc[None, :] + d[None, :] * lp[:, None] + (n - d)[None, :] * lq[:, None] + LOG_W[:, None]
    return np.exp(logsumexp(terms, axis=0))


def convolve_power(f: np.ndarray, t: int) -> np.ndarray:
    """The PMF of a sum of t iid copies of f, truncated to f's length (exact on 0..len(f) - 1)."""
    m = f.size
    out = np.zeros(m)
    out[0] = 1.0
    base = f.copy()
    while t:
        if t & 1:
            out = np.convolve(out, base)[:m]
        t >>= 1
        if t:
            base = np.convolve(base, base)[:m]
    return out


def total_cdf(n: int, t: int, pd: float, rho: float, dmax: int) -> np.ndarray:
    """F(D) = P(S <= D) for D = 0..dmax."""
    return np.cumsum(convolve_power(period_pmf(n, pd, rho, dmax), t))


def expit(u: float) -> float:
    return 1.0 / (1.0 + math.exp(-u))


def logit(p: float) -> float:
    return math.log(p) - math.log1p(-p)


@cache
def bound_independent(n_total: int, d: int, gamma: float) -> float:
    if d >= n_total:
        return 1.0
    return float(beta.ppf(gamma, d + 1, n_total - d))


@cache
def bound_factor(n: int, t: int, d: int, rho_a: float, gamma: float) -> float:
    """(b): the p with P(S <= d | p, rho_a) = 1 - gamma, solved in logit p."""
    target = math.log1p(-gamma)

    def g(u: float) -> float:
        return math.log(max(total_cdf(n, t, expit(u), rho_a, d)[d], 1e-300)) - target  # underflow: far above the root

    lo, hi = logit(1e-12), logit(0.9)
    assert g(lo) > 0 > g(hi), (n, t, d, rho_a, gamma)
    return expit(brentq(g, lo, hi, xtol=SOLVE_XTOL, rtol=4 * np.finfo(float).eps, maxiter=200))


def dmax_for(s: dict) -> int:
    mean = s["periods"] * s["obligors"] * s["pd"]
    return int(4 * mean + 60)


def first_at_least(cdf: np.ndarray, level: float) -> int:
    idx = int(np.searchsorted(cdf, level, side="left"))
    assert idx < cdf.size, "dmax too small"
    return idx


def methods(s: dict) -> list[tuple[str, float]]:
    return [("independent", 0.0), ("oracle", s["rho"]), ("rho_0.12", RHO_FIXED)]


def bound(method: str, rho_a: float, s: dict, d: int, gamma: float) -> float:
    n, t = s["obligors"], s["periods"]
    if method == "independent":
        return bound_independent(n * t, d, gamma)
    return bound_factor(n, t, d, rho_a, gamma)


def rows() -> list[dict]:
    out = []
    for s in SCENARIOS:
        n, t, pd = s["obligors"], s["periods"], s["pd"]
        dmax = dmax_for(s)
        f_true = total_cdf(n, t, pd, s["rho"], dmax)
        d_med = first_at_least(f_true, 0.5)
        half = Z_BAND * math.sqrt(0.25 / R)
        d_lo, d_hi = first_at_least(f_true, 0.5 - half), first_at_least(f_true, 0.5 + half)
        for method, rho_a in methods(s):
            if method == "independent":
                f_assumed = binom.cdf(np.arange(dmax + 1), n * t, pd)
            else:
                f_assumed = total_cdf(n, t, pd, rho_a, dmax)
            for gamma in GAMMAS:
                d_min = first_at_least(f_assumed, 1.0 - gamma)
                coverage = 1.0 - (f_true[d_min - 1] if d_min > 0 else 0.0)
                b_at = bound(method, rho_a, s, d_min, gamma)
                b_below = bound(method, rho_a, s, d_min - 1, gamma) if d_min > 0 else float("nan")
                assert b_at >= pd and (d_min == 0 or b_below < pd), (s, method, gamma)
                c = coverage
                sd = math.sqrt(max(c * (1.0 - c), 0.0) / R)
                out.append(
                    {
                        **s,
                        "method": method,
                        "rho_a": rho_a,
                        "gamma": gamma,
                        "d_min": d_min,
                        "coverage": coverage,
                        "band_lo": max(0.0, c - Z_BAND * sd),
                        "band_hi": min(1.0, c + Z_BAND * sd),
                        "p_zero": f_true[0],
                        "margin_at_logit": logit(b_at) - logit(pd),
                        "margin_below_logit": (logit(b_below) - logit(pd)) if d_min > 0 else float("nan"),
                        "d_med": d_med,
                        "median_ratio": bound(method, rho_a, s, d_med, gamma) / pd,
                        "d_lo": d_lo,
                        "d_hi": d_hi,
                        "ratio_lo": bound(method, rho_a, s, d_lo, gamma) / pd,
                        "ratio_hi": bound(method, rho_a, s, d_hi, gamma) / pd,
                    }
                )
    return out


def check_quadrature() -> float:
    """The trapezoid PMF against scipy.integrate.quad, at the extremes the study uses. Returns the worst relative error."""
    worst = 0.0
    for n, pd, rho, d in [(100, 0.001, 0.24, 0), (100, 0.001, 0.24, 5), (1000, 0.001, 0.02, 3), (1000, 0.001, 0.24, 12),
                          (100, 0.01, 0.12, 9), (100, 0.01, 0.24, 30), (1000, 0.03, 0.12, 40)]:
        lc = gammaln(n + 1.0) - gammaln(d + 1.0) - gammaln(n - d + 1.0)

        def f(z: float, n: int = n, pd: float = pd, rho: float = rho, d: int = d, lc: float = lc) -> float:
            lp, lq = log_cond(pd, rho, np.array([z]))
            return math.exp(lc + d * lp[0] + (n - d) * lq[0] - 0.5 * z * z) / math.sqrt(2.0 * math.pi)

        exact, _ = integrate.quad(f, -12.0, 12.0, epsabs=0.0, epsrel=1e-13, limit=500, points=[-3, 0, 3])
        trap = period_pmf(n, pd, rho, d)[d]
        worst = max(worst, abs(trap - exact) / exact)
    return worst


def check_closed_forms() -> tuple[float, float]:
    """(a) at D = 0 against 1 - (1 - gamma)^(1/N); (b) at rho_a = 1e-10 against (a). Worst relative errors."""
    worst_zero = worst_limit = 0.0
    for n, t in [(100, 20), (1000, 100)]:
        for gamma in GAMMAS:
            closed = -math.expm1(math.log1p(-gamma) / (n * t))
            worst_zero = max(worst_zero, abs(bound_independent(n * t, 0, gamma) - closed) / closed)
            for d in (0, 1, 7):
                a = bound_independent(n * t, d, gamma)
                b = bound_factor(n, t, d, 1e-10, gamma)
                worst_limit = max(worst_limit, abs(b - a) / a)
    return worst_zero, worst_limit


TIE_LOGIT = 1e-11  # two bounds closer than this in logit are a tie: 100x the root solve's tolerance

# The orderings P6 registers, as (lower, upper) laws: "independent" or a rho_a. Each should give lower <= upper.
PAIRS = [("independent", r) for r in RHOS] + [(0.02, 0.12), (0.12, 0.24), (0.02, 0.24)]


def d_top(n: int, t: int) -> int:
    """The largest 0.9999 quantile of the true total over the scenarios with this (n, T)."""
    return max(quantile_total(tuple(sorted(s.items())), 0.9999) for s in SCENARIOS if (s["obligors"], s["periods"]) == (n, t))


@cache
def quantile_total(s_key: tuple, level: float) -> int:
    """The true total's level quantile, widening the truncation until it holds."""
    s = dict(s_key)
    dmax = dmax_for(s)
    while True:
        cdf = total_cdf(s["obligors"], s["periods"], s["pd"], s["rho"], dmax)
        if cdf[-1] >= level:
            return first_at_least(cdf, level)
        dmax *= 2


def bound_vector(job: tuple[int, int, float, object]) -> tuple[tuple[int, int, float, object], list[float]]:
    """One law's bound at every D = 0..d_top (a process-pool job)."""
    n, t, gamma, law = job
    top = d_top(n, t)
    if law == "independent":
        return job, [bound_independent(n * t, d, gamma) for d in range(top + 1)]
    assert isinstance(law, float)
    return job, [bound_factor(n, t, d, law, gamma) for d in range(top + 1)]


def check_ordering(workers: int) -> tuple[int, int, list[dict], dict, dict]:
    """Per panel, the bound as a function of D, for every pair in PAIRS, at every D from 0 to d_top(n, T) and every
    gamma. Returns (pairs checked, ties, the violations, the bounds by (n, T, gamma, law), and per (pair, gamma) the
    smallest separation logit(upper) - logit(lower) with the (n, T, D) where it occurs)."""
    from concurrent.futures import ProcessPoolExecutor

    cells = sorted({(s["obligors"], s["periods"]) for s in SCENARIOS})
    jobs = [(n, t, g, law) for n, t in cells for g in GAMMAS for law in ["independent", *RHOS]]
    jobs.sort(key=lambda j: -d_top(j[0], j[1]))  # longest first
    with ProcessPoolExecutor(max_workers=workers) as ex:
        bounds = dict(ex.map(bound_vector, jobs))
    checked = ties = 0
    violations = []
    closest: dict = {}
    for n, t in cells:
        for g in GAMMAS:
            for lo_law, hi_law in PAIRS:
                lo_b, hi_b = bounds[(n, t, g, lo_law)], bounds[(n, t, g, hi_law)]
                for d, (lo, hi) in enumerate(zip(lo_b, hi_b)):
                    checked += 1
                    gap = logit(lo) - logit(hi)
                    key = (law_name(lo_law), law_name(hi_law), g)
                    if key not in closest or -gap < closest[key][0]:
                        closest[key] = (-gap, n, t, d)
                    if abs(gap) <= TIE_LOGIT:
                        ties += 1
                    elif gap > 0:
                        violations.append({"obligors": n, "periods": t, "gamma": g, "d": d, "lower": lo_law,
                                           "upper": hi_law, "bound_lower": lo, "bound_upper": hi, "gap_logit": gap})
    return checked, ties, violations, bounds, closest


def law_name(law: object) -> str:
    return "independent" if law == "independent" else f"rho_a={law}"


def scenario_exposure(violations: list[dict]) -> list[dict]:
    """Per scenario and gamma, the true probability of a total D at which a pair that scenario's predictions use is
    violated: (independent, oracle) and (oracle, fixed) or (fixed, oracle)."""
    out = []
    for s in SCENARIOS:
        n, t, rho = s["obligors"], s["periods"], s["rho"]
        pmf = np.diff(total_cdf(n, t, s["pd"], rho, dmax_for(s)), prepend=0.0)
        used: set[tuple[object, object]] = {("independent", rho)}
        if rho != RHO_FIXED:
            used.add((min(rho, RHO_FIXED), max(rho, RHO_FIXED)))
        for g in GAMMAS:
            ds = sorted({v["d"] for v in violations
                         if (v["obligors"], v["periods"], v["gamma"]) == (n, t, g) and (v["lower"], v["upper"]) in used})
            if ds:
                out.append({"scenario": s["scenario"], "gamma": g, "d": " ".join(map(str, ds)),
                            "probability": float(sum(pmf[d] for d in ds if d < pmf.size))})
    return out


COLUMNS = [
    "scenario", "pd", "rho", "periods", "obligors", "method", "rho_a", "gamma", "d_min", "coverage", "band_lo",
    "band_hi", "p_zero", "margin_at_logit", "margin_below_logit", "d_med", "median_ratio", "d_lo", "d_hi", "ratio_lo",
    "ratio_hi",
]


def main() -> None:
    ap = argparse.ArgumentParser(description="S-6 exact reference numbers (reference.csv)")
    ap.add_argument("--out", default=str(HERE / "reference.csv"))
    ap.add_argument("--ordering-out", default=str(HERE / "reference_ordering.csv"))
    ap.add_argument("--workers", type=int, default=None, help="processes for the ordering check (default: all)")
    args = ap.parse_args()

    q = check_quadrature()
    zero, limit = check_closed_forms()
    print(f"quadrature: worst relative error against quad {q:.2e}")
    print(f"(a) at D = 0 against the closed form: {zero:.2e}; (b) at rho_a = 1e-10 against (a): {limit:.2e}")
    assert q < 1e-9 and zero < 1e-12 and limit < 1e-6

    checked, ties, violations, _, closest = check_ordering(args.workers)
    exposure = scenario_exposure(violations)
    summary = f"per-D ordering: {len(violations)} violations and {ties} ties (|gap| <= {TIE_LOGIT:g} logit) in {checked} pairs"
    print(summary)
    with open(args.ordering_out, "w", newline="", encoding="utf-8") as fh:
        fh.write("# S-6 per-D ordering check (studies/pluto-tasche/reference.py): every (n, T, gamma, D, pair) at which the\n"
                 "# registered ordering lower <= upper fails, D = 0 to each (n, T)'s largest 0.9999 quantile of the true total.\n"
                 f"# {summary}\n")
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["obligors", "periods", "gamma", "d", "lower", "upper", "bound_lower", "bound_upper", "gap_logit"])
        for v in violations:
            w.writerow([v["obligors"], v["periods"], v["gamma"], v["d"], law_name(v["lower"]), law_name(v["upper"]),
                        f"{v['bound_lower']:.12g}", f"{v['bound_upper']:.12g}", f"{v['gap_logit']:.4g}"])
        fh.write("# Exposure: per scenario and gamma, the true probability of a total at which a pair its predictions use\n"
                 "# (independent vs oracle; oracle vs the fixed arm) is violated.\n")
        w.writerow(["scenario", "gamma", "d", "probability"])
        for e in exposure:
            w.writerow([e["scenario"], e["gamma"], e["d"], f"{e['probability']:.6g}"])
        fh.write("# Closest approach: per pair and gamma, the smallest logit(upper) - logit(lower) over every (n, T, D)\n"
                 "# checked, and where it occurs (negative would be a violation).\n")
        w.writerow(["lower", "upper", "gamma", "min_gap_logit", "obligors", "periods", "d"])
        for (lo_name, hi_name, g), (gap, n, t, d) in sorted(closest.items()):
            w.writerow([lo_name, hi_name, g, f"{gap:.4g}", n, t, d])
    print(f"wrote {len(violations)} violations, {len(exposure)} exposed (scenario, gamma) cells to {args.ordering_out}")

    table = rows()
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        fh.write("# S-6 exact reference (studies/pluto-tasche/reference.py): one row per scenario x method x gamma.\n")
        fh.write(f"# quadrature check {q:.2e}; zero-default closed form {zero:.2e}; rho_a -> 0 limit {limit:.2e}; "
                 f"{summary}\n")
        w = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        for r in table:
            w.writerow({k: (f"{v:.10g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {len(table)} rows to {args.out}")


if __name__ == "__main__":
    main()
