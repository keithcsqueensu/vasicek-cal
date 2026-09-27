# SPDX-License-Identifier: Apache-2.0
"""Cross-check the profile-likelihood interval for the 99.9% conditional PD with plain scipy (S-23).

The quantity (studies/derived-quantity-intervals/PREDICTION.md) is
    q(PD, rho) = Phi((Phi^-1(PD) + sqrt(rho) z) / sqrt(1 - rho)),   z = Phi^-1(0.999),
the PD when the systematic factor is at its 0.1% adverse quantile. On the curve q = c,
    PD_c(rho) = Phi(sqrt(1 - rho) Phi^-1(c) - sqrt(rho) z),
so the profile P_q(c) is a maximisation over rho alone, over the rho for which (PD_c(rho), rho) lies
in the estimation box. The 95% interval is {c : P_q(c) >= l_max - 1.9207...}.

For each recovery replay panel (tests/golden/recovery/replay_panels.csv) this script computes, on
its own:
  - the panel log-likelihood: each period's integral over the factor by the trapezoid rule on a
    window around the integrand's mode, widened until both ends are negligible;
  - l_max, by nested bounded Brent in logit coordinates over the box;
  - each end point the engine solved (tests/golden/recovery/replay.csv), by brentq on logit(c)
    near it, the inner maximum by bounded Brent over the feasible stretch of logit(rho).
It checks the end points agree within TOL_SCIPY_Q_PROFILE_ENDPOINT_S (logit units), that an end
the engine truncated at the limit of q in the box is inside the interval here too, and that the
engine's box-limited flags match where this script's inner maximiser lies.

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/conditional_pd_profile.py [DIR]

DIR holds replay.csv and replay_panels.csv (default: tests/golden/recovery).
"""

import csv
import math
import sys
import tomllib
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import optimize, special, stats

ROOT = Path(__file__).resolve().parents[2]
PD_BOX, RHO_BOX = (1e-4, 0.2), (1e-3, 0.5)  # the recovery grid's box (D-115)
Z = stats.norm.ppf(0.999)
THRESHOLD = stats.chi2.ppf(0.95, 1) / 2.0
LOWER_BOX, UPPER_BOX, NOT_COMPUTED, LOWER_TRUNC, UPPER_TRUNC = 8, 16, 4, 1, 2  # engine flag bits


def logit(v):
    return math.log(v) - math.log1p(-v)


def inv_logit(u):
    return 1.0 / (1.0 + math.exp(-u))


class Panel:
    """Distinct default counts and their multiplicities; the log-likelihood at (PD, rho)."""

    def __init__(self, n, ds):
        counts = Counter(ds)
        self.n = n
        self.d = np.array(sorted(counts), dtype=float)
        self.k = np.array([counts[int(d)] for d in self.d], dtype=float)
        self.lchoose = special.gammaln(n + 1) - special.gammaln(self.d + 1) - special.gammaln(n - self.d + 1)

    def loglik(self, pd, rho):
        c, a, b = stats.norm.ppf(pd), math.sqrt(rho), math.sqrt(1.0 - rho)
        d, n = self.d[:, None], self.n

        def h(z):  # log(binom.pmf(d; n, p(z)) phi(z)), p(z) = Phi((c - a z) / b)
            x = (c - a * z) / b
            return (self.lchoose[:, None] + d * special.log_ndtr(x) + (n - d) * special.log_ndtr(-x)
                    - 0.5 * z * z - 0.5 * math.log(2.0 * math.pi))

        # Gaussian approximation to each integrand: the binomial likelihood in z, times phi(z).
        p0 = (self.d + 0.5) / (n + 1.0)
        x0 = stats.norm.ppf(p0)
        z0 = (c - b * x0) / a
        s_lik = np.sqrt(p0 * (1.0 - p0) / n) / (a / b * stats.norm.pdf(x0))
        prec = 1.0 / s_lik**2 + 1.0
        mode, sd = z0 / s_lik**2 / prec, 1.0 / np.sqrt(prec)
        lo, hi = mode - 14.0 * sd, mode + 14.0 * sd
        for _ in range(40):  # widen each window until its ends are e^-60 below its maximum
            grid = lo[:, None] + (hi - lo)[:, None] * np.linspace(0.0, 1.0, 2001)[None, :]
            v = h(grid)
            top = v.max(axis=1)
            low_ok, high_ok = v[:, 0] < top - 60.0, v[:, -1] < top - 60.0
            if low_ok.all() and high_ok.all():
                break
            width = hi - lo
            lo = np.where(low_ok, lo, lo - width)
            hi = np.where(high_ok, hi, hi + width)
        else:
            raise RuntimeError("integration window did not converge")
        step = (hi - lo) / 2000.0
        w = np.exp(v - top[:, None])
        area = (w.sum(axis=1) - 0.5 * (w[:, 0] + w[:, -1])) * step
        return float(np.sum(self.k * (top + np.log(area))))


def minimise_bounded(f, lo, hi):
    """Bounded Brent on [lo, hi]. It stops about sqrt(eps)|x| short of a bound, so a minimum within
    1e-6 of a bound is also tried on the bound itself. Returns (f, x, on_bound)."""
    r = optimize.minimize_scalar(f, bounds=(lo, hi), method="bounded", options={"xatol": 1e-10})
    fun, x = r.fun, r.x
    on_bound = x - lo < 1e-6 or hi - x < 1e-6
    if on_bound:
        end = lo if x - lo < 1e-6 else hi
        f_end = f(end)
        if f_end < fun:
            fun, x = f_end, end
    return fun, x, on_bound


def pd_on_curve(x, w):
    """PD on the curve Phi^-1(q) = x at logit(rho) = w (array), NaN outside the PD box."""
    rho = 1.0 / (1.0 + np.exp(-w))
    pd = special.ndtr(np.sqrt(1.0 - rho) * x - np.sqrt(rho) * Z)
    return np.where((pd >= PD_BOX[0]) & (pd <= PD_BOX[1]), pd, np.nan)


class Profile:
    """P_q at logit(c) = s, remembering where the last inner maximum was."""

    W_LO, W_HI = logit(RHO_BOX[0]), logit(RHO_BOX[1])

    def __init__(self, panel):
        self.panel, self.last = panel, None

    def __call__(self, s):
        x = stats.norm.ppf(inv_logit(s)) if s < 0 else -stats.norm.ppf(inv_logit(-s))
        scan = np.linspace(self.W_LO, self.W_HI, 4001)
        ok = ~np.isnan(pd_on_curve(x, scan))
        if not ok.any():
            return -math.inf, False
        # The feasible stretches of logit(rho): found on the scan, their inner ends refined by bisection.
        edges = np.flatnonzero(np.diff(ok.astype(int)))
        starts = [0, *(edges + 1)]
        ends = [*edges, len(scan) - 1]
        step = scan[1] - scan[0]
        spans = []
        for i, j in zip(starts, ends):
            if ok[i]:
                lo = scan[i] if i == 0 else self._bisect(x, scan[i], scan[i] - step)
                hi = scan[j] if j == len(scan) - 1 else self._bisect(x, scan[j], scan[j] + step)
                spans.append((lo, hi))

        def neg(w):
            return -self.panel.loglik(float(pd_on_curve(x, np.array([w]))[0]), inv_logit(w))

        # Bounded Brent over each stretch; near the last maximiser first, which is cheaper when it is
        # interior to the narrower bracket.
        results = []
        for lo, hi in spans:
            if self.last is not None and lo <= self.last <= hi:
                a, b = max(lo, self.last - 0.3), min(hi, self.last + 0.3)
                fun, w, _ = minimise_bounded(neg, a, b)
                if not ((w - a < 1e-6 < a - lo) or (b - w < 1e-6 < hi - b)):
                    results.append((fun, w, (w - lo < 1e-6) or (hi - w < 1e-6)))
                    continue
            results.append(minimise_bounded(neg, lo, hi))
        fun, w, on_bound = min(results)
        self.last = w
        return -fun, on_bound

    @staticmethod
    def _bisect(x, inside, outside):
        """The last feasible logit(rho) between a feasible and an infeasible point."""
        for _ in range(100):
            mid = 0.5 * (inside + outside)
            if np.isnan(pd_on_curve(x, np.array([mid])))[0]:
                outside = mid
            else:
                inside = mid
        return inside


def l_max(panel):
    """Nested bounded Brent over the whole box: logit PD outside, logit rho inside. (Nelder-Mead
    can stall on a bound of the box when the maximum lies just inside it.)"""
    w_box, u_box = (logit(RHO_BOX[0]), logit(RHO_BOX[1])), (logit(PD_BOX[0]), logit(PD_BOX[1]))

    def neg_profile(u):
        return minimise_bounded(lambda w: -panel.loglik(inv_logit(u), inv_logit(w)), *w_box)[0]

    return -minimise_bounded(neg_profile, *u_box)[0]


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def main():
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "tests" / "golden" / "recovery"
    with open(ROOT / "tests" / "tolerances.toml", "rb") as f:
        tol = {k: v["value"] for k, v in tomllib.load(f).items()}["TOL_SCIPY_Q_PROFILE_ENDPOINT_S"]
    panels = {(int(r["scenario"]), int(r["replicate"])): Panel(int(r["n"]), [int(d) for d in r["d"].split()])
              for r in read_csv(folder / "replay_panels.csv")}
    worst, where, bad, solved, truncated = 0.0, "", [], 0, 0
    for r in read_csv(folder / "replay.csv"):
        key = (int(r["scenario"]), int(r["replicate"]))
        flags = int(r["q_profile_flags"])
        if flags & NOT_COMPUTED:
            bad.append(f"{key}: not computed by the engine")
            continue
        panel = panels[key]
        level = l_max(panel) - THRESHOLD
        prof = Profile(panel)
        for side, col, trunc_bit, box_bit in (("lower", "q_lo_hex", LOWER_TRUNC, LOWER_BOX),
                                              ("upper", "q_hi_hex", UPPER_TRUNC, UPPER_BOX)):
            s_e = logit(float.fromhex(r[col]))
            if flags & trunc_bit:  # the limit of q in the box: it must be inside the interval here too
                truncated += 1
                inside = prof(s_e)[0] >= level
                if not inside:
                    bad.append(f"{key} {side}: truncated by the engine, but outside the interval here")
                continue
            prof.last = None
            g = lambda s: prof(s)[0] - level  # noqa: E731
            delta = 1e-3
            while delta < 1.0 and g(s_e - delta) * g(s_e + delta) > 0.0:
                delta *= 4.0
            if delta >= 1.0:
                bad.append(f"{key} {side}: no crossing within 1 logit unit of the engine's end point")
                continue
            s_scipy = optimize.brentq(g, s_e - delta, s_e + delta, xtol=1e-12, rtol=4 * np.finfo(float).eps)
            solved += 1
            diff = abs(s_scipy - s_e)
            if diff > worst:
                worst, where = diff, f"scenario {key[0]} replicate {key[1]}, {side} end"
            on_bound = prof(s_scipy)[1]
            if on_bound != bool(flags & box_bit):
                bad.append(f"{key} {side}: box-limited here {on_bound}, engine flag {bool(flags & box_bit)}")
    passed = worst <= tol and not bad
    for b in bad:
        print("  " + b)
    print(f"{len(panels)} panels, {solved} solved end points, {truncated} truncated: worst |logit(c) scipy - "
          f"engine| {worst:.3g} ({where}); tolerance {tol:g}; {'ok' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
