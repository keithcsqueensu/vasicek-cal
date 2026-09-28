# SPDX-License-Identifier: Apache-2.0
"""Replicate the Vasicek-rate MLE with plain scipy (M3; docs/methodology/vasicek_rate_mle.md).

The log-likelihood of an observed default rate r, with x = Phi^-1(r) and c = Phi^-1(PD), is
    l(PD, rho) = 1/2 log(1 - rho) - 1/2 log(rho) + x^2/2 - (sqrt(1 - rho) x - c)^2 / (2 rho);
under the censored treatment a rate of 0 contributes log Phi((sqrt(1 - rho) Phi^-1(delta) - c) / sqrt(rho)),
delta the detection limit (count data: 1/(2n)).

For each panel of tests/golden/vasicek_rate/replay.csv (rate panels from the recovery DGP, and count
panels with zero-default periods under the censored likelihood) this script, on its own:
  - finds the maximum over the box (PD in [1e-4, 0.2], rho in [1e-3, 0.5]) by nested bounded Brent in
    logit coordinates, and, for rate panels, also the closed form (the mean and 1/T variance of x);
  - solves each profile-interval end point the engine solved, by brentq on the profile near it (the
    profile's inner maximum by bounded Brent), and checks that an end the engine truncated at the
    box is inside the interval there.
The 99.9% conditional PD's interval uses S-23's method, cross-checked by conditional_pd_profile.py.

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/vasicek_rate_mle.py
"""

import csv
import math
import sys
from pathlib import Path

import numpy as np
import tomllib
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[2]
PD_BOX, RHO_BOX = (1e-4, 0.2), (1e-3, 0.5)
THRESHOLD = stats.chi2.ppf(0.95, 1) / 2.0
LOWER_TRUNC, UPPER_TRUNC, NOT_COMPUTED = 1, 2, 4  # engine interval flag bits
XTOL = 1e-11


def logit(v):
    return math.log(v) - math.log1p(-v)


def inv_logit(u):
    return 1.0 / (1.0 + math.exp(-u))


class Panel:
    def __init__(self, rates, detect, censored):
        r = np.array(rates)
        self.x = stats.norm.ppf(r[(r > 0) & (r < 1)])
        self.zeros = int(np.sum(r == 0))
        self.ones = int(np.sum(r == 1))
        self.xd = stats.norm.ppf(detect)
        assert censored or (self.zeros == 0 and self.ones == 0)

    def loglik(self, pd, rho):
        c = stats.norm.ppf(pd)
        u = math.sqrt(1 - rho) * self.x - c
        ll = np.sum(0.5 * math.log1p(-rho) - 0.5 * math.log(rho) + 0.5 * self.x**2 - u**2 / (2 * rho))
        s1, s = math.sqrt(1 - rho), math.sqrt(rho)
        if self.zeros:
            ll += self.zeros * stats.norm.logcdf((s1 * self.xd - c) / s)
        if self.ones:
            ll += self.ones * stats.norm.logcdf((s1 * self.xd + c) / s)
        return float(ll)

    def profile(self, axis, v):
        """The maximum over the other parameter, in its box, with parameter `axis` fixed at v."""
        lo, hi = (RHO_BOX if axis == 0 else PD_BOX)
        f = (lambda u: -self.loglik(v, inv_logit(u))) if axis == 0 else (lambda u: -self.loglik(inv_logit(u), v))
        res = optimize.minimize_scalar(f, bounds=(logit(lo), logit(hi)), method="bounded", options={"xatol": XTOL})
        best = -res.fun
        for end in (logit(lo), logit(hi)):  # bounded Brent never evaluates the ends themselves
            best = max(best, -f(end))
        return best

    def maximum(self):
        g = lambda u: -self.profile(0, inv_logit(u))
        res = optimize.minimize_scalar(g, bounds=(logit(PD_BOX[0]), logit(PD_BOX[1])), method="bounded",
                                       options={"xatol": XTOL})
        pd = inv_logit(res.x)
        h = lambda u: -self.loglik(pd, inv_logit(u))
        rr = optimize.minimize_scalar(h, bounds=(logit(RHO_BOX[0]), logit(RHO_BOX[1])), method="bounded",
                                      options={"xatol": XTOL})
        return pd, inv_logit(rr.x), -res.fun

    def closed_form(self):
        mu, v = float(np.mean(self.x)), float(np.var(self.x))
        return stats.norm.cdf(mu / math.sqrt(1 + v)), v / (1 + v)


def main():
    with open(ROOT / "tests" / "tolerances.toml", "rb") as f:
        tol = {k: v["value"] for k, v in tomllib.load(f).items()}
    path = ROOT / "tests" / "golden" / "vasicek_rate" / "replay.csv"
    with open(path, newline="") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    h = float.fromhex
    worst = {"estimate_u": 0.0, "loglik": 0.0, "end_u": 0.0, "closed_form_u": 0.0}
    failures = []
    for row in rows:
        rates = [h(v) for v in row["rates_hex"].split()]
        p = Panel(rates, h(row["detect_hex"]), row["treatment"] == "censor")
        at = f"{row['kind']} {row['source_scenario']}/{row['replicate']}"
        pd_e, rho_e, ll_e = h(row["pd_hex"]), h(row["rho_hex"]), h(row["loglik_max_hex"])
        pd_s, rho_s, ll_s = p.maximum()
        du = max(abs(logit(pd_s) - logit(pd_e)), abs(logit(rho_s) - logit(rho_e)))
        worst["estimate_u"] = max(worst["estimate_u"], du)
        worst["loglik"] = max(worst["loglik"], abs(ll_s - ll_e))
        if du > tol["TOL_SCIPY_RATE_ESTIMATE_U"] or abs(ll_s - ll_e) > tol["TOL_SCIPY_RATE_LOGLIK_ABS"]:
            failures.append(f"{at}: maximum differs (logit {du:.3g}, loglik {abs(ll_s - ll_e):.3g})")
        if row["kind"] == "rates":
            pd_c, rho_c = p.closed_form()
            if PD_BOX[0] < pd_c < PD_BOX[1] and RHO_BOX[0] < rho_c < RHO_BOX[1]:
                dc = max(abs(logit(pd_c) - logit(pd_e)), abs(logit(rho_c) - logit(rho_e)))
                worst["closed_form_u"] = max(worst["closed_form_u"], dc)
                if dc > tol["TOL_SCIPY_RATE_ESTIMATE_U"]:
                    failures.append(f"{at}: closed form differs by {dc:.3g} in logit")
        flags = int(row["interval_flags"])
        target = ll_s - THRESHOLD
        for axis, (lo_col, hi_col) in enumerate((("pd_lo_hex", "pd_hi_hex"), ("rho_lo_hex", "rho_hi_hex"))):
            f_axis = (flags >> (8 * axis)) & 0xFF
            if f_axis & NOT_COMPUTED:
                failures.append(f"{at}: the engine did not compute axis {axis}'s interval")
                continue
            box = PD_BOX if axis == 0 else RHO_BOX
            for end, col, trunc in ((0, lo_col, LOWER_TRUNC), (1, hi_col, UPPER_TRUNC)):
                v_e = h(row[col])
                if f_axis & trunc:
                    # Truncated at the box: the bound itself must lie inside the interval here too.
                    if p.profile(axis, box[end]) < target:
                        failures.append(f"{at}: axis {axis} end {end} truncated by the engine, outside here")
                    continue
                g = lambda u: p.profile(axis, inv_logit(u)) - target  # noqa: B023
                u0, w = logit(v_e), 0.02
                while g(u0 - w) * g(u0 + w) > 0 and w < 1.0:
                    w *= 2
                u_s = optimize.brentq(g, u0 - w, u0 + w, xtol=1e-12)
                worst["end_u"] = max(worst["end_u"], abs(u_s - u0))
                if abs(u_s - u0) > tol["TOL_SCIPY_RATE_PROFILE_U"]:
                    failures.append(f"{at}: axis {axis} end {end} differs by {abs(u_s - u0):.3g} in logit")
    print(f"{len(rows)} panels ({sum(r['kind'] == 'counts' for r in rows)} censored count panels)")
    for k, v in worst.items():
        print(f"  worst {k}: {v:.3g}")
    for line in failures:
        print("FAIL " + line, file=sys.stderr)
    print("scipy agrees with the engine" if not failures else f"{len(failures)} disagreement(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
