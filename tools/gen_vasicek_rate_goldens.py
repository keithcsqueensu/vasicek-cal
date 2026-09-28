# SPDX-License-Identifier: Apache-2.0
"""Golden values for the Vasicek-rate objective (core/objectives/vasicek_rate.hpp; M3, D-033).

Every value is computed in mpmath at 50 significant digits and rounded once to binary64:

- log_density.csv: the log density of the Vasicek rate distribution,
      l = 1/2 log(1 - rho) - 1/2 log(rho) + x^2/2 - (sqrt(1 - rho) x - c)^2 / (2 rho),
  x = Phi^-1(r), c = Phi^-1(PD), over rates from 1e-8 to 1 - 1e-6 and the box's extremes of PD and rho;
- boundary.csv: the zero-rate treatments at detection limits 1/(2n), n = 100 and 10^6: the censored
  terms log P(R <= detect) and log P(R >= 1 - detect), and the substituted densities at detect and
  1 - detect.

Each row carries the size of the terms the double computation adds (scale_hex), for the tolerance
TOL_VASICEK_RATE_LOGLIK_EPS.

    uv run --no-project --with mpmath==1.4.1 python tools/gen_vasicek_rate_goldens.py
    uv run --no-project --with mpmath==1.4.1 python tools/gen_vasicek_rate_goldens.py --check
"""

from __future__ import annotations

import math
import sys

from mpmath import mp, mpf
from tablegen import csv_note, log_phi_mp, probit_mp, row, run, to_double

GENERATOR = "tools/gen_vasicek_rate_goldens.py"
MANIFEST = "tests/golden/vasicek_rate/MANIFEST.json"
OUT = "tests/golden/vasicek_rate"

PDS = ["1e-4", "0.001", "0.01", "0.05", "0.2"]
RHOS = ["1e-3", "0.02", "0.12", "0.24", "0.5"]
RATES = ["1e-8", "1e-5", "0.001", "0.01", "0.05", "0.2", "0.5", "0.8", "0.99", "0.999999"]
NS = [100, 1000000]


def density(x: mpf, c: mpf, rho: mpf) -> mpf:
    u = mp.sqrt(1 - rho) * x - c
    return mp.log(1 - rho) / 2 - mp.log(rho) / 2 + x * x / 2 - u * u / (2 * rho)


def density_scale(x: float, c: float, rho: float) -> float:
    # The size of the terms added in double: |1/2 log(1 - rho)|, |1/2 log rho|, x^2/2 and the square
    # of u's terms (u can cancel, so its error is relative to |sqrt(1 - rho) x| + |c|) over 2 rho.

    return (abs(0.5 * math.log1p(-rho)) + abs(0.5 * math.log(rho)) + 0.5 * x * x
            + (abs(math.sqrt(1 - rho) * x) + abs(c)) ** 2 / (2 * rho))


def d(s: str) -> float:
    """The binary64 input the C++ test uses: goldens are evaluated at it exactly."""
    return float(s)


def build_outputs() -> dict[str, str]:
    dens = [csv_note(GENERATOR),
            "pd_hex,rho_hex,rate_hex,log_density_hex,scale_hex,pd,rho,rate,log_density\n"]
    for pd_s in PDS:
        for rho_s in RHOS:
            for r_s in RATES:
                pd, rho, r = mpf(d(pd_s)), mpf(d(rho_s)), mpf(d(r_s))
                c, x = probit_mp(pd), probit_mp(r)
                v = density(x, c, rho)
                scale = density_scale(to_double(x), to_double(c), to_double(rho))
                dens.append(row(d(pd_s).hex(), d(rho_s).hex(), d(r_s).hex(), to_double(v).hex(), scale.hex(),
                                pd_s, rho_s, r_s, mp.nstr(v, 20)))
    bound = [csv_note(GENERATOR),
             "pd_hex,rho_hex,detect_hex,treatment,rate,value_hex,scale_hex,pd,rho,n,value\n"]
    for pd_s in PDS:
        for rho_s in RHOS:
            for n in NS:
                pd, rho = mpf(d(pd_s)), mpf(d(rho_s))
                det_d = 0.5 / n  # rate_obs's detect, in double
                det = mpf(det_d)
                c, xd = probit_mp(pd), probit_mp(det)
                s1, s = mp.sqrt(1 - rho), mp.sqrt(rho)
                wz, wo = (s1 * xd - c) / s, (s1 * xd + c) / s
                cases = [("censor", 0, log_phi_mp(wz), wz), ("censor", 1, log_phi_mp(wo), wo),
                         ("substitute", 0, density(xd, c, rho), None), ("substitute", 1, density(-xd, c, rho), None)]
                for name, rate, v, w in cases:
                    if w is None:
                        scale = density_scale(to_double(xd), to_double(c), to_double(rho))
                    else:
                        # log Phi(w): w's terms, times the slope of log Phi (about max(1, |w|)).
                        wd = to_double(w)
                        scale = (abs(to_double(s1 * xd)) + abs(to_double(c))) / to_double(s) * max(1.0, abs(wd))
                        scale = max(scale, abs(to_double(v)))
                    bound.append(row(d(pd_s).hex(), d(rho_s).hex(), det_d.hex(), name, rate, to_double(v).hex(),
                                     scale.hex(), pd_s, rho_s, n, mp.nstr(v, 20)))
    return {f"{OUT}/log_density.csv": "".join(dens), f"{OUT}/boundary.csv": "".join(bound)}


if __name__ == "__main__":
    sys.exit(run(GENERATOR, MANIFEST, build_outputs, __doc__))
