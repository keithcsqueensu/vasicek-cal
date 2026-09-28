# SPDX-License-Identifier: Apache-2.0
"""Golden values for the method of moments' joint default probability (engine/moments.hpp; M3, D-033).

    PD_2(PD, rho) = Phi_2(c, c; rho) = integral of Phi((c - sqrt(rho) z) / sqrt(1 - rho))^2 phi(z) dz,

c = Phi^-1(PD), by mpmath's tanh-sinh quadrature at 50 significant digits (split at z = 0 and at the
integrand's peak), independent of the engine's Gauss-Hermite rule. Evaluated at the binary64 inputs
and rounded once. Output: tests/golden/moments/joint_default.csv with log PD_2.

    uv run --no-project --with mpmath==1.4.1 python tools/gen_moments_goldens.py
    uv run --no-project --with mpmath==1.4.1 python tools/gen_moments_goldens.py --check
"""

from __future__ import annotations

import sys

from mpmath import mp, mpf
from tablegen import csv_note, probit_mp, row, run, to_double

GENERATOR = "tools/gen_moments_goldens.py"
MANIFEST = "tests/golden/moments/MANIFEST.json"
OUT = "tests/golden/moments"

PDS = ["1e-4", "0.001", "0.01", "0.05", "0.2"]
RHOS = ["1e-3", "0.02", "0.12", "0.24", "0.5"]


def joint_default(pd: mpf, rho: mpf) -> mpf:
    c = probit_mp(pd)
    s, s1 = mp.sqrt(rho), mp.sqrt(1 - rho)

    def f(z):
        return mp.ncdf((c - s * z) / s1) ** 2 * mp.npdf(z)

    # The integrand is a normal density times Phi^2, which is ~1 for z below c / s and falls away
    # above it: split at 0 and at the knee c / s so tanh-sinh sees smooth pieces.
    knee = c / s
    points = sorted({mpf(0), knee, knee + 8 * s1 / s, mpf(-40), mpf(40)})
    points = [p for p in points if -40 <= p <= 40]
    return mp.quad(f, points)


def build_outputs() -> dict[str, str]:
    out = [csv_note(GENERATOR), "pd_hex,rho_hex,log_pd2_hex,pd,rho,log_pd2\n"]
    for pd_s in PDS:
        for rho_s in RHOS:
            pd, rho = mpf(float(pd_s)), mpf(float(rho_s))
            v = mp.log(joint_default(pd, rho))
            out.append(row(float(pd_s).hex(), float(rho_s).hex(), to_double(v).hex(), pd_s, rho_s, mp.nstr(v, 20)))
    return {f"{OUT}/joint_default.csv": "".join(out)}


if __name__ == "__main__":
    sys.exit(run(GENERATOR, MANIFEST, build_outputs, __doc__))
