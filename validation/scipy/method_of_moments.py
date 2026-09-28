# SPDX-License-Identifier: Apache-2.0
"""Replicate the method of moments with plain scipy (M3; engine/moments.hpp).

Joint-default-probability form: PD-hat = sum d / sum n, PD_2-hat = sum d(d - 1) / sum n(n - 1) for
count panels (for rate series, the mean rate and the mean squared rate), and rho-hat solves
    Phi_2(c, c; rho) = PD_2-hat,   c = Phi^-1(PD-hat),
within [1e-3, 0.5]; at or beyond an end, rho-hat is that end (flagged). Phi_2 at equal arguments is
computed exactly with Owen's T function, Phi_2(h, h; rho) = Phi(h) - 2 T(h, sqrt((1 - rho) / (1 + rho))),
independent of the engine's quadrature, and the root by brentq in logit(rho).

Checks every row of tests/golden/moments/reference.csv: PD-hat and PD_2-hat to
TOL_SCIPY_MOM_MOMENT_REL, rho-hat to TOL_SCIPY_MOM_RHO_U in logit units, and the flags.

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/method_of_moments.py
"""

import csv
import math
import sys
from pathlib import Path

import tomllib
from scipy import optimize, special, stats

ROOT = Path(__file__).resolve().parents[2]
PD_BOX, RHO_BOX = (1e-4, 0.2), (1e-3, 0.5)
REFUSED, FLOOR, CAP, PD_OUTSIDE = 1, 2, 4, 8  # engine/moments.hpp flag bits


def phi2_equal(h, rho):
    return stats.norm.cdf(h) - 2.0 * special.owens_t(h, math.sqrt((1.0 - rho) / (1.0 + rho)))


def logit(v):
    return math.log(v) - math.log1p(-v)


def invert(pd, pd2):
    """rho-hat and flags from PD-hat and PD_2-hat."""
    if not pd > 0:
        return float("nan"), REFUSED
    flags = PD_OUTSIDE if not PD_BOX[0] <= pd <= PD_BOX[1] else 0
    h = stats.norm.ppf(pd)
    target = math.log(pd2) if pd2 > 0 else -math.inf

    def g(u):
        return math.log(phi2_equal(h, 1.0 / (1.0 + math.exp(-u)))) - target

    a, b = logit(RHO_BOX[0]), logit(RHO_BOX[1])
    if not g(a) < 0:
        return RHO_BOX[0], flags | FLOOR
    if not g(b) > 0:
        return RHO_BOX[1], flags | CAP
    u = optimize.brentq(g, a, b, xtol=1e-13, rtol=4 * sys.float_info.epsilon)
    return 1.0 / (1.0 + math.exp(-u)), flags


def main():
    with open(ROOT / "tests" / "tolerances.toml", "rb") as f:
        tol = {k: v["value"] for k, v in tomllib.load(f).items()}
    path = ROOT / "tests" / "golden" / "moments" / "reference.csv"
    with open(path, newline="") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    h = float.fromhex
    worst_m = worst_u = 0.0
    failures = []
    for r in rows:
        at = f"{r['kind']} {r['cell_scenario']}/{r['replicate']}"
        if r["kind"] == "counts":
            n = int(r["n"])
            d = [int(x) for x in r["data"].split()]
            pd = sum(d) / (n * len(d))
            pd2 = sum(x * (x - 1) for x in d) / (len(d) * n * (n - 1))
        else:
            rates = [h(x) for x in r["data"].split()]
            pd = sum(rates) / len(rates)
            pd2 = sum(x * x for x in rates) / len(rates)
        rho, flags = invert(pd, pd2)
        for got, col in ((pd, "pd_hat_hex"), (pd2, "pd2_hat_hex")):
            want = h(r[col])
            rel = abs(got / want - 1.0)
            worst_m = max(worst_m, rel)
            if rel > tol["TOL_SCIPY_MOM_MOMENT_REL"]:
                failures.append(f"{at}: {col} differs by {rel:.3g} relative")
        if flags != int(r["flags"]):
            failures.append(f"{at}: flags {flags}, engine {r['flags']}")
        elif not flags & REFUSED:
            du = abs(logit(rho) - logit(h(r["rho_hat_hex"])))
            worst_u = max(worst_u, du)
            if du > tol["TOL_SCIPY_MOM_RHO_U"]:
                failures.append(f"{at}: rho-hat differs by {du:.3g} in logit")
    print(f"{len(rows)} panels; worst moment difference {worst_m:.3g} relative, worst rho-hat {worst_u:.3g} in logit")
    for line in failures:
        print("FAIL " + line, file=sys.stderr)
    print("scipy agrees with the engine" if not failures else f"{len(failures)} disagreement(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
