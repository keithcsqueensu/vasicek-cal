# SPDX-License-Identifier: Apache-2.0
"""S-34, the mechanism behind its predictions: the large-n limit, from the model alone.

No recovery panel is read. As n -> infinity each period's probit default rate is exact:
x_t = (c - sqrt(rho) z_t) / sqrt(1 - rho), c = probit(PD), z_t ~ N(0, 1) (the engine's convention,
higher Z = better conditions). So x_t ~ N(m, s^2) with m = c / sqrt(1 - rho), s^2 = rho / (1 - rho),
and the MLE is the sample mean mu and the 1/T variance v of the x_t:
    rho-hat = v / (1 + v),  probit(PD-hat) = mu / sqrt(1 + v),  probit(q-hat) = mu + z999 sqrt(v).
A severe period at z = -u (u = probit(0.99) or probit(0.999)) adds x = m + s u. Standard errors are
the original panel's, by the delta method from var(mu) = v / T and var(v) = 2 v^2 / T.

Prints, per T, severity and number of added periods k (medians over draws of T-period histories):
  shift of q-hat and rho-hat in SE units (q's does not depend on PD or rho; rho's barely does);
  P(new q-hat above the original interval's upper end);
  coverage of the unchanged true q by the Wald interval, before and after (T + k periods).
Also the ratio of q-hat's shift for two periods to one, and rho-hat's shift across rho.

    python studies/severe-period-sensitivity/large_n_reference.py
"""

import math
import random
from statistics import NormalDist, median

N = NormalDist()
Z95 = N.inv_cdf(0.975)
Z999 = N.inv_cdf(0.999)
SEVERITY = {"1-in-100": N.inv_cdf(0.99), "1-in-1,000": N.inv_cdf(0.999)}
DRAWS = 20000


def fit(xs):
    mu = sum(xs) / len(xs)
    return mu, sum((x - mu) ** 2 for x in xs) / len(xs)


def estimates(mu, v):
    """rho-hat, probit(PD-hat), probit(q-hat)."""
    return v / (1 + v), mu / math.sqrt(1 + v), mu + Z999 * math.sqrt(v)


def standard_errors(_mu, v, periods):
    var_mu, var_v = v / periods, 2 * v * v / periods
    se_rho = math.sqrt(var_v) / (1 + v) ** 2
    se_q = math.sqrt(var_mu + (0.5 * Z999 / math.sqrt(v)) ** 2 * var_v)
    return se_rho, se_q


def cell(pd, rho, periods, u, k, seed):
    rng = random.Random(seed)
    m, s = N.inv_cdf(pd) / math.sqrt(1 - rho), math.sqrt(rho / (1 - rho))
    true_q = m + Z999 * s
    q_shift, rho_shift = [], []
    above = cover_before = cover_after = 0
    for _ in range(DRAWS):
        xs = [rng.gauss(m, s) for _ in range(periods)]
        before, after = fit(xs), fit(xs + [m + s * u] * k)
        (r0, _, q0), (r1, _, q1) = estimates(*before), estimates(*after)
        (se_r0, se_q0), (_, se_q1) = standard_errors(*before, periods), standard_errors(*after, periods + k)
        q_shift.append((q1 - q0) / se_q0)
        rho_shift.append((r1 - r0) / se_r0)
        above += q1 > q0 + Z95 * se_q0
        cover_before += abs(q0 - true_q) <= Z95 * se_q0
        cover_after += abs(q1 - true_q) <= Z95 * se_q1
    return median(q_shift), median(rho_shift), above / DRAWS, cover_before / DRAWS, cover_after / DRAWS


def main():
    print("| T | Added | q-hat shift (SE) | rho-hat shift (SE) | P(q-hat above old upper end) | q coverage before -> after |")
    print("|---|---|---|---|---|---|")
    shift = {}
    for periods in (20, 40, 100):
        for name, u in SEVERITY.items():
            for k in (1, 2):
                # PD 1%, rho 0.12; q's shift in SE units is the same for every PD and rho.
                q, r, above, c0, c1 = cell(0.01, 0.12, periods, u, k, seed=periods * 10 + k)
                shift[(periods, name, k)] = q
                label = f"{k} × {name}"
                print(f"| {periods} | {label} | {q:.2f} | {r:.2f} | {above:.3f} | {c0:.3f} -> {c1:.3f} |")
    print()
    print("Two periods over one, q-hat's shift: " + ", ".join(
        f"T = {p}, {n}: {shift[(p, n, 2)] / shift[(p, n, 1)]:.2f}" for p in (20, 40, 100) for n in SEVERITY))
    print()
    print("rho-hat's shift in SE units across rho (PD 1%), one period:")
    for periods in (20, 40, 100):
        for name, u in SEVERITY.items():
            vals = [cell(0.01, rho, periods, u, 1, seed=7)[1] for rho in (0.02, 0.12, 0.24)]
            print(f"  T = {periods}, {name}: " + ", ".join(f"{v:.2f}" for v in vals))


if __name__ == "__main__":
    main()
