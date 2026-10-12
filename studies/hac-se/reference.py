# SPDX-License-Identifier: Apache-2.0
"""S-48 (HAC SE under serial dependence): reference numbers from the model and pinned results only.

No panel is simulated or fitted here. The inputs are the model, S-35's score structure
(studies/lrpd-serial-correlation/reference.py: per period, the binomial MLE's score s(d) in (logit PD,
logit rho), its information I = E[s s'], and the lag-k covariance C_k = sum_j a_j a_j' phi^(jk) under the
AR(1) factor), and the pinned phi = 0 Wald coverage of the recovery summary.

For axis a, the estimate's error is about (T I)^-1 sum_t s_t, so it is a sum of the scalar series
x_t = q' s_t, q = I^-1 e_a, whose T x T covariance Sigma is Toeplitz with entries q' C_|t-r| q
(C_0 = I). Then, exactly for these second moments:
    true variance (times T^2)    v = 1' Sigma 1
    iid SE's variance (times T^2) T q' I q              (the information; the Hessian's mean)
    HAC variance (times T^2)     Q = x' A x,  A = M K M, M = I_T - 11'/T (centring), K_tr = k(|t - r|)
                                 the Bartlett weights k(m) = 1 - m/(L + 1) for m <= L, else 0;
                                 E[Q] = tr(A Sigma), and for Gaussian x, Var[Q] = 2 tr((A Sigma)^2).
    the iid inflation f = v / (T q' I q) (S-35's f); the HAC's remaining inflation g = v / E[Q];
    the HAC's effective degrees of freedom nu = 2 E[Q]^2 / Var[Q] (Satterthwaite: Q ~ E[Q] chi2_nu / nu).
The expected SE ratio HAC / iid is sqrt(E[Q] / (T q' I q)) = sqrt(f / g).

Coverage mapping, a heuristic held loosely (as S-35's): the scenario's pinned phi = 0 Wald coverage c0
gives the iid Wald's effective critical value z0 = Phi^-1((1 + c0)/2), which carries the small-T
non-normality of the Wald pivot. The HAC Wald interval's coverage is then
    F_nu(z0 / sqrt(g) - delta) - F_nu(-z0 / sqrt(g) - delta),
F_nu Student's t: the HAC SE's own variability on top of the bias g. delta (rho only) is S-35's extra
downward bias of rho-hat from the dependent mean, in standard deviations of the estimate.

Pooled per group (A-D as S-35) and phi, weighting each scenario by its pinned unflagged count, the
replicates the Wald (and HAC) interval is computed on.

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python studies/hac-se/reference.py [--out FILE]
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import ndtri
from scipy.stats import t as student_t

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "studies" / "lrpd-serial-correlation"))
from reference import (
    HE,
    HERMITE,
    WZ,
    binom_given_z,
    cond_pd,
    f_ar,
    scenario,
    support,
)
from score_m6 import group

PHIS = (0.0, 0.3, 0.6)


def lag_rule(t: int) -> int:
    return math.floor(4.0 * (t / 100.0) ** (2.0 / 9.0))


def score_structure(pd: float, rho: float, n: int):
    """I and the Hermite coefficients a_j of E[s | z], as S-35's binomial_inflation computes them."""
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
    k_dz = binom_given_z(n, cond_pd(pd, rho), d)
    f_d = (k_dz * WZ).sum(axis=1)
    info = (s[:, None, :] * s[None, :, :] * f_d).sum(axis=2)
    h = s @ k_dz
    a = (h[:, None, :] * HE[None, :, :] * WZ).sum(axis=2)
    return info, a[:, 1:]  # a_j for j = 1..HERMITE


def hac_moments(info, a, phi: float, t: int, lag: int, axis: int):
    """f, g, nu and the expected SE ratio for one axis (see the module docstring)."""
    q = np.linalg.solve(info, np.eye(2)[axis])
    sig0 = float(q @ info @ q)
    qa = a.T @ q  # q' a_j, j = 1..J
    lags = np.arange(t)
    acov = np.array([sig0] + [float(np.sum(qa * qa * phi ** (np.arange(1, HERMITE + 1) * k))) for k in lags[1:]])
    Sigma = acov[np.abs(lags[:, None] - lags[None, :])]
    m = np.abs(lags[:, None] - lags[None, :])
    K = np.where(m <= lag, 1.0 - m / (lag + 1.0), 0.0)
    M = np.eye(t) - 1.0 / t
    AS = M @ K @ M @ Sigma
    eq = float(np.trace(AS))
    var_q = 2.0 * float(np.trace(AS @ AS))
    v = float(Sigma.sum())
    f = v / (t * sig0)
    g = v / eq
    return {"f": f, "g": g, "nu": 2.0 * eq * eq / var_q, "se_ratio": math.sqrt(eq / (t * sig0))}


def mapped(c0: float, g: float, nu: float, delta: float = 0.0) -> float:
    c0 = min(max(c0, 1e-6), 1 - 1e-6)
    z0 = float(ndtri((1 + c0) / 2))
    w = z0 / math.sqrt(g)
    return float(student_t.cdf(w - delta, nu) - student_t.cdf(-w - delta, nu))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "studies/hac-se/reference.csv"))
    args = ap.parse_args()
    pinned = {int(r["scenario"]): r for r in read_csv(REPO / "tests/golden/recovery/summary.csv")}
    cache = {}
    rows = []
    for sid in range(81):
        pd, rho, t, n = scenario(sid)
        if (pd, rho, n) not in cache:
            cache[(pd, rho, n)] = score_structure(pd, rho, n)
        info, a = cache[(pd, rho, n)]
        r = pinned[sid]
        unflagged = int(r["unflagged"])
        rule = lag_rule(t)
        for phi in PHIS:
            f_mu, f_s = f_ar(phi, t), f_ar(phi, t, 2)
            delta = ((f_mu - 1) / t) / math.sqrt(2 * f_s / t) if phi > 0 else 0.0
            row = {"scenario": sid, "group": group(sid, t), "pd": pd, "rho": rho, "periods": t, "obligors": n,
                   "phi": phi, "unflagged": unflagged, "lag": rule}
            for axis, name in enumerate(("pd", "rho")):
                c0 = float(r[f"{name}_coverage"]) if unflagged > 0 else float("nan")
                for lag, tag in ((rule, ""), (1, "_l1"), (2 * rule, "_l2x")):
                    mo = hac_moments(info, a, phi, t, lag, axis)
                    row[f"{name}_g{tag}"] = mo["g"]
                    row[f"{name}_nu{tag}"] = mo["nu"]
                    row[f"{name}_se_ratio{tag}"] = mo["se_ratio"]
                    row[f"{name}_hac_pred{tag}"] = (mapped(c0, mo["g"], mo["nu"], delta if axis == 1 else 0.0)
                                                    if unflagged > 0 else float("nan"))
                    if tag == "":
                        row[f"{name}_f"] = mo["f"]
                row[f"{name}_wald_c0"] = c0
                row[f"{name}_iid_pred"] = (mapped(c0, row[f"{name}_f"], math.inf, delta if axis == 1 else 0.0)
                                           if unflagged > 0 else float("nan"))
            rows.append(row)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in row.items()})

    print("Pooled per group and phi (weights: pinned unflagged counts). Wald: the iid Wald mapped with f;")
    print("HAC: mapped with g and nu at the rule's lag (L = 1 and 2 x rule in brackets). SE ratio: median of")
    print("the expected HAC/iid SE ratio over the group's scenarios; pooled RMS: the root of its weighted")
    print("mean square, the expected root mean square of the per-fit ratio.")
    for name in ("pd", "rho"):
        print(f"{name.upper()}:")
        for grp in "ABCD":
            for phi in PHIS:
                sel = [r for r in rows if r["group"] == grp and r["phi"] == phi and r["unflagged"] > 0]
                wsum = sum(r["unflagged"] for r in sel)

                def pool(key, sel=sel, wsum=wsum):
                    return sum(r[key] * r["unflagged"] for r in sel) / wsum

                print(f"  {grp} phi {phi}: Wald {pool(f'{name}_iid_pred'):.3f} | HAC {pool(f'{name}_hac_pred'):.3f} "
                      f"[{pool(f'{name}_hac_pred_l1'):.3f}, {pool(f'{name}_hac_pred_l2x'):.3f}] | "
                      f"f {np.median([r[f'{name}_f'] for r in sel]):.2f}, g {np.median([r[f'{name}_g'] for r in sel]):.2f}"
                      f", nu {np.median([r[f'{name}_nu'] for r in sel]):.1f}, SE ratio "
                      f"{np.median([r[f'{name}_se_ratio'] for r in sel]):.3f} (pooled RMS "
                      f"{math.sqrt(sum(r[f'{name}_se_ratio'] ** 2 * r['unflagged'] for r in sel) / wsum):.3f})")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
