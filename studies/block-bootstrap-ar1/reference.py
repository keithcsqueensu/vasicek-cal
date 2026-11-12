# SPDX-License-Identifier: Apache-2.0
"""S-53 (block bootstraps under an AR(1) factor): reference numbers from the model and pinned results only.

No study panel is simulated or fitted here. The inputs are the model, the pinned phi = 0 iid-bootstrap
coverage (tests/golden/recovery/summary.csv), and S-35's reference machinery
(studies/lrpd-serial-correlation/reference.py), imported unchanged.

What a resampling scheme captures. For the iid binomial MLE the per-period score s_t has lag-k covariance
C_k = sum_j a_j a_j' phi^(jk) (S-35), and the estimate's true variance is I^-1 B_T I^-1 with
B_T = I + 2 sum_{k<T} (1 - k/T) C_k. To first order a bootstrap of periods reproduces
B_w = I + 2 sum_{k<T} w_k C_k, with
  iid           w_k = 0;
  moving block  w_k = 1 - k/l for k < l, else 0 (Kunsch 1989; the Bartlett window), block length l;
  stationary    w_k = (1 - k/T) q^k + (k/T) q^(T-k), q = 1 - 1/b, mean block length b (Politis and
                Romano 1994, circular).
The shortfall of coordinate a is r_a = (I^-1 B_T I^-1)_aa / (I^-1 B_w I^-1)_aa, and the predicted coverage
maps the pinned phi = 0 iid-bootstrap coverage c0 through it with S-35's mapping:
Phi(z0/sqrt r - delta) - Phi(-z0/sqrt r - delta), z0 = Phi^-1((1 + c0)/2), delta S-35's extra downward
bias of rho (rho only). The end effects of a finite T (O(l/T)) are not modelled.

The data-driven block length (Politis and White 2004, with the correction of Patton, Politis and White
2009) is computed from the series x_t = Phi^-1((d_t + 1/2)/(n + 1)), the probit of each period's
continuity-corrected default rate (owner's choice (a), 2026-11-11). Its selector is noisy at small T, so
the reference does not plug in a single block length: it draws the selector's distribution on a Gaussian
AR(1) proxy series with x's own lag-k autocorrelation, R(k)/R(0) = s phi^k approximately (s the share of
Var x carried by its first Hermite term times the factor), 4,000 proxy series per (scenario, phi), numpy
PCG64 seed 0x53353342 ("S53B"). These are model draws, not study panels. Each arm's predicted coverage is
the mean of the mapped coverage over that distribution. The oracle block length (the selector's population
target from x's exact autocovariances) is reported beside it.

    uv run --no-project --with scipy==1.16.2 python studies/block-bootstrap-ar1/reference.py [--out FILE]
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from pathlib import Path

import numpy as np
from scipy.special import ndtr, ndtri

REPO = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("s35", REPO / "studies/lrpd-serial-correlation/reference.py")
s35 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s35)

PHIS = (0.0, 0.3, 0.6)
PROXY_DRAWS = 4000
PROXY_SEED = 0x53353342  # "S53B"
GROUPS = {
    "A": [0, 1, 2, 3, 4, 6, 7, 9, 10, 12, 15, 18, 19, 20, 21, 22, 23, 24, 27, 28, 30, 33, 36, 39, 45, 46, 48, 54,
          57, 60, 72],
    "B": [11, 29, 37, 38, 47, 55, 56, 63, 64, 65, 73, 74],
    "C": [5, 13, 14, 31, 32, 40, 41, 49, 50, 58, 59, 66, 67, 68, 75, 76, 77],
    "D": [8, 16, 17, 25, 26, 34, 35, 42, 43, 44, 51, 52, 53, 61, 62, 69, 70, 71, 78, 79, 80],
}
GROUP_OF = {s: g for g, ss in GROUPS.items() for s in ss}


# --- the Politis-White selector, as the engine will implement it (PREDICTION.md, "The block length") ---------

def politis_white(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(b_SB, b_CB) for each row of x (series of length n), before rounding, capped at b_max."""
    x = np.atleast_2d(x)
    r_count, n = x.shape
    x = x - x.mean(axis=1, keepdims=True)
    kn = max(5, math.ceil(math.sqrt(math.log10(n))))
    m_max = math.ceil(math.sqrt(n)) + kn
    b_max = math.ceil(min(3.0 * math.sqrt(n), n / 3.0))
    c = 2.0 * math.sqrt(math.log10(n) / n)
    acv = np.stack([(x[:, : n - k] * x[:, k:]).sum(axis=1) / n for k in range(m_max + 1)], axis=1)
    rho = acv / np.where(acv[:, :1] > 0.0, acv[:, :1], 1.0)
    insignificant = np.abs(rho) < c  # [r, k]
    m_hat = np.full(r_count, m_max)
    found = np.zeros(r_count, dtype=bool)
    for m in range(m_max - kn + 1):
        ok = insignificant[:, m + 1 : m + kn + 1].all(axis=1) & ~found
        m_hat[ok] = m
        found |= ok
    M = np.minimum(2 * np.maximum(m_hat, 1), m_max)
    b_sb, b_cb = np.ones(r_count), np.ones(r_count)
    for i in range(r_count):
        if acv[i, 0] <= 0.0:
            continue
        k = np.arange(-M[i], M[i] + 1)
        s = np.abs(k) / M[i]
        lam = np.where(s <= 0.5, 1.0, np.where(s <= 1.0, 2.0 * (1.0 - s), 0.0))
        a = acv[i, np.abs(k)]
        g = float((lam * a).sum())
        G = float((lam * np.abs(k) * a).sum())
        if g <= 0.0:
            continue
        b_sb[i] = min(((2.0 * G * G) / (2.0 * g * g)) ** (1.0 / 3.0) * n ** (1.0 / 3.0), b_max)
        b_cb[i] = min(((2.0 * G * G) / ((4.0 / 3.0) * g * g)) ** (1.0 / 3.0) * n ** (1.0 / 3.0), b_max)
    return b_sb, b_cb


def sb_mean_block(b_sb: float) -> float:
    return max(1.0, b_sb)


def mbb_length(b_cb: float) -> int:
    return max(1, int(math.floor(b_cb + 0.5)))


# --- the model -------------------------------------------------------------------------------------------

def score_model(pd: float, rho: float, n: int):
    """I, I^-1 and the Hermite coefficients a_j of E[s | z] for the iid binomial MLE (S-35's construction),
    and the probit-rate series' marginal variance and Hermite coefficients."""
    d = s35.support(n, pd, rho)
    u = np.array([math.log(pd / (1 - pd)), math.log(rho / (1 - rho))])
    step = 1e-4

    def logf(uu):
        p_ = 1 / (1 + math.exp(-uu[0]))
        r_ = 1 / (1 + math.exp(-uu[1]))
        return np.log((s35.binom_given_z(n, s35.cond_pd(p_, r_), d) * s35.WZ).sum(axis=1))

    s = np.empty((2, d.size))
    for a in range(2):
        e = np.zeros(2)
        e[a] = step
        s[a] = (logf(u + e) - logf(u - e)) / (2 * step)
    k_dz = s35.binom_given_z(n, s35.cond_pd(pd, rho), d)
    f_d = (k_dz * s35.WZ).sum(axis=1)
    info = (s[:, None, :] * s[None, :, :] * f_d).sum(axis=2)
    h = s @ k_dz
    a = (h[:, None, :] * s35.HE[None, :, :] * s35.WZ).sum(axis=2)
    x = ndtri((d + 0.5) / (n + 1.0))  # the probit rate of each count
    ex = float((x * f_d).sum())
    var_x = float(((x - ex) ** 2 * f_d).sum())
    gx = x @ k_dz  # E[x | z]
    cx = (gx[None, :] * s35.HE * s35.WZ).sum(axis=1)
    return info, np.linalg.inv(info), a, var_x, cx


def lag_cov(a: np.ndarray, phi: float, k: int) -> np.ndarray:
    return sum(np.outer(a[:, j], a[:, j]) * phi ** (j * k) for j in range(1, s35.HERMITE + 1))


_CKS: dict = {}


def shortfall(info, inv, a, phi: float, t: int, w) -> tuple[float, float]:
    key = (id(a), phi, t)
    if key not in _CKS:
        _CKS[key] = [lag_cov(a, phi, k) if phi > 0 else np.zeros((2, 2)) for k in range(1, t)]
    cks = _CKS[key]
    b_true = info + 2 * sum((1 - k / t) * cks[k - 1] for k in range(1, t))
    b_w = info + 2 * sum(w(k) * cks[k - 1] for k in range(1, t))
    v, vw = inv @ b_true @ inv, inv @ b_w @ inv
    return v[0, 0] / vw[0, 0], v[1, 1] / vw[1, 1]


def w_mbb(l: int):
    return lambda k: (1 - k / l) if k < l else 0.0


def w_sb(b: float, t: int):
    q = 1.0 - 1.0 / b
    return lambda k: (1 - k / t) * q ** k + (k / t) * q ** (t - k)


def x_acv(var_x: float, cx: np.ndarray, phi: float, k: int) -> float:
    if k == 0:
        return var_x
    return float(sum(cx[j] ** 2 * phi ** (j * k) for j in range(1, s35.HERMITE + 1)))


def oracle_lengths(var_x, cx, phi, t) -> tuple[float, float]:
    """The selector's population target from x's exact autocovariances (no lag window), capped at b_max."""
    b_max = math.ceil(min(3.0 * math.sqrt(t), t / 3.0))
    if phi == 0.0:
        return 1.0, 1.0
    r = [x_acv(var_x, cx, phi, k) for k in range(0, 2000)]
    g = r[0] + 2 * sum(r[1:])
    G = 2 * sum(k * r[k] for k in range(1, len(r)))
    b_sb = ((2 * G * G) / (2 * g * g)) ** (1 / 3) * t ** (1 / 3)
    b_cb = ((2 * G * G) / ((4 / 3) * g * g)) ** (1 / 3) * t ** (1 / 3)
    return min(b_sb, b_max), min(b_cb, b_max)


def proxy_lengths(var_x, cx, phi, t, rng) -> tuple[np.ndarray, np.ndarray]:
    """The selector on PROXY_DRAWS Gaussian series with x's lag-1 autocorrelation r1 as an AR(1) signal plus
    white noise: x_t = sqrt(s) y_t + sqrt(1 - s) e_t, y an AR(1) in phi, s = r1 / phi."""
    sig = 0.0 if phi == 0.0 else min(1.0, x_acv(var_x, cx, phi, 1) / var_x / phi)
    e = rng.standard_normal((PROXY_DRAWS, t))
    v = rng.standard_normal((PROXY_DRAWS, t))
    y = np.empty((PROXY_DRAWS, t))
    y[:, 0] = v[:, 0]
    for j in range(1, t):
        y[:, j] = phi * y[:, j - 1] + math.sqrt(1 - phi * phi) * v[:, j]
    return politis_white(math.sqrt(sig) * y + math.sqrt(1 - sig) * e)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "studies/block-bootstrap-ar1/reference.csv"))
    args = ap.parse_args()
    pinned = {int(r["scenario"]): r for r in s35.read_summary(REPO / "tests/golden/recovery/summary.csv")}
    rng = np.random.Generator(np.random.PCG64(PROXY_SEED))
    rows = []
    models = {}
    for sid in range(81):
        pd, rho, t, n = s35.scenario(sid)
        if (pd, rho, n) not in models:
            models[(pd, rho, n)] = score_model(pd, rho, n)
        info, inv, a, var_x, cx = models[(pd, rho, n)]
        r = pinned[sid]
        c0 = (float(r["pd_boot_coverage"]), float(r["rho_boot_coverage"]))
        l_fixed = math.ceil(t ** (1 / 3) - 1e-12)
        for phi in PHIS:
            f_mu, f_s = s35.f_ar(phi, t), s35.f_ar(phi, t, 2)
            delta = ((f_mu - 1) / t) / math.sqrt(2 * f_s / t) if phi > 0 else 0.0
            row = {"scenario": sid, "group": GROUP_OF[sid], "pd": pd, "rho": rho, "periods": t, "obligors": n,
                   "phi": phi, "boot_pd_c0": c0[0], "boot_rho_c0": c0[1], "l_fixed": l_fixed}

            def pred(w):
                sp, sr = shortfall(info, inv, a, phi, t, w)
                return s35.coverage(c0[0], sp), s35.coverage(c0[1], sr, delta), sp, sr

            for name, w in (("iid", lambda k: 0.0), ("mbb_fixed", w_mbb(l_fixed))):
                cp, cr, sp, sr = pred(w)
                row.update({f"{name}_pd": cp, f"{name}_rho": cr, f"{name}_short_pd": sp, f"{name}_short_rho": sr})
            o_sb, o_cb = oracle_lengths(var_x, cx, phi, t)
            row.update({"oracle_b_sb": o_sb, "oracle_b_cb": o_cb})
            b_sb, b_cb = proxy_lengths(var_x, cx, phi, t, rng)
            row.update({"pw_b_sb_median": float(np.median(b_sb)), "pw_b_cb_median": float(np.median(b_cb)),
                        "pw_mbb_l_is_1": float(np.mean([mbb_length(b) == 1 for b in b_cb]))})
            cache_sb, cache_cb = {}, {}
            pd_sb = rho_sb = pd_cb = rho_cb = 0.0
            for bs, bc in zip(b_sb, b_cb):
                kb = round(sb_mean_block(bs), 2)
                if kb not in cache_sb:
                    cache_sb[kb] = pred(w_sb(kb, t)) if kb > 1.0 else pred(lambda k: 0.0)
                lb = mbb_length(bc)
                if lb not in cache_cb:
                    cache_cb[lb] = pred(w_mbb(lb))
                pd_sb += cache_sb[kb][0]
                rho_sb += cache_sb[kb][1]
                pd_cb += cache_cb[lb][0]
                rho_cb += cache_cb[lb][1]
            row.update({"mbb_pw_pd": pd_cb / PROXY_DRAWS, "mbb_pw_rho": rho_cb / PROXY_DRAWS,
                        "sb_pw_pd": pd_sb / PROXY_DRAWS, "sb_pw_rho": rho_sb / PROXY_DRAWS})
            rows.append(row)
        print(f"scenario {sid} done", flush=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in row.items()})
    print("group means (PD / rho) by phi:")
    for g in "ABCD":
        for phi in PHIS:
            sel = [r for r in rows if r["group"] == g and r["phi"] == phi]
            txt = ", ".join(f"{a} {np.mean([r[a + '_pd'] for r in sel]):.3f}/{np.mean([r[a + '_rho'] for r in sel]):.3f}"
                            for a in ("iid", "mbb_fixed", "mbb_pw", "sb_pw"))
            print(f"  {g} phi {phi}: {txt}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
