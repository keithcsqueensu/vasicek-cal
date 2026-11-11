"""Golden values for the multi-grade binomial-mixture period log-likelihood (core/objectives/multigrade.hpp; M6,
D-151, D-411).

likelihood.csv: log E[ prod_k p_k(Z)^d_k (1 - p_k(Z))^(n_k - d_k) ] for K grades on one factor with a shared rho,
p_k(z) = Phi((Phi^-1(PD_k) - sqrt(rho) z) / sqrt(1 - rho)): the period log-likelihood without its binomial
coefficients. Each value is computed in mpmath at 50 significant digits and rounded once to binary64. The integral is
mp.quad on h(z) = log of the integrand plus log phi, split at the mode +- k sigma and at every one-sided grade's own
transition point (where its survival factor is 1/2), so a steep d = 0 cliff beside the peak falls on a split point.

The cases (D-411 sec. 2):
- two-grade interior periods across PD, rho and n;
- the named risk: a d = 0 grade of n up to 10^6 beside interior grades, rho up to 0.9;
- a d = 0 grade beside a d = n grade (interior by the classification);
- one-sided periods, every grade d = 0 or every grade d = n, K up to 5;
- D-151's 17-grade master scale, PD_k = 0.03% x 1.5^(k - 1), 10,000 obligors spread by a discretised normal of
  standard deviation 4 grades, at a benign, an average and an adverse year, rho 0.02, 0.12 and 0.24;
- five grades at rho up to 0.9.

Self-checks before writing: K = 1 equals the single-grade integral of tools/gen_quadrature_tables.py's method, and two
grades with one PD equal the pooled single grade plus the log ratio of binomial coefficients, to 1e-30.

    uv run --no-project --with mpmath==1.4.1 python tools/gen_multigrade_goldens.py
    uv run --no-project --with mpmath==1.4.1 python tools/gen_multigrade_goldens.py --check
"""

from __future__ import annotations

import sys

from mpmath import mp, mpf
from tablegen import csv_note, log_phi_mp, probit_mp, row, run, to_double

GENERATOR = "tools/gen_multigrade_goldens.py"
MANIFEST = "tests/golden/multigrade/MANIFEST.json"
OUT = "tests/golden/multigrade"


class MultiGrade:
    """h(z) = sum_k [d_k log Phi(x_k) + (n_k - d_k) log Phi(-x_k)] - z^2/2, x_k = (c_k - sqrt(rho) z)/sqrt(1-rho)."""

    def __init__(self, pds: list[float], rho: float, grades: list[tuple[int, int]]):
        self.c = [probit_mp(mpf(p)) for p in pds]
        self.sr = mp.sqrt(mpf(rho))
        self.s1r = mp.sqrt(1 - mpf(rho))
        self.g = grades

    def x(self, k: int, z: mpf) -> mpf:
        return (self.c[k] - self.sr * z) / self.s1r

    def h(self, z: mpf) -> mpf:
        v = -z * z / 2
        for k, (n, d) in enumerate(self.g):
            if n == 0:
                continue
            x = self.x(k, z)
            if d:
                v += d * log_phi_mp(x)
            if n - d:
                v += (n - d) * log_phi_mp(-x)
        return v

    def dh(self, z: mpf) -> mpf:
        lam = lambda t: mp.npdf(t) / mp.ncdf(t)
        s = mpf(0)
        for k, (n, d) in enumerate(self.g):
            if n == 0:
                continue
            x = self.x(k, z)
            s += d * lam(x) - (n - d) * lam(-x)
        return -(self.sr / self.s1r) * s - z

    def mode(self) -> mpf:
        lo, hi = mpf(-1), mpf(1)
        while self.dh(lo) < 0:
            lo *= 2
        while self.dh(hi) > 0:
            hi *= 2
        return mp.findroot(self.dh, (lo, hi), solver="anderson")

    def transition_points(self) -> list[mpf]:
        """Each one-sided grade's z where its survival factor is 1/2: (1-p)^n = 1/2 for d = 0, p^n = 1/2 for d = n."""
        out = []
        for k, (n, d) in enumerate(self.g):
            if n == 0 or 0 < d < n:
                continue
            q = -mp.expm1(-mp.log(2) / n)
            x_c = probit_mp(q) if d == 0 else -probit_mp(q)
            out.append((self.c[k] - self.s1r * x_c) / self.sr)
        return out

    def log_integral(self) -> mpf:
        zs = self.mode()
        sigma = 1 / mp.sqrt(-mp.diff(self.h, zs, 2))
        h0 = self.h(zs)
        pts = {zs + k * sigma for k in (-60, -30, -15, -8, -4, -2, 0, 2, 4, 8, 15, 30, 60)}
        for zc in self.transition_points():
            if zs - 60 * sigma < zc < zs + 60 * sigma:
                pts.add(zc)
        pts = [-mp.inf] + sorted(pts) + [mp.inf]
        val, err = mp.quad(lambda z: mp.exp(self.h(z) - h0), pts, error=True)
        assert err <= mpf(10) ** -35 * val, ("quad error", self.g, err)
        return h0 + mp.log(val) - mp.log(2 * mp.pi) / 2


def master_scale() -> tuple[list[float], list[int]]:
    """D-151: 17 grades, PD_k = 0.03% x 1.5^(k-1); 10,000 obligors by a discretised normal of sd 4 grades, centred on
    grade 9, rounded to integers with the remainder on the centre grade."""
    pds = [float(mpf("0.0003") * mpf("1.5") ** k) for k in range(17)]
    w = [mp.exp(-((k - 8) ** 2) / (2 * 16)) for k in range(17)]
    tot = sum(w)
    ns = [int(mp.nint(10000 * wk / tot)) for wk in w]
    ns[8] += 10000 - sum(ns)
    return pds, ns


def cases() -> list[tuple[str, list[float], float, list[tuple[int, int]]]]:
    out = []
    # two-grade interior periods
    for rho in (0.02, 0.12, 0.24, 0.5):
        out.append(("two", [0.004, 0.03], rho, [(1000, 3), (500, 20)]))
        out.append(("two", [1e-4, 0.2], rho, [(10**4, 2), (100, 25)]))
        out.append(("two", [0.01, 0.01], rho, [(10**6, 9000), (10**4, 150)]))
    # the named risk: a d = 0 grade of large n beside interior grades
    for rho in (0.02, 0.12, 0.24, 0.5, 0.9):
        for n0 in (10**4, 10**5, 10**6):
            out.append(("risk", [0.01, 1e-4], rho, [(1000, 10), (n0, 0)]))
            out.append(("risk", [0.02, 0.003, 3e-4], rho, [(500, 12), (5000, 20), (n0, 0)]))
    # d = 0 beside d = n (interior by the classification)
    for rho in (0.12, 0.5):
        out.append(("zero_full", [0.01, 0.6], rho, [(100, 0), (5, 5)]))
    # one-sided periods
    for rho in (0.02, 0.24, 0.9):
        out.append(("one_sided_zero", [1e-3, 1e-2, 3e-4], rho, [(1000, 0), (50, 0), (10**5, 0)]))
        out.append(("one_sided_zero", [0.05, 0.002, 0.01, 1e-4, 0.03], rho,
                    [(20, 0), (2000, 0), (300, 0), (10**6, 0), (10, 0)]))
        out.append(("one_sided_full", [0.6, 0.9], rho, [(3, 3), (10, 10)]))
    # D-151's 17-grade master scale: benign, average and adverse years
    pds, ns = master_scale()
    for rho in (0.02, 0.12, 0.24):
        for stress in ("0.3", "1", "3"):
            grades = [(n, int(mp.nint(n * mpf(p) * mpf(stress)))) for n, p in zip(ns, pds)]
            out.append((f"master_x{stress}", pds, rho, grades))
    # five grades at large rho
    for rho in (0.5, 0.75, 0.9):
        out.append(("five", [3e-4, 0.002, 0.01, 0.05, 0.2], rho, [(5000, 1), (2000, 3), (800, 9), (300, 14), (60, 13)]))
    return out


def self_checks() -> None:
    # K = 1: the same integral as a single grade (one term in h)
    a = MultiGrade([1e-3], 0.12, [(10**6, 1000)]).log_integral()
    b = MultiGrade([1e-3, 0.3], 0.12, [(10**6, 1000), (0, 0)]).log_integral()
    assert abs(a - b) <= mpf(10) ** -30 * abs(a), "K = 1 with an empty grade"
    # pooling: equal PDs, Binom(d1; n1) Binom(d2; n2) = [C(n1,d1) C(n2,d2) / C(n1+n2, d1+d2)] Binom(d1+d2; n1+n2)
    # in the integrand without coefficients the ratio cancels: E[p^(d1+d2) (1-p)^(...)] is the pooled single grade.
    two = MultiGrade([0.02, 0.02], 0.24, [(400, 7), (150, 4)]).log_integral()
    pooled = MultiGrade([0.02], 0.24, [(550, 11)]).log_integral()
    assert abs(two - pooled) <= mpf(10) ** -30 * abs(pooled), "pooling"


def build_outputs() -> dict[str, str]:
    self_checks()
    lines = [csv_note(GENERATOR), "kind,rho_hex,grades,pd_hex,log_integral_hex,rho,log_integral\n"]
    for kind, pds, rho, grades in cases():
        v = MultiGrade(pds, rho, grades).log_integral()
        g = ";".join(f"{n}:{d}" for n, d in grades)
        p = ";".join(float(x).hex() for x in pds)
        lines.append(row(kind, float(rho).hex(), g, p, to_double(v).hex(), repr(rho), mp.nstr(v, 20)))
    return {f"{OUT}/likelihood.csv": "".join(lines)}


if __name__ == "__main__":
    sys.exit(run(GENERATOR, MANIFEST, build_outputs, __doc__))
