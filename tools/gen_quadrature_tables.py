#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Generate the Gauss-Hermite tables and the quadrature golden tables.

Nodes and log-weights of the standard-normal Gauss-Hermite rules for every
supported N (D-080) come from mpmath's physicists' rule (weight exp(-x^2)),
converted in one place: z = sqrt(2) x, w = W / sqrt(pi). Weights are stored as
logs, so the tiny outer weights keep full relative precision. Before anything
is written, the *converted* rule must pass, at 50 digits:

  - symmetry: z_i = -z_(N-1-i) and w_i = w_(N-1-i);
  - moments: sum w_i z_i^(2k) = (2k-1)!! for 2k <= 2N-2, and odd moments vanish;
  - weights: each w_i matches the closed form
      w_i = 2^(N-1) N! / (N^2 H_(N-1)(x_i)^2)   (divided by sqrt(pi) in the conversion),
    an independent derivation aimed at the tiny outer weights.

Golden tables (references for tests/unit/quadrature_test.cpp):
  gh_moments.csv    log (2k-1)!! for k = 0..127, exact moments of N(0,1)
  mixture.csv       log E[p(Z)^d (1-p(Z))^(n-d)] for the one-factor Vasicek p(z):
                    closed forms at n <= 2 (E[p] = PD, E[p^2] = Phi2 via Plackett's
                    identity), and spike integrands up to n = 1e6 by adaptive
                    tanh-sinh around the mode.

Run from the repository root:
    uv run --no-project --with mpmath==1.4.1 python tools/gen_quadrature_tables.py
    uv run --no-project --with mpmath==1.4.1 python tools/gen_quadrature_tables.py --check
"""

from __future__ import annotations

import sys

from mpmath import mp, mpf
from tablegen import (
    csv_note,
    cxx,
    header_banner,
    log_phi_mp,
    probit_mp,
    row,
    run,
    to_double,
)

GENERATOR = "tools/gen_quadrature_tables.py"
GEN_DIR = "core/quadrature/generated"
GOLDEN_DIR = "tests/golden/quadrature"
MANIFEST = f"{GOLDEN_DIR}/MANIFEST.json"
CSV_NOTE = csv_note(GENERATOR)

GH_SIZES = [8, 16, 24, 32, 48, 64, 96, 128, 256]  # D-080, D-088
GL_SIZES = [8, 12, 16, 20, 24, 32]  # composite Gauss-Legendre for one-sided integrands (D-118)
CHECK_REL = mpf(10) ** -40


# --- Gauss-Hermite rules -------------------------------------------------------------


def standard_normal_rule(n: int) -> tuple[list[mpf], list[mpf]]:
    """Nodes z and log-weights of the n-point rule for E[f(Z)], Z ~ N(0,1).

    The only place the physicists' -> standard-normal conversion happens.
    """
    xs, ws = mp.gauss_quadrature(n, "hermite")
    pairs = sorted(zip(xs, ws), key=lambda p: p[0])
    z = [mp.sqrt(2) * x for x, _ in pairs]
    log_w = [mp.log(w) - mp.log(mp.pi) / 2 for _, w in pairs]
    return z, log_w


def check_rule(n: int, z: list[mpf], log_w: list[mpf]) -> None:
    w = [mp.exp(v) for v in log_w]
    for i in range(n):
        assert abs(z[i] + z[n - 1 - i]) <= CHECK_REL * max(1, abs(z[i])), ("node symmetry", n, i)
        assert abs(w[i] / w[n - 1 - i] - 1) <= CHECK_REL, ("weight symmetry", n, i)
    for k in range(n):
        even = mp.fsum(wi * zi ** (2 * k) for wi, zi in zip(w, z))
        assert abs(even / mp.fac2(2 * k - 1) - 1) <= CHECK_REL, ("even moment", n, 2 * k)
        odd = [wi * zi ** (2 * k + 1) for wi, zi in zip(w, z)]
        assert abs(mp.fsum(odd)) <= CHECK_REL * mp.fsum(abs(t) for t in odd), ("odd moment", n, 2 * k + 1)
    for zi, lwi in zip(z, log_w):
        x = zi / mp.sqrt(2)
        h_prev, h = mpf(1), 2 * x  # physicists' H_0, H_1 -> H_(n-2), H_(n-1)
        for k in range(1, n - 1):
            h_prev, h = h, 2 * x * h - 2 * k * h_prev
        log_formula = (n - 1) * mp.log(2) + mp.log(mp.factorial(n)) - 2 * mp.log(n) - 2 * mp.log(abs(h))
        assert abs(lwi - log_formula) <= CHECK_REL * max(1, abs(log_formula)), ("closed-form weight", n, zi)


def render_gauss_hermite(rules: dict[int, tuple[list[float], list[float]]], smallest: dict[int, float]) -> str:
    sizes = ", ".join(str(n) for n in GH_SIZES)
    out = [
        header_banner(GENERATOR, MANIFEST),
        (
            "//\n"
            "// Standard-normal Gauss-Hermite rules: E[f(Z)] ~ sum_i exp(log_weight_i) f(node_i), Z ~ N(0,1),\n"
            "// exact for polynomials of degree <= 2N-1. Nodes ascend; weights are stored as logs so the\n"
            "// tiny outer ones keep full relative precision (D-080). Checked at 50 digits before writing:\n"
            "// symmetry, all moments to degree 2N-1, and each weight against the closed form.\n"
            "//\n"
            "// Host-side tables. Backends hand them to integrators by pointer; the CUDA backend copies\n"
            "// them to device memory (namespace-scope arrays are not usable from device code).\n"
            "#pragma once\n\n"
            "namespace vcal::quadrature::generated {\n\n"
            f"inline constexpr int kGaussHermiteSizes[{len(GH_SIZES)}] = {{{sizes}}};\n\n"
        ),
    ]
    for n in GH_SIZES:
        z, lw = rules[n]
        out.append(f"// N = {n}: smallest weight {smallest[n]:.3e}\n")
        out.append(f"inline constexpr double kGaussHermiteNode{n}[{n}] = {{\n")
        out.extend(f"    {cxx(v)},\n" for v in z)
        out.append("};\n")
        out.append(f"inline constexpr double kGaussHermiteLogWeight{n}[{n}] = {{\n")
        out.extend(f"    {cxx(v)},\n" for v in lw)
        out.append("};\n\n")
    out.append(
        "// Sets *node and *log_weight for a supported n and returns true; otherwise returns false.\n"
        "inline bool gauss_hermite_lookup(int n, const double** node, const double** log_weight) {\n"
        "    switch (n) {\n"
    )
    for n in GH_SIZES:
        out.append(
            f"        case {n}: *node = kGaussHermiteNode{n}; *log_weight = kGaussHermiteLogWeight{n}; return true;\n"
        )
    out.append("        default: return false;\n    }\n}\n\n}  // namespace vcal::quadrature::generated\n")
    return "".join(out)


# --- Gauss-Legendre rules (D-118) -------------------------------------------------------------


def legendre_rule(n: int) -> tuple[list[mpf], list[mpf]]:
    """Nodes and weights on [-1, 1], weight 1; checked before use."""
    xs, ws = mp.gauss_quadrature(n, "legendre")
    pairs = sorted(zip(xs, ws), key=lambda p: p[0])
    x = [a for a, _ in pairs]
    w = [b for _, b in pairs]
    for i in range(n):
        assert abs(x[i] + x[n - 1 - i]) <= CHECK_REL, ("GL node symmetry", n, i)
        assert abs(w[i] / w[n - 1 - i] - 1) <= CHECK_REL, ("GL weight symmetry", n, i)
    for k in range(n):  # exact for x^m, m <= 2n - 1: integral of x^(2k) on [-1, 1] is 2/(2k+1)
        even = mp.fsum(wi * xi ** (2 * k) for wi, xi in zip(w, x))
        assert abs(even * (2 * k + 1) / 2 - 1) <= CHECK_REL, ("GL even moment", n, 2 * k)
        odd = [wi * xi ** (2 * k + 1) for wi, xi in zip(w, x)]
        assert abs(mp.fsum(odd)) <= CHECK_REL * mp.fsum(abs(t) for t in odd), ("GL odd moment", n, 2 * k + 1)
    return x, w


def render_gauss_legendre(rules: dict[int, tuple[list[float], list[float]]]) -> str:
    sizes = ", ".join(str(n) for n in GL_SIZES)
    out = [
        header_banner(GENERATOR, MANIFEST),
        (
            "//\n"
            "// Gauss-Legendre rules on [-1, 1]: integral of f ~ sum_i weight_i f(node_i), exact for\n"
            "// polynomials of degree <= 2M-1. Checked at 50 digits before writing: symmetry and every\n"
            "// moment. Host-side tables, handed to integrators by pointer (as for Gauss-Hermite).\n"
            "#pragma once\n\n"
            "namespace vcal::quadrature::generated {\n\n"
            f"inline constexpr int kGaussLegendreSizes[{len(GL_SIZES)}] = {{{sizes}}};\n\n"
        ),
    ]
    for n in GL_SIZES:
        x, w = rules[n]
        out.append(f"inline constexpr double kGaussLegendreNode{n}[{n}] = {{\n")
        out.extend(f"    {cxx(v)},\n" for v in x)
        out.append("};\n")
        out.append(f"inline constexpr double kGaussLegendreWeight{n}[{n}] = {{\n")
        out.extend(f"    {cxx(v)},\n" for v in w)
        out.append("};\n\n")
    out.append(
        "// Sets *node and *weight for a supported n and returns true; otherwise returns false.\n"
        "inline bool gauss_legendre_lookup(int n, const double** node, const double** weight) {\n"
        "    switch (n) {\n"
    )
    for n in GL_SIZES:
        out.append(f"        case {n}: *node = kGaussLegendreNode{n}; *weight = kGaussLegendreWeight{n}; return true;\n")
    out.append("        default: return false;\n    }\n}\n\n}  // namespace vcal::quadrature::generated\n")
    return "".join(out)


# --- golden: moments -----------------------------------------------------------------


def golden_moments() -> str:
    lines = [CSV_NOTE, "k,log_moment_hex,log_moment\n"]  # log E[Z^(2k)] = log (2k-1)!!
    for k in range(max(GH_SIZES)):
        e = to_double(mp.log(mp.fac2(2 * k - 1)))
        lines.append(row(k, e.hex(), repr(e)))
    return "".join(lines)


# --- golden: binomial-mixture integrals ----------------------------------------------


class Mixture:
    """h(z) = d log Phi(x) + (n-d) log Phi(-x) - z^2/2, x = (c - sqrt(rho) z) / sqrt(1-rho)."""

    def __init__(self, pd: float, rho: float, n: int, d: int):
        self.c = probit_mp(mpf(pd))
        self.sr = mp.sqrt(mpf(rho))
        self.s1r = mp.sqrt(1 - mpf(rho))
        self.n, self.d = n, d

    def x(self, z: mpf) -> mpf:
        return (self.c - self.sr * z) / self.s1r

    def h(self, z: mpf) -> mpf:
        x = self.x(z)
        v = -z * z / 2
        if self.d:
            v += self.d * log_phi_mp(x)
        if self.n - self.d:
            v += (self.n - self.d) * log_phi_mp(-x)
        return v

    def dh(self, z: mpf) -> mpf:
        x = self.x(z)
        lam = lambda t: mp.npdf(t) / mp.ncdf(t)
        return -(self.sr / self.s1r) * (self.d * lam(x) - (self.n - self.d) * lam(-x)) - z

    def mode(self) -> mpf:
        # h is strictly concave, so dh is decreasing: bracket the root, then bisect/secant.
        lo, hi = mpf(-1), mpf(1)
        while self.dh(lo) < 0:
            lo *= 2
        while self.dh(hi) > 0:
            hi *= 2
        return mp.findroot(self.dh, (lo, hi), solver="anderson")

    def log_integral(self) -> mpf:
        """log of integral exp(h(z)) dz / sqrt(2 pi) = log E[p^d (1-p)^(n-d)]."""
        zs = self.mode()
        curvature = -mp.diff(self.h, zs, 2)
        sigma = 1 / mp.sqrt(curvature)
        h0 = self.h(zs)
        pts = [-mp.inf] + [zs + k * sigma for k in (-60, -30, -15, -8, -4, -2, 0, 2, 4, 8, 15, 30, 60)] + [mp.inf]
        val, err = mp.quad(lambda z: mp.exp(self.h(z) - h0), pts, error=True)
        assert err <= mpf(10) ** -35 * val, ("quad error", self.n, self.d, err)
        return h0 + mp.log(val) - mp.log(2 * mp.pi) / 2


def phi2_same_threshold(c: mpf, rho: mpf) -> mpf:
    """Phi2(c, c; rho) by Plackett's identity: dPhi2/drho = phi2(c, c; rho)."""
    density = lambda r: mp.exp(-c * c / (1 + r)) / (2 * mp.pi * mp.sqrt(1 - r * r))
    return mp.ncdf(c) ** 2 + mp.quad(density, [0, rho])


MIXTURE_CLOSED_PD = [1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.2, 0.5]
MIXTURE_CLOSED_RHO = [1e-4, 0.01, 0.05, 0.12, 0.24, 0.5, 0.9, 0.95, 0.99]
MIXTURE_SPIKES = [  # (pd, rho, n, d)
    (0.02, 0.2, 100, 0),
    (0.02, 0.2, 100, 2),
    (0.02, 0.2, 100, 15),
    (1e-3, 0.12, 10**4, 10),
    (0.05, 0.05, 10**4, 600),
    (1e-3, 0.12, 10**6, 0),
    (1e-3, 0.12, 10**6, 1000),
    (1e-3, 0.12, 10**6, 3000),
    (0.01, 0.24, 10**6, 5000),
    (0.01, 1e-4, 10**6, 10000),
    (1e-4, 0.12, 10**6, 10**6),
    # low-default portfolios: zero or few defaults among tens to hundreds of obligors
    (0.01, 0.5, 50, 0),
    (0.01, 0.9, 50, 0),
    (0.01, 0.95, 50, 0),
    (1e-3, 0.9, 10, 0),
    (1e-3, 0.9, 200, 0),
    (1e-3, 0.95, 200, 0),
    (0.05, 0.9, 50, 5),
    (0.05, 0.9, 50, 50),
    # zero defaults where few are expected (tiny PD, large n): found by the M1.7 cross-reference
    (1e-6, 0.24, 10**6, 0),
    (1e-5, 0.5, 10**6, 0),
    (1e-4, 0.35, 10**4, 0),
    (0.01, 0.5, 1000, 0),
    (0.02, 0.9, 10**5, 0),
]


def golden_mixture() -> str:
    lines = [CSV_NOTE, "kind,pd_hex,rho_hex,n,d,log_integral_hex,pd,rho,log_integral\n"]

    def add(kind: str, pd: float, rho: float, n: int, d: int, value: mpf) -> None:
        e = to_double(value)
        lines.append(row(kind, pd.hex(), rho.hex(), n, d, e.hex(), repr(pd), repr(rho), repr(e)))

    for pd in MIXTURE_CLOSED_PD:
        for rho in MIXTURE_CLOSED_RHO:
            c, r = probit_mp(mpf(pd)), mpf(rho)
            p2 = phi2_same_threshold(c, r)
            add("closed", pd, rho, 1, 1, mp.log(mpf(pd)))  # E[p] = PD
            add("closed", pd, rho, 1, 0, mp.log1p(-mpf(pd)))  # E[1-p] = 1 - PD
            add("closed", pd, rho, 2, 2, mp.log(p2))  # E[p^2] = Phi2(c, c; rho)
            add("closed", pd, rho, 2, 1, mp.log(mpf(pd) - p2))  # E[p(1-p)] = PD - Phi2
    for pd, rho, n, d in MIXTURE_SPIKES:
        add("spike", pd, rho, n, d, Mixture(pd, rho, n, d).log_integral())
    return "".join(lines)


def self_checks() -> None:
    # Plackett's identity and the generic mixture integral must agree where both apply.
    for pd, rho in ((0.01, 0.12), (1e-4, 0.5)):
        c, r = probit_mp(mpf(pd)), mpf(rho)
        a = mp.log(phi2_same_threshold(c, r))
        b = Mixture(pd, rho, 2, 2).log_integral()
        assert abs(a - b) <= mpf(10) ** -30 * abs(a), ("Plackett vs mixture integral", pd, rho)
        assert abs(Mixture(pd, rho, 1, 1).log_integral() - mp.log(mpf(pd))) <= mpf(10) ** -30, ("E[p]", pd, rho)


# --- driver -------------------------------------------------------------------------


def build_outputs() -> dict[str, str]:
    self_checks()
    rules, smallest = {}, {}
    for n in GH_SIZES:
        z, log_w = standard_normal_rule(n)
        check_rule(n, z, log_w)
        rules[n] = ([to_double(v) for v in z], [to_double(v) for v in log_w])
        smallest[n] = float(mp.exp(min(log_w)))
    return {
        f"{GEN_DIR}/gauss_hermite.hpp": render_gauss_hermite(rules, smallest),
        f"{GEN_DIR}/gauss_legendre.hpp": render_gauss_legendre(
            {n: tuple([to_double(v) for v in part] for part in legendre_rule(n)) for n in GL_SIZES}
        ),
        f"{GOLDEN_DIR}/gh_moments.csv": golden_moments(),
        f"{GOLDEN_DIR}/mixture.csv": golden_mixture(),
    }


if __name__ == "__main__":
    sys.exit(run(GENERATOR, MANIFEST, build_outputs, __doc__))
