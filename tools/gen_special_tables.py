#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Generate the special-function coefficient headers and golden tables.

Every constant in core/special/generated/ and every expected value in
tests/golden/special/ comes from this script (D-033, D-067); nothing there is
typed by hand. Values are computed in mpmath at 50 significant digits and
rounded once to IEEE binary64. Files store them as exact hex floats (decimal
columns in the CSVs are for readers only). Shared machinery: tools/tablegen.py.

Run from the repository root:
    uv run --no-project --with mpmath==1.4.1 python tools/gen_special_tables.py
    uv run --no-project --with mpmath==1.4.1 python tools/gen_special_tables.py --check
"""

from __future__ import annotations

import math
import sys

from mpmath import mp, mpf
from tablegen import (
    csv_note,
    cxx,
    erfcx_mp,
    header_banner,
    log_phi_mp,
    pow10,
    probit_mp,
    row,
    run,
    to_double,
)

GENERATOR = "tools/gen_special_tables.py"
GEN_DIR = "core/special/generated"
GOLDEN_DIR = "tests/golden/special"
MANIFEST = f"{GOLDEN_DIR}/MANIFEST.json"
HEADER_BANNER = header_banner(GENERATOR, MANIFEST)
CSV_NOTE = csv_note(GENERATOR)

# erfcx on x >= 0: a Chebyshev series on each [a, b] in t = (2x - (a+b)) / (b-a),
# then x*erfcx(x) as a series in t = 2*(X0/x) - 1 for x > X0.
ERFCX_SEGMENTS = [(1, 0, 1), (2, 1, 3), (3, 3, 8)]  # (id, a, b)
ERFCX_TAIL_X0 = 8
CHEB_NODES = 60
CHEB_DROP = mpf(2) ** -58  # drop coefficients below this fraction of the largest

# Stirling error delta(n) = log n! - [(n + 1/2) log n - n + log sqrt(2 pi)] (Loader 2000).
STIRLING_TABLE_MAX = 30
STIRLING_SERIES_TERMS = 5


# --- reference functions (mpmath) ----------------------------------------------------


def lbinom_mp(n: int, k: int) -> mpf:
    if k == 0 or k == n:
        return mpf(0)
    return mp.loggamma(n + 1) - mp.loggamma(k + 1) - mp.loggamma(n - k + 1)


def stirling_error_mp(n: int) -> mpf:
    m = mpf(n)
    return mp.loggamma(m + 1) - ((m + mpf(1) / 2) * mp.log(m) - m + mp.log(mp.sqrt(2 * mp.pi)))


def inverse_mills_mp(x: mpf) -> mpf:
    return mp.npdf(x) / mp.ncdf(x)


def self_checks() -> None:
    for n, k in ((10, 3), (40, 20), (1000, 7)):
        assert abs(lbinom_mp(n, k) - mp.log(mp.binomial(n, k))) < mpf(10) ** -40, ("lbinom", n, k)


# --- Chebyshev fitting ---------------------------------------------------------------


def chebyshev_coefficients(f) -> list[float]:
    """Coefficients c_k of f(t) ~ sum c_k T_k(t) on [-1, 1], with c_0 already halved."""
    n = CHEB_NODES
    theta = [mp.pi * (j + mpf(1) / 2) / n for j in range(n)]
    vals = [f(mp.cos(t)) for t in theta]
    cs = [2 * mp.fsum(v * mp.cos(k * t) for v, t in zip(vals, theta)) / n for k in range(n)]
    cs[0] /= 2
    biggest = max(abs(c) for c in cs)
    keep = max(k for k, c in enumerate(cs) if abs(c) > CHEB_DROP * biggest) + 1
    assert keep < n - 10, "series did not converge within the node count"
    return [to_double(c) for c in cs[:keep]]


def clenshaw(c: list[float], t: float) -> float:
    """Binary64 evaluation in exactly the operation order of core/special/chebyshev.hpp."""
    b1 = b2 = 0.0
    for ck in reversed(c[1:]):
        b1, b2 = ck + 2.0 * t * b1 - b2, b1
    return c[0] + t * b1 - b2


def fit_erfcx():
    segments = []
    for sid, a, b in ERFCX_SEGMENTS:
        lo, hi = mpf(a), mpf(b)
        coeffs = chebyshev_coefficients(lambda t, lo=lo, hi=hi: erfcx_mp((hi - lo) * (t + 1) / 2 + lo))
        segments.append((sid, float(a), float(b), coeffs))

    x0 = mpf(ERFCX_TAIL_X0)

    def tail(t: mpf) -> mpf:
        w = (t + 1) / 2
        if w == 0:
            return 1 / mp.sqrt(mp.pi)
        x = x0 / w
        return x * erfcx_mp(x)

    tail_coeffs = chebyshev_coefficients(tail)

    # Measure the whole binary64 pipeline (mapping + Clenshaw), as the C++ evaluates it.
    worst = mpf(0)
    for i in range(1, 4001):
        x = i * 0.005 if i <= 2000 else ERFCX_TAIL_X0 * 1.01 ** (i - 2000)
        approx = erfcx_pipeline(segments, tail_coeffs, x)
        worst = max(worst, abs(mpf(approx) / erfcx_mp(mpf(x)) - 1))
    return segments, tail_coeffs, float(worst)


def erfcx_pipeline(segments, tail_coeffs, x: float) -> float:
    for _, a, b, c in segments:
        if x <= b:
            return clenshaw(c, (2.0 * x - (a + b)) / (b - a))
    return clenshaw(tail_coeffs, 2.0 * (ERFCX_TAIL_X0 / x) - 1.0) / x


# --- generated headers --------------------------------------------------------------


def render_constants() -> str:
    consts = [
        ("kLn2", mp.log(2), "log(2)"),
        ("kInvSqrt2", 1 / mp.sqrt(2), "1/sqrt(2)"),
        ("kLnSqrt2Pi", mp.log(mp.sqrt(2 * mp.pi)), "log(sqrt(2*pi))"),
        ("kInvSqrt2Pi", 1 / mp.sqrt(2 * mp.pi), "1/sqrt(2*pi)"),
        ("kSqrt2OverPi", mp.sqrt(2 / mp.pi), "sqrt(2/pi)"),
    ]
    out = [HEADER_BANNER, "#pragma once\n\nnamespace vcal::special::constants {\n\n"]
    for name, value, what in consts:
        v = to_double(value)
        out.append(f"inline constexpr double {name} = {cxx(v)};  // {what} = {v!r}\n")
    out.append("\n}  // namespace vcal::special::constants\n")
    return "".join(out)


def render_erfcx(segments, tail_coeffs, worst: float) -> str:
    preamble = (
        "//\n"
        "// erfcx(x) = exp(x^2) erfc(x) for x >= 0 as Chebyshev series; see core/special/erfcx.hpp.\n"
        "// Max relative error of this binary64 pipeline, measured by the generator on 4000 points\n"
        f"// in [0.005, {ERFCX_TAIL_X0}*1.01^2000]: {worst:.3e} ({worst / 2.0**-53:.2f} x 2^-53).\n"
        "// Arrays are function-local: namespace-scope arrays are not usable from device code.\n"
        "#pragma once\n\n"
        '#include "core/precision.hpp"\n'
        '#include "core/special/chebyshev.hpp"\n\n'
        "namespace vcal::special::generated {\n\n"
    )
    out = [HEADER_BANNER, preamble]

    def series(name: str, comment: str, c: list[float]) -> None:
        out.append(f"// {comment}; {len(c)} terms\n")
        out.append(f"VCAL_HD double {name}(double t) {{\n    constexpr double c[{len(c)}] = {{\n")
        out.extend(f"        {cxx(v)},\n" for v in c)
        out.append("    };\n    return chebyshev_eval(c, t);\n}\n\n")

    for sid, a, b, c in segments:
        series(f"erfcx_segment{sid}", f"x in [{a:g}, {b:g}], t = (2x - {a + b:g}) / {b - a:g}", c)
    series("erfcx_tail", f"x > {ERFCX_TAIL_X0}: x*erfcx(x) in t = 2*({ERFCX_TAIL_X0}/x) - 1", tail_coeffs)

    out.append("// Requires x >= 0 and finite.\nVCAL_HD double erfcx_nonnegative(double x) {\n")
    for sid, a, b, _ in segments:
        out.append(f"    if (x <= {b!r}) return erfcx_segment{sid}((2.0 * x - {a + b!r}) / {b - a!r});\n")
    out.append(f"    return erfcx_tail(2.0 * ({float(ERFCX_TAIL_X0)!r} / x) - 1.0) / x;\n}}\n\n")
    out.append("}  // namespace vcal::special::generated\n")
    return "".join(out)


def render_stirling() -> str:
    table = [to_double(stirling_error_mp(n)) for n in range(1, STIRLING_TABLE_MAX + 1)]
    # delta(n) ~ sum_{k>=1} B_2k / (2k (2k-1) n^(2k-1))
    coef = [mp.bernoulli(2 * k) / (2 * k * (2 * k - 1)) for k in range(1, STIRLING_SERIES_TERMS + 2)]
    series = [to_double(c) for c in coef[:-1]]
    n_first = STIRLING_TABLE_MAX + 1
    omitted = to_double(abs(coef[-1]) / mpf(n_first) ** (2 * STIRLING_SERIES_TERMS + 1))
    expr = cxx(series[-1])
    for c in reversed(series[:-1]):
        expr = f"{cxx(c)} + r2 * ({expr})"
    preamble = (
        "//\n"
        "// Stirling error delta(n) = log n! - [(n + 1/2) log n - n + log sqrt(2 pi)] (Loader 2000):\n"
        f"// tabulated for 1 <= n <= {STIRLING_TABLE_MAX}; beyond that, {STIRLING_SERIES_TERMS} terms of\n"
        "// sum_k B_2k / (2k (2k-1) n^(2k-1)). "
        f"First omitted term at n = {n_first}: {omitted:.2e} (absolute).\n"
        "#pragma once\n\n"
        '#include "core/precision.hpp"\n\n'
        "namespace vcal::special::generated {\n\n"
        f"inline constexpr int kStirlingTableMax = {STIRLING_TABLE_MAX};\n\n"
        "// Requires 1 <= n <= kStirlingTableMax.\n"
        "VCAL_HD double stirling_error_table(int n) {\n"
        f"    constexpr double t[{STIRLING_TABLE_MAX}] = {{\n"
    )
    out = [HEADER_BANNER, preamble]
    out.extend(f"        {cxx(v)},  // n = {i}\n" for i, v in enumerate(table, start=1))
    out.append("    };\n    return t[n - 1];\n}\n\n")
    out.append("// Requires n > kStirlingTableMax; r = 1/n.\nVCAL_HD double stirling_error_series(double r) {\n")
    out.append(f"    const double r2 = r * r;\n    return r * ({expr});\n}}\n\n")
    out.append("}  // namespace vcal::special::generated\n")
    return "".join(out)


# --- golden tables ------------------------------------------------------------------


def golden_erfcx() -> str:
    xs = {to_double(mpf(i) / 20) for i in range(-520, 0)}  # [-26, 0)
    xs |= {to_double(mpf(i) / 100) for i in range(3001)}  # [0, 30]
    for _, a, b in ERFCX_SEGMENTS:  # segment edges and both neighbours
        for e in (float(a), float(b)):
            xs |= {e, math.nextafter(e, -math.inf), math.nextafter(e, math.inf)}
    xs |= {pow10(mpf(i) / 4) for i in range(6, 1201)}  # 10^1.5 .. 10^300
    xs |= {pow10(-mpf(i) / 2) for i in range(2, 601)}  # 10^-1 .. 10^-300
    xs |= {-26.5, -26.6, -27.0, -30.0}  # around and past overflow of exp(x^2)
    lines = [CSV_NOTE, "x_hex,expected_hex,x,expected\n"]
    for x in sorted(xs):
        e = to_double(erfcx_mp(mpf(x)))
        lines.append(row(x.hex(), e.hex(), repr(x), repr(e)))
    return "".join(lines)


def golden_log_phi() -> str:
    xs = {to_double(mpf(i) / 20) for i in range(-800, 801)}  # [-40, 40]
    xs |= {-pow10(mpf(i) / 20) for i in range(33, 201)}  # -10^1.65 .. -10^10
    xs |= {-1e20, -1e50, -1e100, -1e150, -1e154, -1e155, -1e300}
    xs |= {to_double(mpf(i) / 100) for i in range(3700, 3951)}  # 37 .. 39.5: log Phi turns subnormal
    for k in range(1, 301):
        xs |= {10.0**-k, -(10.0**-k)}
    lines = [CSV_NOTE, "x_hex,expected_hex,x,expected\n"]
    for x in sorted(xs):
        e = to_double(log_phi_mp(mpf(x)))
        lines.append(row(x.hex(), e.hex(), repr(x), repr(e)))
    return "".join(lines)


def golden_probit() -> str:
    ps = {pow10(-mpf(i) / 4) for i in range(2, 1201)}  # 10^-0.5 .. 10^-300
    ps |= {i / 1000 for i in range(1, 1000)}
    for k in range(2, 54):
        ps |= {0.5 - 2.0**-k, 0.5 + 2.0**-k, 1.0 - 2.0**-k}
    ps |= {1.0 - 10.0**-k for k in range(1, 16)}
    ps |= {1e-310, 5e-324}
    lines = [CSV_NOTE, "p_hex,expected_hex,p,expected\n"]
    for p in sorted(ps):
        e = to_double(probit_mp(mpf(p)))
        lines.append(row(p.hex(), e.hex(), repr(p), repr(e)))
    return "".join(lines)


def golden_lbinom() -> str:
    ns = [1, 2, 3, 5, 10, 15, 16, 29, 30, 31, 50, 59, 60, 61, 66, 67, 100, 1000]
    ns += [10**4, 10**5, 10**6, 10**7, 10**9, 2**40, 2**52]
    lines = [CSV_NOTE, "n,k,expected_hex,expected\n"]
    for n in ns:
        ks = {0, 1, 2, 3, 5, 10, 100, n // 1000, n // 10, n // 3, n // 2, n - 1, n - 2, n}
        for k in sorted(k for k in ks if 0 <= k <= n):
            e = to_double(lbinom_mp(n, k))
            lines.append(row(n, k, e.hex(), repr(e)))
    return "".join(lines)


def golden_log_add_exp() -> str:
    base = [-1e300, -745.0, -30.0, -1.0, -1e-10, 0.0, 1e-10, 1.0, 2.5, 30.0, 700.0, 709.0]
    base += [to_double(mp.log(2)), to_double(mp.log(3))]
    lines = [CSV_NOTE, "a_hex,b_hex,expected_hex,a,b,expected\n"]
    for a in base:
        for b in base:
            hi, lo = max(mpf(a), mpf(b)), min(mpf(a), mpf(b))
            e = to_double(hi + mp.log1p(mp.exp(lo - hi)))
            lines.append(row(a.hex(), b.hex(), e.hex(), repr(a), repr(b), repr(e)))
    return "".join(lines)


def golden_inverse_mills() -> str:
    xs = {to_double(mpf(i) / 20) for i in range(-800, 801)}  # [-40, 40]
    xs |= {-pow10(mpf(i) / 20) for i in range(33, 201)}  # -10^1.65 .. -10^10
    for k in range(1, 301):
        xs |= {10.0**-k, -(10.0**-k)}
    lines = [CSV_NOTE, "x_hex,expected_hex,x,expected\n"]
    for x in sorted(xs):
        e = to_double(inverse_mills_mp(mpf(x)))
        lines.append(row(x.hex(), e.hex(), repr(x), repr(e)))
    return "".join(lines)


def golden_stirling() -> str:
    ns = set(range(1, 201)) | {10**k for k in range(3, 16)} | {2**k for k in range(8, 53)}
    lines = [CSV_NOTE, "n,expected_hex,expected\n"]
    for n in sorted(ns):
        e = to_double(stirling_error_mp(n))
        lines.append(row(n, e.hex(), repr(e)))
    return "".join(lines)


# --- driver -------------------------------------------------------------------------


def build_outputs() -> dict[str, str]:
    self_checks()
    segments, tail_coeffs, worst = fit_erfcx()
    return {
        f"{GEN_DIR}/constants.hpp": render_constants(),
        f"{GEN_DIR}/erfcx_chebyshev.hpp": render_erfcx(segments, tail_coeffs, worst),
        f"{GEN_DIR}/stirling_error.hpp": render_stirling(),
        f"{GOLDEN_DIR}/erfcx.csv": golden_erfcx(),
        f"{GOLDEN_DIR}/inverse_mills.csv": golden_inverse_mills(),
        f"{GOLDEN_DIR}/log_phi.csv": golden_log_phi(),
        f"{GOLDEN_DIR}/probit.csv": golden_probit(),
        f"{GOLDEN_DIR}/lbinom.csv": golden_lbinom(),
        f"{GOLDEN_DIR}/log_add_exp.csv": golden_log_add_exp(),
        f"{GOLDEN_DIR}/stirling_error.csv": golden_stirling(),
    }


if __name__ == "__main__":
    sys.exit(run(GENERATOR, MANIFEST, build_outputs, __doc__))
