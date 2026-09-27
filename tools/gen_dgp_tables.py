#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Generate the DGP constants, accuracy tables and reference panels (M1.6).

Two independent jobs (D-052..D-057, D-061, D-109, D-110):

1. mpmath: the constants in dgp/generated/constants.hpp, and reference values for measuring
   the accuracy of the DGP's in-house log, exp and Phi (Phi^-1 reuses special/probit.csv).

2. A pure-Python mirror of dgp/: Philox4x32-10, det_log, det_exp, det_ncdf, det_probit, the
   uniform conversion, the Marsaglia polar factor draw and the Bernoulli sum. It uses only
   operations that are exact or correctly rounded in IEEE binary64, so they match C++ bit for
   bit: + - * /, math.sqrt, math.floor, math.ldexp, and int -> float for integers below 2^53.
   It never calls math.log, math.exp or math.erfc, which go through the platform libm. Every
   function has the same operation order as the C++. The mirror must reproduce Random123's
   known answers. It produces the reference panels and their SHA-256, which every platform's
   C++ build must match bit for bit (the cross-OS check, D-057).

3. The same mirror for resampling draws (resample/weights.hpp, M2b, D-132, D-133): uniforms and
   indices in the resampling key domain (seed XOR "RSMPBOOT"), for the iid and moving-block
   bootstrap, which the C++ must reproduce bit for bit.

Run from the repository root:
    uv run --no-project --with mpmath==1.4.1 python tools/gen_dgp_tables.py
    uv run --no-project --with mpmath==1.4.1 python tools/gen_dgp_tables.py --check
"""

from __future__ import annotations

import hashlib
import math
import struct
import sys

from mpmath import mp, mpf
from tablegen import ROOT, csv_note, cxx, header_banner, pow10, row, run, to_double

GENERATOR = "tools/gen_dgp_tables.py"
GEN_DIR = "dgp/generated"
GOLDEN_DIR = "tests/golden/dgp"
MANIFEST = f"{GOLDEN_DIR}/MANIFEST.json"
CSV_NOTE = csv_note(GENERATOR)
KAT_FILE = f"{GOLDEN_DIR}/random123_philox4x32_kat.txt"

INF = math.inf
NAN = math.nan
DBL_MIN = 2.0**-1022
TWO54 = math.ldexp(1.0, 54)
TWO_M52 = math.ldexp(1.0, -52)
MILLS_ITERATIONS = 128  # fixed counts, worst case measured: 106 steps at t = 2 (D-113)
SERIES_TERMS = 32  # worst case measured: 23 terms at |x| = 2

LOG_TERMS = 10  # log m = 2s + 2s z sum_{k=1..10} z^(k-1)/(2k+1), z = s^2 <= 0.0295
EXP_TERMS = 16  # exp r = sum_{n=0..16} r^n / n!, |r| <= 0.347
NCDF_TAIL = 2.0  # continued fraction beyond |x| = 2, series inside
FACTOR_BLOCKS = 1 << 16  # Philox blocks 0 .. 2^16-1 of each period: the factor draw (D-109)
RESAMPLE_KEY_TAG = 0x52534D50424F4F54  # "RSMPBOOT": the resampling key domain (D-132)
SCHEME_IID, SCHEME_BLOCK = 1, 2
RESAMPLE_PERIODS, RESAMPLE_REPLICATES, RESAMPLE_BLOCK = 23, 4, 5  # odd T exercises both word pairs


# --- mpmath constants -----------------------------------------------------------------------


def ln2_split() -> tuple[float, float]:
    """ln 2 = hi + lo, hi with 32 significant bits so k * hi is exact for |k| < 2^21."""
    ln2 = mp.log(2)
    bits = struct.unpack("<Q", struct.pack("<d", to_double(ln2)))[0]
    hi = struct.unpack("<d", struct.pack("<Q", bits & ~((1 << 21) - 1)))[0]
    return hi, to_double(ln2 - mpf(hi))


def constants() -> dict[str, float]:
    hi, lo = ln2_split()
    return {
        "kLn2Hi": hi,
        "kLn2Lo": lo,
        "kInvLn2": to_double(1 / mp.log(2)),
        "kSqrt2": to_double(mp.sqrt(2)),
        "kInvSqrt2Pi": to_double(1 / mp.sqrt(2 * mp.pi)),
    }


# --- the pure-Python mirror (keep in step with dgp/*.cpp, operation for operation) ------------

MASK32 = 0xFFFFFFFF
PHILOX_M0, PHILOX_M1 = 0xD2511F53, 0xCD9E8D57
PHILOX_W0, PHILOX_W1 = 0x9E3779B9, 0xBB67AE85


def philox4x32(ctr: tuple[int, int, int, int], key: tuple[int, int], rounds: int = 10) -> tuple[int, int, int, int]:
    c0, c1, c2, c3 = ctr
    k0, k1 = key
    for r in range(rounds):
        if r:
            k0 = (k0 + PHILOX_W0) & MASK32
            k1 = (k1 + PHILOX_W1) & MASK32
        p0 = PHILOX_M0 * c0
        p1 = PHILOX_M1 * c2
        c0, c1, c2, c3 = (p1 >> 32) ^ c1 ^ k0, p1 & MASK32, (p0 >> 32) ^ c3 ^ k1, p0 & MASK32
    return c0, c1, c2, c3


class Mirror:
    def __init__(self, c: dict[str, float]):
        self.ln2_hi, self.ln2_lo = c["kLn2Hi"], c["kLn2Lo"]
        self.inv_ln2, self.sqrt2, self.inv_sqrt_2pi = c["kInvLn2"], c["kSqrt2"], c["kInvSqrt2Pi"]
        self.log_coef = [0.0] + [1.0 / float(2 * k + 1) for k in range(1, LOG_TERMS + 1)]
        self.exp_coef = [1.0]
        for n in range(1, EXP_TERMS + 1):
            self.exp_coef.append(self.exp_coef[-1] / float(n))

    def log(self, x: float) -> float:
        if math.isnan(x):
            return x
        if x < 0.0:
            return NAN
        if x == 0.0:
            return -INF
        if x == INF:
            return INF
        e_adj = 0
        if x < DBL_MIN:
            x = x * TWO54
            e_adj = -54
        bits = struct.unpack("<Q", struct.pack("<d", x))[0]
        e = ((bits >> 52) & 0x7FF) - 1023 + e_adj
        m = struct.unpack("<d", struct.pack("<Q", (bits & 0x000FFFFFFFFFFFFF) | 0x3FF0000000000000))[0]
        if m > self.sqrt2:
            m = m * 0.5
            e = e + 1
        f = m - 1.0
        s = f / (2.0 + f)
        z = s * s
        poly = self.log_coef[LOG_TERMS]
        for k in range(LOG_TERMS - 1, 0, -1):
            poly = poly * z + self.log_coef[k]
        r = z * poly
        log_m = 2.0 * s + 2.0 * (s * r)
        fe = float(e)
        return (fe * self.ln2_lo + log_m) + fe * self.ln2_hi

    def exp(self, x: float) -> float:
        if math.isnan(x):
            return x
        if x > 800.0:
            return INF
        if x < -800.0:
            return 0.0
        k = math.floor(x * self.inv_ln2 + 0.5)
        fk = float(k)
        r = (x - fk * self.ln2_hi) - fk * self.ln2_lo
        p = self.exp_coef[EXP_TERMS]
        for n in range(EXP_TERMS - 1, -1, -1):
            p = p * r + self.exp_coef[n]
        return math.ldexp(p, k)

    def phi(self, x: float) -> float:
        return self.inv_sqrt_2pi * self.exp(-0.5 * (x * x))

    @staticmethod
    def mills(t: float) -> float:
        tiny = 1e-300
        f = t
        c = t
        d = 0.0
        for j in range(1, MILLS_ITERATIONS + 1):
            fj = float(j)
            d = t + fj * d
            c = t + fj / c
            if d == 0.0:
                d = tiny
            if c == 0.0:
                c = tiny
            d = 1.0 / d
            delta = c * d
            f = f * delta
        return 1.0 / f

    def ncdf(self, x: float) -> float:
        if math.isnan(x):
            return x
        if x == -INF:
            return 0.0
        if x == INF:
            return 1.0
        if x < -NCDF_TAIL:
            return self.phi(x) * self.mills(-x)
        if x > NCDF_TAIL:
            return 1.0 - self.phi(x) * self.mills(x)
        x2 = x * x
        term = x
        s = x
        for n in range(1, SERIES_TERMS + 1):
            term = (term * x2) / float(2 * n + 1)
            s = s + term
        return 0.5 + self.phi(x) * s

    def probit(self, p: float) -> float:
        """Wichura AS241 (PPND16), same coefficients and Horner order as dgp/det_math.cpp."""
        if not (0.0 <= p <= 1.0):
            return NAN
        if p == 0.0:
            return -INF
        if p == 1.0:
            return INF
        q = p - 0.5
        if abs(q) <= 0.425:
            r = 0.180625 - q * q
            num = (((((((2.5090809287301226727e3 * r + 3.3430575583588128105e4) * r + 6.7265770927008700853e4) * r
                       + 4.5921953931549871457e4) * r + 1.3731693765509461125e4) * r + 1.9715909503065514427e3) * r
                    + 1.3314166789178437745e2) * r + 3.3871328727963666080e0)
            den = (((((((5.2264952788528545610e3 * r + 2.8729085735721942674e4) * r + 3.9307895800092710610e4) * r
                       + 2.1213794301586595867e4) * r + 5.3941960214247511077e3) * r + 6.8718700749205790830e2) * r
                    + 4.2313330701600911252e1) * r + 1.0)
            return q * num / den
        r = math.sqrt(-self.log(p if q < 0.0 else 1.0 - p))
        if r <= 5.0:
            r = r - 1.6
            x = ((((((((7.74545014278341407640e-4 * r + 2.27238449892691845833e-2) * r + 2.41780725177450611770e-1) * r
                       + 1.27045825245236838258e0) * r + 3.64784832476320460504e0) * r + 5.76949722146069140550e0) * r
                     + 4.63033784615654529590e0) * r + 1.42343711074968357734e0)
                 / (((((((1.05075007164441684324e-9 * r + 5.47593808499534494600e-4) * r + 1.51986665636164571966e-2) * r
                         + 1.48103976427480074590e-1) * r + 6.89767334985100004550e-1) * r + 1.67638483018380384940e0) * r
                       + 2.05319162663775882187e0) * r + 1.0))
        else:
            r = r - 5.0
            x = ((((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r + 1.24266094738807843860e-3) * r
                       + 2.65321895265761230930e-2) * r + 2.96560571828504891230e-1) * r + 1.78482653991729133580e0) * r
                     + 5.46378491116411436990e0) * r + 6.65790464350110377720e0)
                 / (((((((2.04426310338993978564e-15 * r + 1.42151175831644588870e-7) * r + 1.84631831751005468180e-5) * r
                         + 7.86869131145613259100e-4) * r + 1.48753612908506148525e-2) * r + 1.36929880922735805310e-1) * r
                       + 5.99832206555887937690e-1) * r + 1.0))
        return -x if q < 0.0 else x

    @staticmethod
    def uniform(a: int, b: int) -> float:
        """(k + 1/2) * 2^-52 with k = top 26 bits of a, then top 26 bits of b: in (0, 1), exact."""
        k = ((a >> 6) << 26) | (b >> 6)
        return (float(k) + 0.5) * TWO_M52

    def factor(self, key: tuple[int, int], scenario: int, replicate: int, period: int) -> float:
        for block in range(FACTOR_BLOCKS):
            w = philox4x32((scenario, replicate, period, block), key)
            v1 = 2.0 * self.uniform(w[0], w[1]) - 1.0
            v2 = 2.0 * self.uniform(w[2], w[3]) - 1.0
            s = v1 * v1 + v2 * v2
            if s > 0.0 and s < 1.0:
                return v1 * math.sqrt((-2.0 * self.log(s)) / s)
        raise RuntimeError("factor stream exhausted")

    def panel(self, seed: int, scenario: int, replicate: int, pd: float, rho: float, ns: list[int]):
        key = (seed & MASK32, (seed >> 32) & MASK32)
        c = self.probit(pd)
        sr = math.sqrt(rho)
        s1 = math.sqrt(1.0 - rho)
        out = []
        for t, n in enumerate(ns):
            z = self.factor(key, scenario, replicate, t)
            p = self.ncdf((c - sr * z) / s1)
            d = 0
            w = (0, 0, 0, 0)
            for i in range(n):
                if i % 2 == 0:
                    w = philox4x32((scenario, replicate, t, FACTOR_BLOCKS + i // 2), key)
                u = self.uniform(w[0], w[1]) if i % 2 == 0 else self.uniform(w[2], w[3])
                if u < p:
                    d += 1
            out.append((n, d, z))
        return out


# --- known answers ----------------------------------------------------------------------------


def check_kat() -> None:
    lines = [ln for ln in (ROOT / KAT_FILE).read_text(encoding="utf-8").splitlines() if ln.startswith("philox4x32 ")]
    assert len(lines) == 6, "expected six philox4x32 known-answer lines"
    for ln in lines:
        f = ln.split()
        rounds = int(f[1])
        h = [int(v, 16) for v in f[2:12]]
        got = philox4x32((h[0], h[1], h[2], h[3]), (h[4], h[5]), rounds)
        assert got == (h[6], h[7], h[8], h[9]), ("Random123 known answer", ln)


# --- outputs ------------------------------------------------------------------------------------

REFERENCE_SEED = 0x5641_5349_4345_4B31  # arbitrary fixed 64-bit seed
REFERENCE_PANELS = [  # (panel, scenario, replicate, pd, rho, n per period)
    (0, 0, 0, 0.01, 0.12, [10, 50, 100, 500, 1000, 2000, 5000, 1000]),
    (1, 0, 1, 0.01, 0.12, [10, 50, 100, 500, 1000, 2000, 5000, 1000]),
    (2, 0, 7, 0.01, 0.12, [10, 50, 100, 500, 1000, 2000, 5000, 1000]),
    (3, 1, 0, 1e-4, 0.3, [1000, 1000, 1000, 1000]),
    (4, 2, 0, 0.2, 0.05, [100, 100, 100]),
    (5, 3, 0, 0.5, 0.9, [0, 1, 2, 3]),
]


def render_constants(c: dict[str, float]) -> str:
    out = [
        header_banner(GENERATOR, MANIFEST),
        (
            "//\n// Constants for the DGP's own deterministic maths (dgp/det_math.cpp). kLn2Hi keeps 32\n"
            "// significant bits, so k * kLn2Hi is exact for |k| < 2^21 (Cody-Waite reduction).\n"
            "#pragma once\n\nnamespace vcal::dgp::generated {\n\n"
        ),
    ]
    for name, v in c.items():
        out.append(f"inline constexpr double {name} = {cxx(v)};  // {v!r}\n")
    out.append("\n}  // namespace vcal::dgp::generated\n")
    return "".join(out)


def golden_accuracy() -> dict[str, str]:
    xs_log = {pow10(mpf(k) / 2) for k in range(-614, 617)}  # 1e-307 .. 1e308
    xs_log |= {1.0 + i / 512 for i in range(-256, 513)}  # [0.5, 2] finely
    xs_log |= {1.0 + 2.0**-k for k in range(1, 53)} | {1.0 - 2.0**-k for k in range(1, 54)}
    xs_log |= {5e-324, 1e-310, 2.2250738585072014e-308}
    log_rows = [CSV_NOTE, "x_hex,expected_hex,x,expected\n"]
    for x in sorted(xs_log):
        e = to_double(mp.log(mpf(x)))
        log_rows.append(row(x.hex(), e.hex(), repr(x), repr(e)))

    xs_exp = {to_double(mpf(i) / 8) for i in range(-5960, 5678)}  # [-745, 709.75]
    xs_exp |= {pow10(-mpf(k) / 2) for k in range(2, 60)} | {-pow10(-mpf(k) / 2) for k in range(2, 60)}
    exp_rows = [CSV_NOTE, "x_hex,expected_hex,x,expected\n"]
    for x in sorted(xs_exp):
        e = to_double(mp.exp(mpf(x)))
        exp_rows.append(row(x.hex(), e.hex(), repr(x), repr(e)))

    xs_ncdf = {to_double(mpf(i) / 100) for i in range(-3800, 801)}  # [-38, 8]
    ncdf_rows = [CSV_NOTE, "x_hex,expected_hex,x,expected\n"]
    for x in sorted(xs_ncdf):
        e = to_double(mp.ncdf(mpf(x)))
        ncdf_rows.append(row(x.hex(), e.hex(), repr(x), repr(e)))
    return {
        f"{GOLDEN_DIR}/log.csv": "".join(log_rows),
        f"{GOLDEN_DIR}/exp.csv": "".join(exp_rows),
        f"{GOLDEN_DIR}/ncdf.csv": "".join(ncdf_rows),
    }


def golden_mirror(m: Mirror) -> dict[str, str]:
    """Values the C++ must reproduce bit for bit (not accuracy: identity with the mirror)."""
    xs = [5e-324, 1e-300, 1e-12, 0.1, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 10.0, 1e10, 1e300]
    rows = [CSV_NOTE, "function,x_hex,value_hex,x,value\n"]
    for x in xs:
        rows.append(row("log", x.hex(), m.log(x).hex(), repr(x), repr(m.log(x))))
    for x in [-745.0, -700.0, -20.5, -1.0, -1e-12, 0.0, 1e-12, 0.3466, 1.0, 20.5, 700.0, 709.7]:
        rows.append(row("exp", x.hex(), m.exp(x).hex(), repr(x), repr(m.exp(x))))
    for x in [-37.5, -8.0, -3.0, -2.0000000000000004, -2.0, -1.0, 0.0, 0.5, 2.0, 2.0000000000000004, 5.0, 8.5]:
        rows.append(row("ncdf", x.hex(), m.ncdf(x).hex(), repr(x), repr(m.ncdf(x))))
    for x in [5e-324, 1e-300, 1e-10, 0.01, 0.075, 0.3, 0.5, 0.7, 0.925, 0.99, 0.999999]:
        rows.append(row("probit", x.hex(), m.probit(x).hex(), repr(x), repr(m.probit(x))))
    for a, b in [(0, 0), (MASK32, MASK32), (0x12345678, 0x9ABCDEF0)]:
        u = m.uniform(a, b)
        rows.append(row("uniform", f"{a:08x}{b:08x}", u.hex(), f"{a:#x}/{b:#x}", repr(u)))
    return {f"{GOLDEN_DIR}/mirror_values.csv": "".join(rows)}


def golden_panels(m: Mirror) -> dict[str, str]:
    rows = [CSV_NOTE, "panel,seed_hex,scenario,replicate,pd_hex,rho_hex,period,n,d,z_hex,pd,rho,z\n"]
    blob = bytearray()
    for panel, scenario, replicate, pd, rho, ns in REFERENCE_PANELS:
        for t, (n, d, z) in enumerate(m.panel(REFERENCE_SEED, scenario, replicate, pd, rho, ns)):
            rows.append(row(panel, f"{REFERENCE_SEED:016x}", scenario, replicate, pd.hex(), rho.hex(), t, n, d,
                            z.hex(), repr(pd), repr(rho), repr(z)))
            blob += struct.pack("<q", d) + struct.pack("<d", z)
    digest = hashlib.sha256(bytes(blob)).hexdigest()
    return {
        f"{GOLDEN_DIR}/reference_panels.csv": "".join(rows),
        f"{GOLDEN_DIR}/reference_panels.sha256": digest + "\n",
    }


def golden_sha256_vectors() -> dict[str, str]:
    msgs = [b"", b"abc", b"a" * 55, b"a" * 56, b"a" * 63, b"a" * 64, b"a" * 65, bytes(range(256)) * 4]
    rows = [CSV_NOTE, "message_hex,sha256_hex\n"]
    for msg in msgs:
        rows.append(row(msg.hex(), hashlib.sha256(msg).hexdigest()))
    return {f"{GOLDEN_DIR}/sha256_vectors.csv": "".join(rows)}


# --- resampling draws (D-132, D-133) ---------------------------------------------------------


def resample_uniform(seed: int, scheme: int, replicate: int, j: int) -> float:
    k = seed ^ RESAMPLE_KEY_TAG
    w = philox4x32((scheme, replicate, j // 2, 0), (k & 0xFFFFFFFF, (k >> 32) & 0xFFFFFFFF))
    return Mirror.uniform(w[0], w[1]) if j % 2 == 0 else Mirror.uniform(w[2], w[3])


def resample_index(seed: int, scheme: int, replicate: int, j: int, m: int) -> int:
    i = math.floor(resample_uniform(seed, scheme, replicate, j) * float(m))
    return min(i, m - 1)


def golden_resample() -> dict[str, str]:
    """Every draw of the reference iid and moving-block bootstrap: its uniform, and the period
    index (iid) or block start (moving block) it becomes; plus the resulting index rows."""
    seed, T, B, ell = REFERENCE_SEED, RESAMPLE_PERIODS, RESAMPLE_REPLICATES, RESAMPLE_BLOCK
    rows = [CSV_NOTE, f"# seed {seed:016x}, T = {T}, B = {B}, block length {ell} (moving block)\n",
            "scheme,replicate,draw,u_hex,value,position,index\n"]
    for b in range(B):
        for j in range(T):  # iid: draw j is position j
            u = resample_uniform(seed, SCHEME_IID, b, j)
            idx = resample_index(seed, SCHEME_IID, b, j, T)
            rows.append(row("iid", b, j, u.hex(), idx, j, idx))
        pos, i = 0, 0
        while pos < T:  # moving block: draw i is block i's start
            u = resample_uniform(seed, SCHEME_BLOCK, b, i)
            start = resample_index(seed, SCHEME_BLOCK, b, i, T - ell + 1)
            for s_ in range(ell):
                if pos >= T:
                    break
                rows.append(row("block", b, i, u.hex(), start, pos, start + s_))
                pos += 1
            i += 1
    return {f"{GOLDEN_DIR}/resample_reference.csv": "".join(rows)}


def build_outputs() -> dict[str, str]:
    check_kat()
    c = constants()
    m = Mirror(c)
    out = {f"{GEN_DIR}/constants.hpp": render_constants(c)}
    out.update(golden_accuracy())
    out.update(golden_mirror(m))
    out.update(golden_panels(m))
    out.update(golden_sha256_vectors())
    out.update(golden_resample())
    return out


if __name__ == "__main__":
    sys.exit(run(GENERATOR, MANIFEST, build_outputs, __doc__))
