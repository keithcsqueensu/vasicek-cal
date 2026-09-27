// SPDX-License-Identifier: Apache-2.0
//
// The DGP's own deterministic maths (D-052, D-053, D-061, D-110). Built only from operations
// IEEE 754 makes exact or correctly rounded (+ - * /, sqrt, floor, ldexp, bit manipulation),
// so results are identical on every platform and to the Python mirror in
// tools/gen_dgp_tables.py. Keep the two in step, operation for operation, including
// parenthesisation. The accuracy target is 1e-12 relative (D-052); the measured error against
// mpmath is in the tolerance register.
//
// Compiled with strict FP semantics and FP contraction off (D-056): a fused multiply-add would
// round differently from the mirror.
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

#include "dgp/dgp.hpp"
#include "dgp/generated/constants.hpp"

namespace vcal::dgp {
namespace {

namespace gc = generated;

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kDblMin = std::numeric_limits<double>::min();
constexpr double kTwo54 = 0x1p54;
// Fixed iteration counts (D-113), chosen for the worst case: the continued fraction needs at most
// 106 steps (at t = 2), the series at most 23 terms (at |x| = 2). No data-dependent exits.
constexpr int kMillsIterations = 128;
constexpr int kSeriesTerms = 32;
constexpr int kLogTerms = 10;
constexpr int kExpTerms = 16;
constexpr double kNcdfTail = 2.0;

// 1/(2k+1) and 1/n!, computed by IEEE division exactly as the mirror does.
struct Coefficients {
    double log_coef[kLogTerms + 1];
    double exp_coef[kExpTerms + 1];
    Coefficients() {
        log_coef[0] = 0.0;
        for (int k = 1; k <= kLogTerms; ++k) log_coef[k] = 1.0 / static_cast<double>(2 * k + 1);
        exp_coef[0] = 1.0;
        for (int n = 1; n <= kExpTerms; ++n) exp_coef[n] = exp_coef[n - 1] / static_cast<double>(n);
    }
};

const Coefficients& coef() {
    static const Coefficients c;
    return c;
}

double phi(double x) { return gc::kInvSqrt2Pi * det_exp(-0.5 * (x * x)); }

// Mills ratio R(t) = Q(t)/phi(t), t > 2: continued fraction 1/(t + 1/(t + 2/(t + ...))) by
// the modified Lentz algorithm, exactly kMillsIterations steps.
double mills(double t) {
    constexpr double tiny = 1e-300;
    double f = t;
    double c = t;
    double d = 0.0;
    for (int j = 1; j <= kMillsIterations; ++j) {
        const double fj = static_cast<double>(j);
        d = t + fj * d;
        c = t + fj / c;
        if (d == 0.0) d = tiny;
        if (c == 0.0) c = tiny;
        d = 1.0 / d;
        const double delta = c * d;
        f = f * delta;
    }
    return 1.0 / f;
}

}  // namespace

// log x = e ln 2 + log m, x = m 2^e with m in (sqrt(1/2), sqrt 2];
// log m = 2 atanh(s) = 2s + 2s z sum_{k>=1} z^(k-1)/(2k+1), s = (m-1)/(m+1), z = s^2 <= 0.0295.
double det_log(double x) {
    if (std::isnan(x)) return x;
    if (x < 0.0) return kNaN;
    if (x == 0.0) return -kInf;
    if (x == kInf) return kInf;
    int e_adj = 0;
    if (x < kDblMin) {
        x = x * kTwo54;
        e_adj = -54;
    }
    std::uint64_t bits = 0;
    std::memcpy(&bits, &x, sizeof bits);
    int e = static_cast<int>((bits >> 52) & 0x7FF) - 1023 + e_adj;
    const std::uint64_t mbits = (bits & 0x000FFFFFFFFFFFFFull) | 0x3FF0000000000000ull;
    double m = 0.0;
    std::memcpy(&m, &mbits, sizeof m);
    if (m > gc::kSqrt2) {
        m = m * 0.5;
        e = e + 1;
    }
    const double f = m - 1.0;
    const double s = f / (2.0 + f);
    const double z = s * s;
    const auto& c = coef();
    double poly = c.log_coef[kLogTerms];
    for (int k = kLogTerms - 1; k > 0; --k) poly = poly * z + c.log_coef[k];
    const double r = z * poly;
    const double log_m = 2.0 * s + 2.0 * (s * r);
    const double fe = static_cast<double>(e);
    return (fe * gc::kLn2Lo + log_m) + fe * gc::kLn2Hi;
}

// e^x = 2^k e^r, k = floor(x / ln 2 + 1/2), r = (x - k ln2_hi) - k ln2_lo (Cody-Waite),
// e^r by its Taylor series to r^16 (|r| <= 0.347), scaled by ldexp (correctly rounded).
double det_exp(double x) {
    if (std::isnan(x)) return x;
    if (x > 800.0) return kInf;
    if (x < -800.0) return 0.0;
    const double kd = std::floor(x * gc::kInvLn2 + 0.5);
    const double r = (x - kd * gc::kLn2Hi) - kd * gc::kLn2Lo;
    const auto& c = coef();
    double p = c.exp_coef[kExpTerms];
    for (int n = kExpTerms - 1; n >= 0; --n) p = p * r + c.exp_coef[n];
    return std::ldexp(p, static_cast<int>(kd));
}

// Phi(x): x < -2: phi(x) R(-x); x > 2: 1 - phi(x) R(x); otherwise the everywhere-convergent
// series Phi(x) = 1/2 + phi(x) (x + x^3/3 + x^5/15 + ...) (terms x^(2n+1)/(2n+1)!!, n <= 32).
double det_ncdf(double x) {
    if (std::isnan(x)) return x;
    if (x == -kInf) return 0.0;
    if (x == kInf) return 1.0;
    if (x < -kNcdfTail) return phi(x) * mills(-x);
    if (x > kNcdfTail) return 1.0 - phi(x) * mills(x);
    const double x2 = x * x;
    double term = x;
    double s = x;
    for (int n = 1; n <= kSeriesTerms; ++n) {
        term = (term * x2) / static_cast<double>(2 * n + 1);
        s = s + term;
    }
    return 0.5 + phi(x) * s;
}

// Wichura (1988) AS241 PPND16, on det_log and sqrt.
double det_probit(double p) {
    if (!(p >= 0.0 && p <= 1.0)) return kNaN;
    if (p == 0.0) return -kInf;
    if (p == 1.0) return kInf;
    const double q = p - 0.5;
    if (std::fabs(q) <= 0.425) {
        const double r = 0.180625 - q * q;
        const double num = (((((((2.5090809287301226727e3 * r + 3.3430575583588128105e4) * r + 6.7265770927008700853e4) * r +
                                4.5921953931549871457e4) * r + 1.3731693765509461125e4) * r + 1.9715909503065514427e3) * r +
                              1.3314166789178437745e2) * r + 3.3871328727963666080e0);
        const double den = (((((((5.2264952788528545610e3 * r + 2.8729085735721942674e4) * r + 3.9307895800092710610e4) * r +
                                2.1213794301586595867e4) * r + 5.3941960214247511077e3) * r + 6.8718700749205790830e2) * r +
                              4.2313330701600911252e1) * r + 1.0);
        return q * num / den;
    }
    double r = std::sqrt(-det_log(q < 0.0 ? p : 1.0 - p));
    double x;
    if (r <= 5.0) {
        r = r - 1.6;
        x = ((((((((7.74545014278341407640e-4 * r + 2.27238449892691845833e-2) * r + 2.41780725177450611770e-1) * r +
                  1.27045825245236838258e0) * r + 3.64784832476320460504e0) * r + 5.76949722146069140550e0) * r +
                4.63033784615654529590e0) * r + 1.42343711074968357734e0) /
             (((((((1.05075007164441684324e-9 * r + 5.47593808499534494600e-4) * r + 1.51986665636164571966e-2) * r +
                   1.48103976427480074590e-1) * r + 6.89767334985100004550e-1) * r + 1.67638483018380384940e0) * r +
                 2.05319162663775882187e0) * r + 1.0));
    } else {
        r = r - 5.0;
        x = ((((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r + 1.24266094738807843860e-3) * r +
                  2.65321895265761230930e-2) * r + 2.96560571828504891230e-1) * r + 1.78482653991729133580e0) * r +
                5.46378491116411436990e0) * r + 6.65790464350110377720e0) /
             (((((((2.04426310338993978564e-15 * r + 1.42151175831644588870e-7) * r + 1.84631831751005468180e-5) * r +
                   7.86869131145613259100e-4) * r + 1.48753612908506148525e-2) * r + 1.36929880922735805310e-1) * r +
                 5.99832206555887937690e-1) * r + 1.0));
    }
    return q < 0.0 ? -x : x;
}

double uniform_from_words(std::uint32_t a, std::uint32_t b) {
    const std::uint64_t k = (static_cast<std::uint64_t>(a >> 6) << 26) | static_cast<std::uint64_t>(b >> 6);
    return (static_cast<double>(k) + 0.5) * 0x1p-52;
}

bool fp_environment_ok() {
    // Flush-to-zero turns DBL_MIN / 2 into 0; denormals-are-zero turns a subnormal input into 0.
    volatile double tiny = kDblMin;
    volatile double half = tiny / 2.0;
    volatile double back = half * 2.0;
    return half != 0.0 && back == kDblMin;
}

}  // namespace vcal::dgp
