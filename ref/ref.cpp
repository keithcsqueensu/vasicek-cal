// SPDX-License-Identifier: Apache-2.0
#include "ref/ref.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

namespace vcalref {
namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kEps = std::numeric_limits<double>::epsilon();
const double kPi = std::acos(-1.0);
const double kLogSqrt2Pi = 0.5 * std::log(2.0 * kPi);
const double kSqrt2 = std::sqrt(2.0);

// Compensated (Neumaier) running sum.
class Sum {
public:
    void add(double x) {
        const double s = sum_ + x;
        comp_ += std::fabs(sum_) >= std::fabs(x) ? (sum_ - s) + x : (x - s) + sum_;
        sum_ = s;
    }
    double value() const { return sum_ + comp_; }

private:
    double sum_ = 0.0;
    double comp_ = 0.0;
};

// Mills ratio R(t) = Q(t) / phi(t), t > 2, from the continued fraction
//     R(t) = 1 / (t + 1/(t + 2/(t + 3/(t + ...))))
// evaluated with the modified Lentz algorithm.
double mills_ratio(double t) {
    constexpr double tiny = 1e-300;
    double f = t;
    double c = t;
    double d = 0.0;
    for (int j = 1; j < 10000; ++j) {
        d = t + j * d;
        c = t + j / c;
        if (d == 0.0) d = tiny;
        if (c == 0.0) c = tiny;
        d = 1.0 / d;
        const double delta = c * d;
        f *= delta;
        if (std::fabs(delta - 1.0) <= kEps / 2) break;
    }
    return 1.0 / f;
}

// Phi^-1(p) for 0 < p <= 1/2 by Newton on log Phi(x) = log p, safeguarded by the bracket
// [-40, 0] (log Phi(-40) ~ -804 < log(5e-324) ~ -744) with bisection whenever Newton leaves it.
double ncdf_inv_lower(double p) {
    if (p == 0.5) return 0.0;
    const double target = std::log(p);
    double lo = -40.0;
    double hi = 0.0;
    double x = std::max(-std::sqrt(-2.0 * target), -39.0);
    for (int it = 0; it < 500; ++it) {
        const double lp = log_ncdf(x);
        const double f = lp - target;
        if (f == 0.0) return x;
        if (f > 0.0) {
            hi = x;
        } else {
            lo = x;
        }
        const double slope = std::exp(-0.5 * x * x - kLogSqrt2Pi - lp);  // d/dx log Phi = phi / Phi
        double next = x - f / slope;
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        const bool done = std::fabs(next - x) <= 2.0 * kEps * std::fabs(next) || hi - lo <= 2.0 * kEps * std::fabs(lo);
        x = next;
        if (done) break;
    }
    return x;
}

// Phi^-1(1/2 + q) for |q| <= 1/4 by Newton on erf(x / sqrt 2) / 2 = q. Near the centre the
// log-space solver's absolute noise (~1e-16) would swamp tiny results; erf keeps relative accuracy.
double ncdf_inv_central(double q) {
    if (q == 0.0) return 0.0;
    double x = q * std::sqrt(2.0 * kPi);
    for (int it = 0; it < 100; ++it) {
        const double g = 0.5 * std::erf(x / kSqrt2) - q;
        const double step = g / (std::exp(-0.5 * x * x) / std::sqrt(2.0 * kPi));
        x -= step;
        if (std::fabs(step) <= kEps * std::fabs(x)) break;
    }
    return x;
}

// --- Gauss-Legendre rules, nodes by Newton on the Legendre recurrence (computed, not typed) --

struct GaussLegendre {
    std::vector<double> node;
    std::vector<double> weight;
};

GaussLegendre make_gauss_legendre(int n) {
    GaussLegendre r{std::vector<double>(static_cast<std::size_t>(n)), std::vector<double>(static_cast<std::size_t>(n))};
    for (int i = 0; i < (n + 1) / 2; ++i) {
        double z = std::cos(kPi * (i + 0.75) / (n + 0.5));
        double dp = 0.0;
        for (int it = 0; it < 100; ++it) {
            double p1 = 1.0;
            double p2 = 0.0;
            for (int j = 1; j <= n; ++j) {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            dp = n * (z * p1 - p2) / (z * z - 1.0);
            const double z_prev = z;
            z = z_prev - p1 / dp;
            if (std::fabs(z - z_prev) <= kEps) break;
        }
        const auto lo = static_cast<std::size_t>(i);
        const auto hi = static_cast<std::size_t>(n - 1 - i);
        r.node[lo] = -z;
        r.node[hi] = z;
        r.weight[lo] = r.weight[hi] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
    return r;
}

const GaussLegendre& gl10() {
    static const GaussLegendre r = make_gauss_legendre(10);
    return r;
}
const GaussLegendre& gl20() {
    static const GaussLegendre r = make_gauss_legendre(20);
    return r;
}

double apply(const GaussLegendre& rule, const std::function<double(double)>& f, double a, double b) {
    const double mid = 0.5 * (a + b);
    const double half = 0.5 * (b - a);
    Sum s;
    for (std::size_t i = 0; i < rule.node.size(); ++i) s.add(rule.weight[i] * f(mid + half * rule.node[i]));
    return half * s.value();
}

// Adaptive bisection: accept an interval when the 10- and 20-point rules agree within tol, or
// within the integrand's own relative noise floor times the interval's value. A tolerance
// below that noise could never be met and would split to the depth cap everywhere.
double adaptive(const std::function<double(double)>& f, double a, double b, double tol, double noise_rel,
                int depth) {
    const double q20 = apply(gl20(), f, a, b);
    const double q10 = apply(gl10(), f, a, b);
    const double diff = std::fabs(q20 - q10);
    if (diff <= tol || diff <= noise_rel * std::fabs(q20) || depth >= 30) return q20;
    const double m = 0.5 * (a + b);
    return adaptive(f, a, m, 0.5 * tol, noise_rel, depth + 1) + adaptive(f, m, b, 0.5 * tol, noise_rel, depth + 1);
}

// --- one-dimensional maximisation ------------------------------------------------------------

struct Max1 {
    double x;
    double value;
};

// Golden-section search for the maximum of a unimodal f on [a, b], to |interval| <= tol.
Max1 golden_max(const std::function<double(double)>& f, double a, double b, double tol) {
    const double r = (std::sqrt(5.0) - 1.0) / 2.0;
    double x1 = b - r * (b - a);
    double x2 = a + r * (b - a);
    double f1 = f(x1);
    double f2 = f(x2);
    while (b - a > tol) {
        if (f1 >= f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - r * (b - a);
            f1 = f(x1);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + r * (b - a);
            f2 = f(x2);
        }
    }
    return f1 >= f2 ? Max1{x1, f1} : Max1{x2, f2};
}

double logit(double v) { return std::log(v) - std::log1p(-v); }
double inv_logit(double u) {
    const double e = std::exp(-std::fabs(u));
    return u >= 0.0 ? 1.0 / (1.0 + e) : e / (1.0 + e);
}

}  // namespace

// Regions: x <= -20 continued fraction; -20 < x <= 2 std::erfc; x > 2 log1p(-Q) with the upper
// tail Q(x) = phi(x) R(x) from the continued fraction. Beyond x = 2 the rounding of x / sqrt 2
// would be amplified ~2 x^2 times by erfc, so there phi(x) uses the exact square
// x^2 = hi + lo (lo = fma(x, x, -hi), an error-free product).
double log_ncdf(double x) {
    if (std::isnan(x)) return x;
    if (x == -kInf) return -kInf;
    if (x > 2.0) {
        if (x == kInf) return 0.0;
        const double hi = x * x;
        const double lo = std::fma(x, x, -hi);
        const double phi = std::exp(-0.5 * hi) * std::exp(-0.5 * lo) / std::sqrt(2.0 * kPi);
        return std::log1p(-phi * mills_ratio(x));
    }
    if (x > 0.0) return std::log1p(-0.5 * std::erfc(x / kSqrt2));
    if (x > -20.0) return std::log(0.5 * std::erfc(-x / kSqrt2));
    const double t = -x;
    return -0.5 * t * t - kLogSqrt2Pi + std::log(mills_ratio(t));
}

double ncdf_inv(double p) {
    if (!(p >= 0.0 && p <= 1.0)) return kNaN;
    if (p == 0.0) return -kInf;
    if (p == 1.0) return kInf;
    if (p >= 0.25 && p <= 0.75) return ncdf_inv_central(p - 0.5);  // p - 1/2 is exact here
    if (p > 0.5) return -ncdf_inv_lower(1.0 - p);                  // 1 - p is exact for p >= 1/2
    return ncdf_inv_lower(p);
}

double lchoose(std::int64_t n, std::int64_t k) {
    if (n < 0 || k < 0 || k > n) return kNaN;
    k = std::min(k, n - k);
    if (k > 100000000) return kNaN;
    Sum s;
    for (std::int64_t i = 1; i <= k; ++i) {
        s.add(std::log(static_cast<double>(n - k + i) / static_cast<double>(i)));
    }
    return s.value();
}

double log_mixture(double pd, double rho, std::int64_t n, std::int64_t d) {
    if (!(pd > 0.0 && pd < 1.0 && rho > 0.0 && rho < 1.0 && n >= 1 && d >= 0 && d <= n)) return kNaN;
    const double c = ncdf_inv(pd);
    const double sr = std::sqrt(rho);
    const double s1 = std::sqrt(1.0 - rho);
    const double dd = static_cast<double>(d);
    const double ss = static_cast<double>(n - d);
    const auto h = [&](double z) {  // log integrand including the N(0,1) kernel (up to 1/sqrt(2 pi))
        const double x = (c - sr * z) / s1;
        double v = -0.5 * z * z;
        if (dd > 0.0) v += dd * log_ncdf(x);
        if (ss > 0.0) v += ss * log_ncdf(-x);
        return v;
    };

    // h is concave (log Phi is), so it is unimodal: bracket its maximum, then golden-section.
    // After expansion h(lo) < h(inner) >= h(outer), so the maximum lies in [lo, hi].
    double lo = -1.0;
    double hi = 1.0;
    if (h(1.0) > h(0.0)) {
        lo = 0.0;
        double step = 1.0;
        while (h(hi + step) > h(hi) && step < 1e6) {
            lo = hi;
            hi += step;
            step *= 2.0;
        }
        hi += step;
    } else if (h(-1.0) > h(0.0)) {
        hi = 0.0;
        double step = 1.0;
        while (h(lo - step) > h(lo) && step < 1e6) {
            hi = lo;
            lo -= step;
            step *= 2.0;
        }
        lo -= step;
    }
    const Max1 peak = golden_max(h, lo, hi, 1e-12 * (1.0 + std::fabs(lo) + std::fabs(hi)));
    const double zm = peak.x;
    const double hm = peak.value;

    // One-sided widths where h has dropped by 1/2, i.e. one sigma on that side of a Gaussian.
    const auto half_drop = [&](double sign) {
        double w_in = 0.0;
        double w_out = 1e-6;
        while (h(zm + sign * w_out) > hm - 0.5 && w_out < 1e6) {
            w_in = w_out;
            w_out *= 2.0;
        }
        for (int it = 0; it < 100; ++it) {
            const double w = 0.5 * (w_in + w_out);
            (h(zm + sign * w) > hm - 0.5 ? w_in : w_out) = w;
        }
        return w_out;
    };
    const double w_left = half_drop(-1.0);
    const double w_right = half_drop(1.0);

    // Panels at the mode +- width * 2^j, out to where the integrand has fallen below e^-60.
    std::vector<double> edges{zm};
    for (double w = w_right; ; w *= 2.0) {
        edges.push_back(zm + w);
        if (h(zm + w) < hm - 60.0 || w > 1e6) break;
    }
    for (double w = w_left; ; w *= 2.0) {
        edges.insert(edges.begin(), zm - w);
        if (h(zm - w) < hm - 60.0 || w > 1e6) break;
    }
    // f = exp(h - hm) is only as accurate as h, whose absolute error is about eps * |hm| (its terms
    // share a sign near the mode). That relative noise is the floor for accepting a panel.
    const auto f = [&](double z) { return std::exp(h(z) - hm); };
    const double tol = 1e-16 * (w_left + w_right) / static_cast<double>(edges.size() - 1);
    const double noise_rel = 64.0 * kEps * (1.0 + std::fabs(hm));
    Sum total;
    for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
        total.add(adaptive(f, edges[i], edges[i + 1], tol, noise_rel, 0));
    }
    return hm + std::log(total.value()) - kLogSqrt2Pi;
}

double period_loglik(double pd, double rho, const Period& y) {
    return lchoose(y.n, y.d) + log_mixture(pd, rho, y.n, y.d);
}

double loglik(double pd, double rho, const std::vector<Period>& panel) {
    Sum s;
    for (const auto& y : panel) s.add(period_loglik(pd, rho, y));
    return s.value();
}

Fit fit(const std::vector<Period>& panel, double pd_lo, double pd_hi, double rho_lo, double rho_hi) {
    const double tol = 1e-8;  // logit units
    const double a_pd = logit(pd_lo);
    const double b_pd = logit(pd_hi);
    const double a_rho = logit(rho_lo);
    const double b_rho = logit(rho_hi);
    const auto ll = [&](double u_pd, double u_rho) { return loglik(inv_logit(u_pd), inv_logit(u_rho), panel); };
    const auto inner = [&](double u_rho) { return golden_max([&](double u) { return ll(u, u_rho); }, a_pd, b_pd, tol); };
    const Max1 outer = golden_max([&](double u_rho) { return inner(u_rho).value; }, a_rho, b_rho, tol);
    const double u_rho = outer.x;
    const double u_pd = inner(u_rho).x;

    Fit r{inv_logit(u_pd), inv_logit(u_rho), ll(u_pd, u_rho), kNaN, kNaN, kNaN, false};
    r.on_boundary = std::fabs(u_pd - a_pd) < 10 * tol || std::fabs(u_pd - b_pd) < 10 * tol ||
                    std::fabs(u_rho - a_rho) < 10 * tol || std::fabs(u_rho - b_rho) < 10 * tol;

    // Observed information in logit coordinates by central differences, then the delta method.
    const double s = 1e-3;
    const double f0 = r.loglik;
    const double h00 = (ll(u_pd + s, u_rho) - 2 * f0 + ll(u_pd - s, u_rho)) / (s * s);
    const double h11 = (ll(u_pd, u_rho + s) - 2 * f0 + ll(u_pd, u_rho - s)) / (s * s);
    const double h01 = (ll(u_pd + s, u_rho + s) - ll(u_pd + s, u_rho - s) - ll(u_pd - s, u_rho + s) +
                        ll(u_pd - s, u_rho - s)) /
                       (4 * s * s);
    const double det = h00 * h11 - h01 * h01;
    if (h00 < 0.0 && det > 0.0) {
        const double j0 = r.pd * (1.0 - r.pd);
        const double j1 = r.rho * (1.0 - r.rho);
        r.se_pd = j0 * std::sqrt(-h11 / det);
        r.se_rho = j1 * std::sqrt(-h00 / det);
        r.corr = (h01 / det) / std::sqrt((h11 / det) * (h00 / det));
    }
    return r;
}

ProfileInterval profile_interval(const std::vector<Period>& panel, const Fit& est, int param, double pd_lo,
                                 double pd_hi, double rho_lo, double rho_hi, double threshold) {
    const double tol = 1e-10;  // logit units, both searches
    // Distinct periods with their counts, so each is integrated once per parameter point.
    std::vector<Period> distinct;
    std::vector<double> count;
    for (const auto& y : panel) {
        std::size_t j = 0;
        while (j < distinct.size() && !(distinct[j].n == y.n && distinct[j].d == y.d)) ++j;
        if (j == distinct.size()) {
            distinct.push_back(y);
            count.push_back(1.0);
        } else {
            count[j] += 1.0;
        }
    }
    const auto ll = [&](double u_pd, double u_rho) {
        Sum s;
        for (std::size_t j = 0; j < distinct.size(); ++j) {
            s.add(count[j] * period_loglik(inv_logit(u_pd), inv_logit(u_rho), distinct[j]));
        }
        return s.value();
    };
    const double box[2][2] = {{logit(pd_lo), logit(pd_hi)}, {logit(rho_lo), logit(rho_hi)}};
    const int other = 1 - param;
    const auto profile = [&](double u) {
        return golden_max([&](double w) { return param == 0 ? ll(u, w) : ll(w, u); }, box[other][0], box[other][1], tol)
            .value;
    };
    const double level = est.loglik - threshold;
    const double centre = logit(param == 0 ? est.pd : est.rho);
    ProfileInterval r{kNaN, kNaN, false, false};
    for (int side = 0; side < 2; ++side) {
        const double bound = box[param][side];
        double inside = centre;
        double outside = bound;
        bool at_bound = profile(bound) >= level;
        if (!at_bound) {
            while (std::fabs(outside - inside) > tol) {
                const double mid = 0.5 * (inside + outside);
                (profile(mid) >= level ? inside : outside) = mid;
            }
        }
        const double v = at_bound ? inv_logit(bound) : inv_logit(0.5 * (inside + outside));
        (side == 0 ? r.lo : r.hi) = v;
        (side == 0 ? r.lo_at_bound : r.hi_at_bound) = at_bound;
    }
    return r;
}

}  // namespace vcalref
