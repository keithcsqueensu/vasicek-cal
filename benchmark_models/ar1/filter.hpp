// SPDX-License-Identifier: Apache-2.0
//
// The binomial-mixture log-likelihood under an AR(1) factor (M6; D-195; docs/methodology/ar1_binomial_mle.md):
//     Z_1 ~ N(0, 1),  Z_t | Z_{t-1} ~ N(phi Z_{t-1}, 1 - phi^2),  p_t = p(Z_t; PD, rho),
//     L = integral of prod_t C(n_t, d_t) p_t^d_t (1 - p_t)^(n_t - d_t)  dP(Z_1, ..., Z_T),
// by a forward filter over each period's own quadrature nodes.
//
// Period t's integral is  f_t(z) pi_t(z) dz,  with f_t the binomial term and pi_t the predictive density of
// Z_t given d_1..d_{t-1} (pi_1 = N(0, 1)). Its nodes are placed for a Gaussian proxy q_t = N(mu_t, sigma_t^2):
// with z = mu_t + sigma_t y,
//     integral f_t(z) q_t(z) dz = integral f_t(mu_t + sigma_t y) phi(y) dy,
// and f_t(mu + sigma y) is again a Vasicek binomial term, with c' = c - sqrt(rho) mu and sqrt(rho)' =
// sqrt(rho) sigma. So the per-period rule the binomial likelihood uses (parity: adaptive Gauss-Hermite, or
// composite Gauss-Legendre for d in {0, n}, D-118) integrates it unchanged, hint included. The exact
// predictive enters through the ratio pi_t / q_t at each node, folded into the rule's accumulator:
//     M(z) = LSE_j(lhat_j + log K(z_j -> z)) - LSE_j(lhat_j) - log q_t(z),
// over the previous period's nodes z_j with normalised log-weights lhat_j. Any proxy gives the same integral
// up to quadrature error; the proxy only decides where the nodes go. Everything is in log space, and lhat is
// normalised after each period so the difference never cancels large magnitudes.
//
// Two placements (D-195):
//   - predictive: q_t moment-matched to pi_t from the filtered weights, mu_t = phi m_{t-1}, sigma_t^2 =
//     phi^2 s_{t-1}^2 + 1 - phi^2. One forward pass. It still misses mass when a later period selects an
//     earlier period's tail (a zero-default year before a crisis year at n = 10,000, phi = 0.95: 1e-4 nats),
//     since a forward filter's nodes cannot see the future.
//   - neighbours (the estimate): q_t the law of Z_t given its neighbours, at their smoothed moments m^, s^2
//     from a backward pass over the predictive placement's nodes:
//         mu_t = phi (m^_{t-1} + m^_{t+1}) / (1 + phi^2),
//         sigma_t^2 = (1 - phi^2) / (1 + phi^2) + phi^2 (s^2_{t-1} + s^2_{t+1}) / (1 + phi^2)^2,
//     one neighbour at either end (mu = phi m^, sigma^2 = 1 - phi^2 + phi^2 s^2). Forward, backward, forward.
// The difference between the two is the per-run placement check: the doubled node count cannot see a
// misplaced centre, since both rules share it. tests/unit/ar1_filter_test.cpp holds the hard cases against an
// independent dense-grid filter.
//
// phi = 0 gives the per-period likelihood bit for bit under both placements, with no special case: mu = 0
// and sigma^2 = 1 exactly (every phi term is 0 * finite), so the transformed model's fields equal the
// original's and z = 0 + 1 y = y; the log-ratio below is then exactly 0 (0.5 y^2 - y^2 * 0.5 and
// 0.5 log 1 - 0.5 log 1); M = LSE(lhat) - LSE(lhat) = 0 exactly; each period's term is
// BinomialMixture::log_contrib's; and the periods are summed by engine::weighted_sum in period order, as
// calibrate does.
#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

#include "core/model/binomial_mixture.hpp"
#include "core/model/vasicek.hpp"
#include "core/quadrature/integrator.hpp"
#include "core/special/lbinom.hpp"
#include "engine/surface.hpp"

namespace vcal::benchmark_models::ar1 {

inline constexpr double kPhiBound = 0.95;  // |phi| <= 0.95 (D-195): 1 - phi^2 >= 0.0975

struct Obs {
    std::int64_t n;
    std::int64_t d;
};

struct Theta {
    double pd;
    double rho;
    double phi;
};

// Neighbours is the estimate; Predictive is its first pass and the placement check's comparator. Prior (nodes
// placed for N(0, 1) in every period, the first design) is kept only to reproduce D-195's finding in the tests:
// with it the doubled-node check agrees to 1e-11 while both rules are 0.9 nats off.
enum class Placement { Prior, Predictive, Neighbours };

// One period's nodes from the last forward pass.
struct PeriodNodes {
    std::vector<double> z;     // abscissae
    std::vector<double> y;     // the same in the proxy's coordinates, z = mu + sigma y
    std::vector<double> v;     // the rule's own log-terms (weight, binomial term, proxy density)
    std::vector<double> lhat;  // normalised filtered log-weights
    double mu = 0.0, sigma = 1.0;
    double lhat_lse = 0.0;     // LSE of lhat in node order
};

struct FilterWork {
    std::vector<PeriodNodes> nodes;
    std::vector<double> contrib;           // per-period terms l_t = log p(d_t | d_1..d_{t-1})
    std::vector<double> mean, var;         // smoothed moments of Z_t (from smooth())
    std::vector<double> mu, sigma;         // proxies for a Neighbours pass
};

namespace detail {

// log K(u -> z) - log q(z) at a node y of the proxy (z = mu + sigma y), without the common log sqrt(2 pi):
//     0.5 y^2 + log sigma - (z - phi u)^2 / (2 (1 - phi^2)) - 0.5 log(1 - phi^2).
// At phi = 0 (mu = 0, sigma = 1, z = y) every term is exact and the result is 0.
struct ProxyTransition {
    double phi;
    double inv_two_var;   // 1 / (2 (1 - phi^2))
    double half_log_var;  // 0.5 log(1 - phi^2)
    double half_log_s2;   // 0.5 log(sigma^2) = log sigma

    double log_ratio(double u, double z, double y) const {
        const double a = z - phi * u;
        return (0.5 * (y * y) + half_log_s2) - ((a * a) * inv_two_var + half_log_var);
    }
    // The first period's: log phi(z) - log q(z) = 0.5 y^2 + log sigma - 0.5 z^2, exactly 0 when q = N(0, 1).
    double log_ratio_first(double z, double y) const { return (0.5 * (y * y) + half_log_s2) - 0.5 * (z * z); }
};

// The rules' accumulator, with the predictive-to-proxy ratio folded into each node's term. In the first period
// the predictive is the prior N(0, 1); with the prior as proxy (the Predictive placement, and every placement
// at phi = 0) the ratio is exactly 0 and the accumulator is OnlineLogSumExp exactly.
class FilterAcc {
public:
    FilterAcc(const ProxyTransition& tr, const PeriodNodes* prev, PeriodNodes& out, std::vector<double>& w)
        : tr_(tr), prev_(prev), out_(out), w_(w) {
        out_.z.clear();
        out_.y.clear();
        out_.v.clear();
        w_.clear();
    }

    void add(double v, double y) {
        const double z = out_.mu + out_.sigma * y;
        double w;
        if (prev_ != nullptr) {
            quadrature::detail::OnlineLogSumExp m;
            const std::size_t np = prev_->z.size();
            for (std::size_t j = 0; j < np; ++j) m.add(prev_->lhat[j] + tr_.log_ratio(prev_->z[j], z, y));
            w = v + (m.result() - prev_->lhat_lse);
        } else {
            w = v + tr_.log_ratio_first(z, y);
        }
        acc_.add(w);
        out_.z.push_back(z);
        out_.y.push_back(y);
        out_.v.push_back(v);
        w_.push_back(w);
    }

    double result() const { return acc_.result(); }

private:
    const ProxyTransition& tr_;
    const PeriodNodes* prev_;
    PeriodNodes& out_;
    std::vector<double>& w_;
    quadrature::detail::OnlineLogSumExp acc_;
};

// Mean and variance of z under the weights exp(lw + extra) (extra optional), normalised. Two passes, so the
// variance is a sum of non-negative terms: no floor or clamp, which would put a kink into the likelihood
// through the node placement (D-195's smoothness condition).
inline void weighted_moments(const std::vector<double>& z, const std::vector<double>& lw,
                             const std::vector<double>* extra, double& mean, double& var) {
    const std::size_t n = z.size();
    double mx = -HUGE_VAL;
    for (std::size_t k = 0; k < n; ++k) mx = std::fmax(mx, lw[k] + (extra ? (*extra)[k] : 0.0));
    double sw = 0.0, szw = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        const double wk = std::exp(lw[k] + (extra ? (*extra)[k] : 0.0) - mx);
        sw += wk;
        szw += wk * z[k];
    }
    mean = szw / sw;
    double sdd = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        const double wk = std::exp(lw[k] + (extra ? (*extra)[k] : 0.0) - mx);
        const double e = z[k] - mean;
        sdd += wk * (e * e);
    }
    var = sdd / sw;
}

// One forward pass. Predictive: the proxies come from the filtered moments as the pass goes; Neighbours: from
// work.mu and work.sigma; Prior: N(0, 1) throughout.
template <class Integrator>
double forward(const Obs* obs, std::int64_t periods, const Theta& th, const Integrator& integrator,
               Placement placement, FilterWork& work) {
    const auto base = model::make_vasicek1f(th.pd, th.rho);
    const double phi = th.phi;
    const double var = 1.0 - phi * phi;
    ProxyTransition tr{phi, 1.0 / (2.0 * var), 0.5 * std::log(var), 0.0};
    work.nodes.resize(static_cast<std::size_t>(periods));
    work.contrib.assign(static_cast<std::size_t>(periods), 0.0);
    std::vector<double> w;
    double m = 0.0, s2 = 1.0;  // filtered mean and variance of the previous period's Z
    for (std::int64_t t = 0; t < periods; ++t) {
        const auto ts = static_cast<std::size_t>(t);
        PeriodNodes& cur = work.nodes[ts];
        double sigma2 = 1.0;
        cur.mu = 0.0;
        if (placement == Placement::Neighbours) {
            cur.mu = work.mu[ts];
            sigma2 = work.sigma[ts] * work.sigma[ts];
        } else if (placement == Placement::Predictive && t > 0) {
            cur.mu = phi * m;
            sigma2 = phi * phi * s2 + var;
        }
        cur.sigma = std::sqrt(sigma2);
        tr.half_log_s2 = 0.5 * std::log(sigma2);
        const model::Vasicek1F shifted{base.c - base.sqrt_rho * cur.mu, base.sqrt_rho * cur.sigma,
                                       base.sqrt_one_minus_rho};
        const auto f = model::make_binomial_log_integrand(shifted, obs[t].n, obs[t].d);
        FilterAcc acc(tr, t == 0 ? nullptr : &work.nodes[ts - 1], cur, w);
        const double l = special::lbinom(obs[t].n, obs[t].d) +
                         integrator.log_integrate_with(f, model::binomial_mixture_hint(f), acc);
        work.contrib[ts] = l;
        if (!std::isfinite(l)) {
            for (std::int64_t u = t + 1; u < periods; ++u) work.contrib[static_cast<std::size_t>(u)] = l;
            return l;
        }
        // Normalise this period's weights for the next step. Their own log-sum-exp, folded in node order as the
        // next step's numerator folds them, is the reference that makes M exact at phi = 0.
        const double lse_w = acc.result();
        cur.lhat.resize(w.size());
        quadrature::detail::OnlineLogSumExp s;
        for (std::size_t k = 0; k < w.size(); ++k) {
            cur.lhat[k] = w[k] - lse_w;
            s.add(cur.lhat[k]);
        }
        cur.lhat_lse = s.result();
        weighted_moments(cur.z, cur.lhat, nullptr, m, s2);
    }
    const std::vector<double> ones(static_cast<std::size_t>(periods), 1.0);
    return engine::weighted_sum(work.contrib.data(), periods, 1, ones.data(), 0);
}

// Smoothed means and variances of each Z_t over the last forward pass's nodes, by a backward pass:
//     beta_T = 0,  beta_t(j) = LSE_k(log K(z_tj -> z_{t+1,k}) - log q_{t+1}(z_{t+1,k}) + v_{t+1,k} + beta_{t+1}(k)),
// smoothed weights exp(lhat_t + beta_t), normalised.
inline void smooth(const Theta& th, FilterWork& work) {
    const std::size_t T = work.nodes.size();
    const double phi = th.phi;
    const double var = 1.0 - phi * phi;
    ProxyTransition tr{phi, 1.0 / (2.0 * var), 0.5 * std::log(var), 0.0};
    work.mean.assign(T, 0.0);
    work.var.assign(T, 1.0);
    std::vector<double> beta_next, beta;
    for (std::size_t tt = T; tt-- > 0;) {
        const PeriodNodes& cur = work.nodes[tt];
        beta.assign(cur.z.size(), 0.0);
        if (tt + 1 < T) {
            const PeriodNodes& nx = work.nodes[tt + 1];
            tr.half_log_s2 = std::log(nx.sigma);
            for (std::size_t j = 0; j < cur.z.size(); ++j) {
                quadrature::detail::OnlineLogSumExp b;
                for (std::size_t k = 0; k < nx.z.size(); ++k) {
                    b.add(tr.log_ratio(cur.z[j], nx.z[k], nx.y[k]) + nx.v[k] + beta_next[k]);
                }
                beta[j] = b.result();
            }
            // Shift for range only; the smoothed weights are normalised below.
            quadrature::detail::OnlineLogSumExp s;
            for (const double x : beta) s.add(x);
            const double shift = s.result();
            for (double& x : beta) x -= shift;
        }
        weighted_moments(cur.z, cur.lhat, &beta, work.mean[tt], work.var[tt]);
        beta_next.swap(beta);
    }
}

// The proxies of the Neighbours placement from the smoothed moments.
inline void neighbour_proxies(const Theta& th, FilterWork& work) {
    const std::size_t T = work.mean.size();
    const double phi = th.phi;
    const double var = 1.0 - phi * phi;
    const double den = 1.0 + phi * phi;
    work.mu.assign(T, 0.0);
    work.sigma.assign(T, 1.0);
    for (std::size_t t = 0; t < T; ++t) {
        const bool has_prev = t > 0;
        const bool has_next = t + 1 < T;
        double mu, s2;
        if (has_prev && has_next) {
            mu = phi * (work.mean[t - 1] + work.mean[t + 1]) / den;
            s2 = var / den + phi * phi * (work.var[t - 1] + work.var[t + 1]) / (den * den);
        } else if (has_prev || has_next) {
            const std::size_t u = has_prev ? t - 1 : t + 1;
            mu = phi * work.mean[u];
            s2 = var + phi * phi * work.var[u];
        } else {
            mu = 0.0 * phi;
            s2 = 1.0;
        }
        work.mu[t] = mu;
        work.sigma[t] = std::sqrt(s2);
    }
}

}  // namespace detail

// The score (d/dPD, d/drho, d/dphi) of the log-likelihood, by Fisher's identity: the expectation, under the
// smoothed law of Z given the data, of the gradient of the complete-data log-density. Over the last forward
// pass's nodes:
//     d/d(PD, rho):  sum_t E[g_t,theta(Z_t)], the binomial terms' parameter derivatives (core's
//                    BinomialLogIntegrandDerivs) under each period's smoothed node weights;
//     d/dphi:        sum_{t >= 2} E[d log K(Z_{t-1} -> Z_t) / dphi] under the smoothed pairwise weights,
//                    d log K / dphi = a u / v - phi a^2 / v^2 + phi / v,  a = z - phi u,  v = 1 - phi^2.
// Call right after log_likelihood at the same theta (it reads work.nodes). It differs from the derivative of the
// quadrature approximation by the nodes' own movement with theta, a term of the size of the quadrature error;
// tests/unit/ar1_filter_test.cpp checks it against Richardson-extrapolated differences.
inline void score(const Obs* obs, std::int64_t periods, const Theta& th, FilterWork& work, double (&grad)[3]) {
    const std::size_t T = static_cast<std::size_t>(periods);
    const double phi = th.phi;
    const double var = 1.0 - phi * phi;
    detail::ProxyTransition tr{phi, 1.0 / (2.0 * var), 0.5 * std::log(var), 0.0};
    grad[0] = grad[1] = grad[2] = 0.0;
    std::vector<double> beta_next, beta, pair;
    for (std::size_t tt = T; tt-- > 0;) {
        const PeriodNodes& cur = work.nodes[tt];
        const std::size_t nc = cur.z.size();
        beta.assign(nc, 0.0);
        if (tt + 1 < T) {
            const PeriodNodes& nx = work.nodes[tt + 1];
            const std::size_t nn = nx.z.size();
            tr.half_log_s2 = std::log(nx.sigma);
            // Pairwise smoothed weights of (Z_t, Z_{t+1}): lhat_t(j) + log R(j -> k) + v_{t+1}(k) + beta_{t+1}(k).
            pair.resize(nc * nn);
            double mx = -HUGE_VAL;
            for (std::size_t j = 0; j < nc; ++j) {
                quadrature::detail::OnlineLogSumExp b;
                for (std::size_t k = 0; k < nn; ++k) {
                    const double e = tr.log_ratio(cur.z[j], nx.z[k], nx.y[k]) + nx.v[k] + beta_next[k];
                    b.add(e);
                    pair[j * nn + k] = cur.lhat[j] + e;
                    mx = std::fmax(mx, pair[j * nn + k]);
                }
                beta[j] = b.result();
            }
            double sw = 0.0, sd = 0.0;
            for (std::size_t j = 0; j < nc; ++j) {
                const double u = cur.z[j];
                for (std::size_t k = 0; k < nn; ++k) {
                    const double wjk = std::exp(pair[j * nn + k] - mx);
                    if (wjk == 0.0) continue;
                    const double a = nx.z[k] - phi * u;
                    sw += wjk;
                    sd += wjk * (a * u / var - phi * (a * a) / (var * var) + phi / var);
                }
            }
            grad[2] += sd / sw;
            quadrature::detail::OnlineLogSumExp s;
            for (const double x : beta) s.add(x);
            const double shift = s.result();
            for (double& x : beta) x -= shift;
        }
        // This period's smoothed node weights, and the binomial term's parameter derivatives at its nodes.
        const auto df = model::make_binomial_log_integrand_derivs(th.pd, th.rho, obs[tt].n, obs[tt].d);
        double mx = -HUGE_VAL;
        for (std::size_t k = 0; k < nc; ++k) mx = std::fmax(mx, cur.lhat[k] + beta[k]);
        double sw = 0.0, s0 = 0.0, s1 = 0.0;
        for (std::size_t k = 0; k < nc; ++k) {
            const double wk = std::exp(cur.lhat[k] + beta[k] - mx);
            if (wk == 0.0) continue;
            double f[5];
            df(cur.z[k], f);
            sw += wk;
            s0 += wk * f[0];
            s1 += wk * f[1];
        }
        grad[0] += s0 / sw;
        grad[1] += s1 / sw;
        beta_next.swap(beta);
    }
}

struct LogLikelihood {
    double value;      // the estimate (Neighbours placement)
    double placement;  // |Neighbours - Predictive|: the per-run placement check
};

// The log-likelihood at theta (|phi| < 1; the search keeps |phi| <= kPhiBound): forward (Predictive),
// backward, forward (Neighbours). work.contrib holds the Neighbours pass's per-period terms. integrator is the
// per-period rule (quadrature::parity_rule(); parity_rule(true) for the doubled check). -inf when the data are
// impossible at theta, NaN for a numeric failure.
template <class Integrator>
LogLikelihood log_likelihood(const Obs* obs, std::int64_t periods, const Theta& th, const Integrator& integrator,
                             FilterWork& work) {
    const double first = detail::forward(obs, periods, th, integrator, Placement::Predictive, work);
    if (!std::isfinite(first)) return {first, 0.0};
    detail::smooth(th, work);
    detail::neighbour_proxies(th, work);
    const double second = detail::forward(obs, periods, th, integrator, Placement::Neighbours, work);
    return {second, std::fabs(second - first)};
}

}  // namespace vcal::benchmark_models::ar1
