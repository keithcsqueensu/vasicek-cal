// SPDX-License-Identifier: Apache-2.0
//
// Grid-Bayesian estimation (M3; studies/bayes-coverage/PREDICTION.md defines the estimator, D-161).
//
// Everything is in the grid's scaled (logit) coordinates u = (u_PD, u_rho), where grid points are
// equally spaced. At grid point k the log posterior density in u is
//     log pi_u(u_k) + sum_t l_t(theta_k),
// with the prior's density taken in u:
//   Flat      flat on the natural scale over the box: pi_u = PD (1 - PD) rho (1 - rho), the Jacobian;
//   Jeffreys  sqrt det I_u(PD, rho), I_u the Fisher information in u of ONE period with the panel's n
//             (T periods multiply I by T, which cancels). I_u is tabulated once per n on a table grid
//             (jeffreys_table) and interpolated bilinearly in log at other points.
// Each point's mass is its density times its cell's width in u, spread uniformly over the cell
// (half-cells at a grid's ends): the cell-uniform convention. Marginals are the masses summed over the
// other axis, so each marginal CDF is piecewise linear and continuous.
//
// The resolution rule (D-161, as amended in studies/bayes-coverage/PREDICTION.md, 2026-09-28): the grid
// must have at least kPosteriorPointsPerSd points per posterior SD along each axis. The estimator
// evaluates the posterior on the given (parity) grid; while the rule fails, it re-evaluates it on a
// local grid. Per axis its extent is the hull of the posterior mean +- 8 posterior SDs and the current
// edges of the cells holding the marginal's kPosteriorTail and 1 - kPosteriorTail quantiles (so a
// skewed posterior's long tail, such
// as rho's towards the box's floor at small n, stays inside), at least +-2 of the current spacings,
// clipped to the box; the SD and quantiles are measured on the current grid. Its spacing is at most
// s/5, s the smaller of that SD and, on the first refinement, the MLE's Hessian SE in u (when given),
// so the new grid is fine enough even where the coarse grid overstates the SD. The fit is refused and
// flagged after three local grids that still fail the rule, if a local grid would need more than
// kPosteriorMaxPoints points on an axis, or if the current grid's posterior has more than
// kPosteriorMassOutside of its mass outside the next grid.
//
// Intervals, 95%: equal-tailed (the marginal CDF's 2.5% and 97.5% points) and HPD (the shortest
// interval holding 95% of the cell-uniform marginal, a partial cell at each end). Results are on the
// natural scale.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "core/grid.hpp"
#include "engine/refine.hpp"
#include "engine/surface.hpp"

namespace vcal::engine {

enum class Prior { Flat, Jeffreys };

inline constexpr double kPosteriorPointsPerSd = 4.0;
// A local grid is built finer than the rule requires (5 points per SD), so that the SD measured on it,
// which shrinks slightly with the cells (the cell-uniform variance h^2/12), still meets the rule.
inline constexpr double kPosteriorRefinePointsPerSd = 5.0;
inline constexpr double kPosteriorSpanSd = 8.0;
inline constexpr int kPosteriorMaxRefinements = 3;
inline constexpr double kPosteriorMassOutside = 1e-6;
inline constexpr double kPosteriorTail = 5e-8;      // each tail kept inside a local grid, per axis
inline constexpr std::int32_t kPosteriorMaxPoints = 2001;  // per axis of a local grid
inline constexpr double kPosteriorLevel = 0.95;
inline constexpr double kJeffreysStep = 1e-4;   // central-difference step in u for the score
inline constexpr double kJeffreysMassFloor = 1e-14;  // counts d with P(d) below this x max are skipped

enum : std::uint32_t {
    kPosteriorRefined = 1u << 0,       // the rule failed on the given grid; the result is on a local grid
    kPosteriorRefused = 1u << 1,       // the rule still failed after kPosteriorMaxRefinements local grids
    kPosteriorMassOutsideLocal = 1u << 2,  // too much parity-grid mass outside the local grid: refused
    kPosteriorNumeric = 1u << 3,       // no finite posterior mass
};

// Jeffreys prior table: log sqrt det I_u at every point of `grid`.
struct JeffreysTable {
    Grid<2> grid;
    std::vector<double> log_prior;
};

// Bilinear interpolation of log sqrt det I_u at (u0, u1); outside the table, the nearest edge value.
inline double jeffreys_log_prior(const JeffreysTable& jt, double u0, double u1) {
    const Grid<2>& g = jt.grid;
    double f[2];
    std::int32_t i[2];
    const double u[2] = {u0, u1};
    for (int a = 0; a < 2; ++a) {
        const double x = (u[a] - g.axis[a].scaled_lo()) / g.axis[a].step();
        const double c = std::min(std::max(x, 0.0), static_cast<double>(g.axis[a].n - 1));
        i[a] = std::min(static_cast<std::int32_t>(c), g.axis[a].n - 2);
        f[a] = c - static_cast<double>(i[a]);
    }
    const auto at = [&](std::int32_t a, std::int32_t b) {
        const std::int32_t idx[2] = {a, b};
        return jt.log_prior[static_cast<std::size_t>(g.flatten(idx))];
    };
    return (1 - f[0]) * ((1 - f[1]) * at(i[0], i[1]) + f[1] * at(i[0], i[1] + 1)) +
           f[0] * ((1 - f[1]) * at(i[0] + 1, i[1]) + f[1] * at(i[0] + 1, i[1] + 1));
}

// The Jeffreys table for periods of n obligors, on `grid` (logit axes). I_u = sum_d P(d) g_d g_d^T,
// g_d the score in u by central differences of log P(d) (step kJeffreysStep), over the counts whose
// probability is at least kJeffreysMassFloor times the largest. log P(d) is the objective's l_t.
template <class Backend, class Objective, class Integrator>
JeffreysTable jeffreys_table(const Backend& backend, const Objective& objective, const Integrator& integrator,
                             std::int64_t n, const Grid<2>& grid) {
    JeffreysTable jt{grid, std::vector<double>(static_cast<std::size_t>(grid.size()))};
    backend.parallel_for(grid.size(), [&](std::int64_t k) {
        std::int32_t idx[2];
        grid.unflatten(k, idx);
        const double u0 = grid.axis[0].scaled_at(idx[0]), u1 = grid.axis[1].scaled_at(idx[1]);
        const auto l = [&](std::int64_t d, double a, double b) {
            const double v[2] = {grid::from_scaled(grid.axis[0].scale, a), grid::from_scaled(grid.axis[1].scale, b)};
            return objective.log_contrib(typename Objective::Obs{n, d}, Objective::theta(v), integrator);
        };
        std::vector<double> lp(static_cast<std::size_t>(n + 1));
        double top = -HUGE_VAL;
        for (std::int64_t d = 0; d <= n; ++d) {
            lp[static_cast<std::size_t>(d)] = l(d, u0, u1);
            top = std::fmax(top, lp[static_cast<std::size_t>(d)]);
        }
        const double h = kJeffreysStep, floor = top + std::log(kJeffreysMassFloor);
        double i00 = 0.0, i01 = 0.0, i11 = 0.0;
        for (std::int64_t d = 0; d <= n; ++d) {
            if (!(lp[static_cast<std::size_t>(d)] >= floor)) continue;
            const double p = std::exp(lp[static_cast<std::size_t>(d)]);
            const double g0 = (l(d, u0 + h, u1) - l(d, u0 - h, u1)) / (2.0 * h);
            const double g1 = (l(d, u0, u1 + h) - l(d, u0, u1 - h)) / (2.0 * h);
            i00 += p * g0 * g0;
            i01 += p * g0 * g1;
            i11 += p * g1 * g1;
        }
        const double det = i00 * i11 - i01 * i01;
        jt.log_prior[static_cast<std::size_t>(k)] = det > 0.0 ? 0.5 * std::log(det) : -HUGE_VAL;
    });
    return jt;
}

// One axis's marginal on a grid: point masses (normalised) and the cell edges in u.
struct Marginal {
    std::vector<double> mass;   // per axis point, sums to 1
    std::vector<double> lo, hi;  // cell edges in u (half-cells at the ends)
};

struct PosteriorResult {
    double mean[2];           // posterior means in u
    double sd[2];             // posterior SDs in u (cell-uniform: within-cell variance included)
    double et_lo[2], et_hi[2];    // equal-tailed interval, natural scale
    double hpd_lo[2], hpd_hi[2];  // HPD interval, natural scale
    std::uint32_t flags;
    int refinements;          // local grids used
    Grid<2> grid;             // the grid the result is on
};

namespace posterior_detail {

inline double log_prior_at(Prior prior, const JeffreysTable* jt, const Grid<2>& g, double u0, double u1) {
    if (prior == Prior::Jeffreys) return jeffreys_log_prior(*jt, u0, u1);
    const double v0 = grid::from_scaled(g.axis[0].scale, u0), v1 = grid::from_scaled(g.axis[1].scale, u1);
    return std::log(grid::dvalue_dscaled(g.axis[0].scale, u0)) + std::log(grid::dvalue_dscaled(g.axis[1].scale, u1)) +
           0.0 * (v0 + v1);
}

inline double cell_width(const Axis& a, std::int32_t i) {
    return (i == 0 || i == a.n - 1) ? 0.5 * a.step() : a.step();
}

// Point masses on grid g from the log likelihood ll[k]: density x cell area, normalised. Returns false
// if no mass is finite.
inline bool masses(Prior prior, const JeffreysTable* jt, const Grid<2>& g, const std::vector<double>& ll,
                   std::vector<double>& m) {
    const std::int64_t K = g.size();
    std::vector<double> lp(static_cast<std::size_t>(K));
    double top = -HUGE_VAL;
    for (std::int64_t k = 0; k < K; ++k) {
        std::int32_t i[2];
        g.unflatten(k, i);
        const double u0 = g.axis[0].scaled_at(i[0]), u1 = g.axis[1].scaled_at(i[1]);
        const double v = ll[static_cast<std::size_t>(k)] + log_prior_at(prior, jt, g, u0, u1) +
                         std::log(cell_width(g.axis[0], i[0]) * cell_width(g.axis[1], i[1]));
        lp[static_cast<std::size_t>(k)] = v;
        if (v > top) top = v;
    }
    if (!std::isfinite(top)) return false;
    m.assign(static_cast<std::size_t>(K), 0.0);
    double sum = 0.0, comp = 0.0;
    for (std::int64_t k = 0; k < K; ++k) {
        const double x = std::exp(lp[static_cast<std::size_t>(k)] - top);
        m[static_cast<std::size_t>(k)] = std::isfinite(x) ? x : 0.0;
        const double s = sum + m[static_cast<std::size_t>(k)];
        comp += std::fabs(sum) >= m[static_cast<std::size_t>(k)] ? (sum - s) + m[static_cast<std::size_t>(k)] : (m[static_cast<std::size_t>(k)] - s) + sum;
        sum = s;
    }
    const double total = sum + comp;
    for (auto& x : m) x /= total;
    return true;
}

inline Marginal marginal(const Grid<2>& g, const std::vector<double>& m, int a) {
    const Axis& ax = g.axis[a];
    Marginal out{std::vector<double>(static_cast<std::size_t>(ax.n), 0.0), {}, {}};
    for (std::int64_t k = 0; k < g.size(); ++k) {
        std::int32_t i[2];
        g.unflatten(k, i);
        out.mass[static_cast<std::size_t>(i[a])] += m[static_cast<std::size_t>(k)];
    }
    for (std::int32_t i = 0; i < ax.n; ++i) {
        const double u = ax.scaled_at(i), h = ax.step();
        out.lo.push_back(i == 0 ? u : u - 0.5 * h);
        out.hi.push_back(i == ax.n - 1 ? u : u + 0.5 * h);
    }
    return out;
}

// Mean and SD in u under the cell-uniform convention.
inline void moments(const Marginal& mg, double& mean, double& sd) {
    double m1 = 0.0, m2 = 0.0;
    for (std::size_t i = 0; i < mg.mass.size(); ++i) {
        const double a = mg.lo[i], b = mg.hi[i], c = 0.5 * (a + b), w = b - a;
        m1 += mg.mass[i] * c;
        m2 += mg.mass[i] * (c * c + w * w / 12.0);
    }
    mean = m1;
    sd = std::sqrt(std::fmax(m2 - m1 * m1, 0.0));
}

// The u at which the cell-uniform marginal CDF reaches p.
inline double quantile(const Marginal& mg, double p) {
    double cum = 0.0;
    for (std::size_t i = 0; i < mg.mass.size(); ++i) {
        if (cum + mg.mass[i] >= p && mg.mass[i] > 0.0) {
            return mg.lo[i] + (p - cum) / mg.mass[i] * (mg.hi[i] - mg.lo[i]);
        }
        cum += mg.mass[i];
    }
    return mg.hi.back();
}

// The marginal's cumulative mass at u (cell-uniform, piecewise linear), and its inverse.
inline double cdf_at(const Marginal& mg, double u) {
    double cum = 0.0;
    for (std::size_t i = 0; i < mg.mass.size(); ++i) {
        if (u >= mg.hi[i]) {
            cum += mg.mass[i];
        } else {
            if (u > mg.lo[i] && mg.hi[i] > mg.lo[i]) cum += mg.mass[i] * (u - mg.lo[i]) / (mg.hi[i] - mg.lo[i]);
            break;
        }
    }
    return cum;
}

// The edge of the cell in which the cell-uniform marginal CDF reaches p: its lower edge (lower = true)
// or its upper edge. A grid extent snapped to it keeps that cell's point, and so its mass, inside.
inline double quantile_cell_edge(const Marginal& mg, double p, bool lower) {
    double cum = 0.0;
    for (std::size_t i = 0; i < mg.mass.size(); ++i) {
        if (cum + mg.mass[i] >= p && mg.mass[i] > 0.0) return lower ? mg.lo[i] : mg.hi[i];
        cum += mg.mass[i];
    }
    return lower ? mg.lo.front() : mg.hi.back();
}

// HPD: the shortest interval holding `level` of the cell-uniform marginal, which is the highest-
// density interval of that piecewise-constant density, with a partial cell at each end. Its length is
// piecewise linear in its lower end between the points where either end crosses a cell edge, so the
// minimum is at one of those points: every cell edge is tried as the lower end and as the upper end.
//
// Near its minimum that length is almost flat, so candidates up to half a cell apart can differ in
// length by rounding alone, and the last bits of a platform's libm would pick between them. The choice
// is therefore canonical: every candidate within kHpdTieRelative of the shortest counts as shortest,
// and among those the one whose centre is closest to the marginal's median is taken (then the lowest).
inline constexpr double kHpdTieRelative = 1e-9;

inline void hpd(const Marginal& mg, double level, double& lo, double& hi) {
    struct Candidate {
        double a, b;
    };
    std::vector<Candidate> cand;
    for (std::size_t i = 0; i < mg.mass.size(); ++i) {
        for (const double e : {mg.lo[i], mg.hi[i]}) {
            const double c = cdf_at(mg, e);
            if (c + level <= 1.0) cand.push_back({e, quantile(mg, c + level)});  // e as the lower end
            if (c >= level) cand.push_back({quantile(mg, c - level), e});         // e as the upper end
        }
    }
    lo = mg.lo.front();
    hi = mg.hi.back();
    if (cand.empty()) return;
    double shortest = HUGE_VAL;
    for (const auto& c : cand) shortest = std::fmin(shortest, c.b - c.a);
    const double median = quantile(mg, 0.5);
    double best_off = HUGE_VAL;
    for (const auto& c : cand) {
        if (!(c.b - c.a <= shortest * (1.0 + kHpdTieRelative))) continue;
        const double off = std::fabs(0.5 * (c.a + c.b) - median);
        if (off < best_off || (off == best_off && c.a < lo)) {
            best_off = off;
            lo = c.a;
            hi = c.b;
        }
    }
}

}  // namespace posterior_detail

// The grid posterior of a panel. L is the panel's T x K surface on `parity` (from calibrate); se_u,
// when non-null, is the MLE's Hessian SE in u per axis (seeds the local grid's SD). jt is required
// for Prior::Jeffreys.
template <class Backend, class Objective, class Integrator>
PosteriorResult grid_posterior(const Backend& backend, const Objective& objective, const Integrator& integrator,
                               const typename Objective::Obs* obs, std::int64_t periods, const Grid<2>& parity,
                               const std::vector<double>& L, Prior prior, const JeffreysTable* jt,
                               const double* se_u = nullptr) {
    namespace pd = posterior_detail;
    PosteriorResult out{};
    out.grid = parity;
    const std::vector<double> ones(static_cast<std::size_t>(periods), 1.0);
    std::vector<double> ll(static_cast<std::size_t>(parity.size()));
    for (std::int64_t k = 0; k < parity.size(); ++k) {
        ll[static_cast<std::size_t>(k)] = weighted_sum(L.data(), periods, parity.size(), ones.data(), k);
    }
    std::vector<double> m;
    if (!pd::masses(prior, jt, parity, ll, m)) {
        out.flags |= kPosteriorNumeric;
        return out;
    }
    Grid<2> g = parity;
    Marginal mg[2] = {pd::marginal(g, m, 0), pd::marginal(g, m, 1)};
    for (int a = 0; a < 2; ++a) pd::moments(mg[a], out.mean[a], out.sd[a]);
    std::vector<double> Lloc;
    for (;;) {
        const bool ok = out.sd[0] >= kPosteriorPointsPerSd * g.axis[0].step() &&
                        out.sd[1] >= kPosteriorPointsPerSd * g.axis[1].step();
        if (ok) break;
        if (out.refinements == kPosteriorMaxRefinements) {
            out.flags |= kPosteriorRefused;
            return out;
        }
        // The local grid: the hull of mean +- 8 SD and the marginal's tail quantiles, at least +-2 of the
        // current spacings, inside the box.
        Grid<2> next = g;
        bool too_fine = false;
        double span_lo[2], span_hi[2];
        for (int a = 0; a < 2; ++a) {
            double fine = out.sd[a];
            if (out.refinements == 0 && se_u != nullptr && std::isfinite(se_u[a]) && se_u[a] > 0.0) {
                fine = std::fmin(fine, se_u[a]);
            }
            const double half = std::fmax(kPosteriorSpanSd * out.sd[a], 2.0 * g.axis[a].step());
            const double lo = std::fmax(std::fmin(out.mean[a] - half, pd::quantile_cell_edge(mg[a], kPosteriorTail, true)),
                                        parity.axis[a].scaled_lo());
            const double hi = std::fmin(std::fmax(out.mean[a] + half, pd::quantile_cell_edge(mg[a], 1.0 - kPosteriorTail, false)),
                                        parity.axis[a].scaled_hi());
            const double want = std::ceil((hi - lo) / (fine / kPosteriorRefinePointsPerSd)) + 1.0;
            if (!(want <= static_cast<double>(kPosteriorMaxPoints))) too_fine = true;
            span_lo[a] = lo;
            span_hi[a] = hi;
            next.axis[a].lo = grid::from_scaled(parity.axis[a].scale, lo);
            next.axis[a].hi = grid::from_scaled(parity.axis[a].scale, hi);
            next.axis[a].n = too_fine ? 3 : std::max<std::int32_t>(static_cast<std::int32_t>(want), 3);
        }
        if (too_fine) {
            out.flags |= kPosteriorRefused;
            return out;
        }
        // The current grid's mass outside the next grid.
        double outside = 0.0;
        for (std::int64_t k = 0; k < g.size(); ++k) {
            std::int32_t i[2];
            g.unflatten(k, i);
            for (int a = 0; a < 2; ++a) {
                const double u = g.axis[a].scaled_at(i[a]);
                if (u < span_lo[a] || u > span_hi[a]) {
                    outside += m[static_cast<std::size_t>(k)];
                    break;
                }
            }
        }
        if (outside > kPosteriorMassOutside) {
            out.flags |= kPosteriorMassOutsideLocal | kPosteriorRefused;
            return out;
        }
        g = next;
        Lloc.assign(static_cast<std::size_t>(periods * g.size()), 0.0);
        evaluate_surface(backend, objective, integrator, obs, periods, g, Lloc.data());
        ll.assign(static_cast<std::size_t>(g.size()), 0.0);
        for (std::int64_t k = 0; k < g.size(); ++k) {
            ll[static_cast<std::size_t>(k)] = weighted_sum(Lloc.data(), periods, g.size(), ones.data(), k);
        }
        if (!pd::masses(prior, jt, g, ll, m)) {
            out.flags |= kPosteriorNumeric;
            return out;
        }
        mg[0] = pd::marginal(g, m, 0);
        mg[1] = pd::marginal(g, m, 1);
        for (int a = 0; a < 2; ++a) pd::moments(mg[a], out.mean[a], out.sd[a]);
        ++out.refinements;
        out.flags |= kPosteriorRefined;
    }
    out.grid = g;
    const double tail = 0.5 * (1.0 - kPosteriorLevel);
    for (int a = 0; a < 2; ++a) {
        const AxisScale s = g.axis[a].scale;
        out.et_lo[a] = grid::from_scaled(s, pd::quantile(mg[a], tail));
        out.et_hi[a] = grid::from_scaled(s, pd::quantile(mg[a], 1.0 - tail));
        double lo = 0.0, hi = 0.0;
        pd::hpd(mg[a], kPosteriorLevel, lo, hi);
        out.hpd_lo[a] = grid::from_scaled(s, lo);
        out.hpd_hi[a] = grid::from_scaled(s, hi);
    }
    return out;
}

}  // namespace vcal::engine
