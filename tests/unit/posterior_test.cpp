// SPDX-License-Identifier: Apache-2.0
// engine/posterior.hpp (M3): the grid posterior, its resolution rule and intervals on a posterior
// known exactly, the rule's refusal, thread-count determinism, and the Jeffreys table.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "engine/calibrate.hpp"
#include "engine/posterior.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/recovery/posterior_panels.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace e = vcal::engine;
namespace tol = vcal::tol;
using Backend = vcal::backends::CpuBackend;

vcal::Grid<2> box_grid() {
    return {{{1e-4, 0.2, 61, vcal::AxisScale::Logit}, {1e-3, 0.5, 41, vcal::AxisScale::Logit}}};
}

// A test objective whose single period's likelihood, times the flat prior's density in u, is a
// normal density in u = (logit PD, logit rho): the posterior is then exactly normal with means mu and
// SDs sigma (and, far from the box's bounds, untruncated). A second mode can be added to make it
// bimodal.
struct NormalInU {
    struct Obs {
        int unused;
        friend bool operator==(const Obs&, const Obs&) { return true; }
    };
    struct Theta {
        double pd, rho;
    };
    static constexpr int n_params = 2;
    double mu[2], sigma[2];
    double second_mode_offset = 0.0;  // in u_PD; 0 = unimodal

    static Theta theta(const double (&v)[2]) { return {v[0], v[1]}; }
    static const char* panel_error(const Obs*, std::int64_t) { return nullptr; }
    static double rounding_scale(const Obs&, double l) { return std::fabs(l) + 1.0; }

    template <class I>
    double log_contrib(const Obs&, const Theta& th, const I&) const {
        const double u0 = std::log(th.pd) - std::log1p(-th.pd), u1 = std::log(th.rho) - std::log1p(-th.rho);
        const auto q = [&](double c0) {
            const double z0 = (u0 - c0) / sigma[0], z1 = (u1 - mu[1]) / sigma[1];
            return -0.5 * (z0 * z0 + z1 * z1);
        };
        double lq = q(mu[0]);
        if (second_mode_offset != 0.0) {
            const double b = q(mu[0] + second_mode_offset);
            lq = std::fmax(lq, b) + std::log1p(std::exp(-std::fabs(lq - b)));
        }
        // Divide out the flat prior's density in u (the Jacobian), so the posterior in u is lq alone.
        return lq - std::log(th.pd * (1.0 - th.pd)) - std::log(th.rho * (1.0 - th.rho));
    }
};

e::PosteriorResult posterior_of(const NormalInU& obj, const Backend& backend = Backend{1},
                                const double* se_u = nullptr, bool rule = true) {
    const auto integrator = vcal::quadrature::parity_rule();
    const auto g = box_grid();
    const NormalInU::Obs obs[1] = {{0}};
    std::vector<double> L(static_cast<std::size_t>(g.size()));
    e::evaluate_surface(backend, obj, integrator, obs, 1, g, L.data());
    return e::grid_posterior(backend, obj, integrator, obs, 1, g, L, e::Prior::Flat, nullptr, se_u, rule);
}

double logit(double v) { return std::log(v) - std::log1p(-v); }

std::string sci(double v) {
    char b[32];
    std::snprintf(b, sizeof b, "%.3g", v);
    return b;
}

}  // namespace

// On an exactly normal posterior the intervals are mu +- 1.959964 sigma. Equal-tailed ends, the
// mean and the SD are checked to TOL_POSTERIOR_ET_END_SD; HPD ends to TOL_POSTERIOR_HPD_END_SD, since
// the HPD of the cell-uniform (step) density can sit up to half a cell from the smooth one. With SDs from 0.02 to 0.8
// in u (0.15 to 6 parity spacings), the mean at least 4 SDs inside the box so the box does not
// truncate it, the rule refines where it must, and every interval end is within
// TOL_POSTERIOR_NORMAL_END_SD standard deviations of the exact one.
VCAL_TEST(posterior_intervals_are_exact_on_a_normal_posterior) {
    const double z = 1.959963984540054;
    double worst = 0.0, worst_hpd = 0.0;
    int refined = 0, cases = 0;
    for (const double s0 : {0.02, 0.06, 0.15, 0.4, 0.8}) {
        for (const double s1 : {0.03, 0.2, 0.45}) {
            NormalInU obj{{logit(0.01), logit(0.12)}, {s0, s1}};
            const auto r = posterior_of(obj);
            VCAL_CHECK((r.flags & (e::kPosteriorRefused | e::kPosteriorNumeric)) == 0u);
            refined += (r.flags & e::kPosteriorRefined) ? 1 : 0;
            ++cases;
            const double sg[2] = {s0, s1};
            for (int a = 0; a < 2; ++a) {
                const double lo = obj.mu[a] - z * sg[a], hi = obj.mu[a] + z * sg[a];
                for (const double got : {logit(r.et_lo[a]) - lo, logit(r.et_hi[a]) - hi}) {
                    worst = std::fmax(worst, std::fabs(got) / sg[a]);
                    VCAL_CHECK(std::fabs(got) / sg[a] <= tol::TOL_POSTERIOR_ET_END_SD);
                }
                for (const double got : {logit(r.hpd_lo[a]) - lo, logit(r.hpd_hi[a]) - hi}) {
                    worst_hpd = std::fmax(worst_hpd, std::fabs(got) / sg[a]);
                    VCAL_CHECK(std::fabs(got) / sg[a] <= tol::TOL_POSTERIOR_HPD_END_SD);
                }
                VCAL_CHECK(std::fabs(r.mean[a] - obj.mu[a]) / sg[a] <= tol::TOL_POSTERIOR_ET_END_SD);
                VCAL_CHECK(std::fabs(r.sd[a] / sg[a] - 1.0) <= tol::TOL_POSTERIOR_ET_END_SD);
            }
        }
    }
    // One case wide enough on both axes that the parity grid already meets the rule.
    {
        NormalInU obj{{logit(0.005), logit(0.03)}, {0.8, 0.8}};
        const auto r = posterior_of(obj);
        VCAL_CHECK(r.flags == 0u);
        ++cases;
        for (int a = 0; a < 2; ++a) {
            const double lo = obj.mu[a] - z * 0.8, hi = obj.mu[a] + z * 0.8;
            for (const double got : {logit(r.et_lo[a]) - lo, logit(r.et_hi[a]) - hi}) {
                worst = std::fmax(worst, std::fabs(got) / 0.8);
                VCAL_CHECK(std::fabs(got) / 0.8 <= tol::TOL_POSTERIOR_ET_END_SD);
            }
            for (const double got : {logit(r.hpd_lo[a]) - lo, logit(r.hpd_hi[a]) - hi}) {
                worst_hpd = std::fmax(worst_hpd, std::fabs(got) / 0.8);
                VCAL_CHECK(std::fabs(got) / 0.8 <= tol::TOL_POSTERIOR_HPD_END_SD);
            }
        }
    }
    vcal::test::note(std::to_string(refined) + " of " + std::to_string(cases) +
                     " refined; worst end error: equal-tailed " + sci(worst) + " SD, HPD " + sci(worst_hpd) + " SD");
    VCAL_CHECK(refined > 0 && refined < cases);  // both branches of the rule were exercised
}

// A bimodal posterior with a narrow Hessian SE at one mode: the local grid's extent follows the
// posterior SD, so it keeps both modes and the fit succeeds (a grid of +-8 SEs would have dropped one;
// that was the registered rule's flaw, amended 2026-09-28). The refinement's spacing follows the SE.
VCAL_TEST(posterior_local_grid_keeps_a_skewed_or_bimodal_posterior) {
    NormalInU obj{{logit(0.003), logit(0.12)}, {0.03, 0.03}};
    obj.second_mode_offset = 3.0;
    const double se_u[2] = {0.03, 0.03};
    const auto r = posterior_of(obj, Backend{1}, se_u);
    VCAL_CHECK((r.flags & e::kPosteriorRefused) == 0u);
    VCAL_CHECK(r.flags & e::kPosteriorRefined);
    // Each mode holds half the mass: the equal-tailed PD interval spans both.
    VCAL_CHECK(logit(r.et_lo[0]) < obj.mu[0] && logit(r.et_hi[0]) > obj.mu[0] + 3.0);
}

// The rule switched off (S-9's and S-15's diagnostic arm): a posterior far narrower than a parity
// cell stays on the parity grid, neither refined nor refused, and its mass sits in about one cell, so
// the equal-tailed interval is about a cell wide. marginal_cdf is the inverse of the interval's
// quantiles: 0.025 and 0.975 at its ends, on either grid.
VCAL_TEST(posterior_without_the_rule_stays_on_the_parity_grid) {
    NormalInU obj{{logit(0.01), logit(0.12)}, {0.01, 0.01}};
    const auto off = posterior_of(obj, Backend{1}, nullptr, false);
    VCAL_CHECK_EQ(off.flags, 0u);
    VCAL_CHECK_EQ(off.refinements, 0);
    VCAL_CHECK(off.grid.axis[0].n == 61 && off.grid.axis[1].n == 41);
    const auto on = posterior_of(obj);
    VCAL_CHECK(on.flags & e::kPosteriorRefined);
    const double step = box_grid().axis[0].step();
    VCAL_CHECK(logit(off.et_hi[0]) - logit(off.et_lo[0]) > 0.5 * step);
    // Refined: 2 x 1.959964 sigma, to the normal case's own tolerance.
    VCAL_CHECK(std::fabs(logit(on.et_hi[0]) - logit(on.et_lo[0]) - 2.0 * 1.959963984540054 * obj.sigma[0]) <
               2.0 * tol::TOL_POSTERIOR_ET_END_SD * obj.sigma[0]);
    std::printf("rule off: PD interval %.3f parity spacings; rule on: %.3f\n",
                (logit(off.et_hi[0]) - logit(off.et_lo[0])) / step, (logit(on.et_hi[0]) - logit(on.et_lo[0])) / step);
    double worst = 0.0;
    for (const auto* r : {&off, &on}) {
        for (int a = 0; a < 2; ++a) {
            const double lo = e::marginal_cdf(r->marginal[a], logit(r->et_lo[a]));
            const double hi = e::marginal_cdf(r->marginal[a], logit(r->et_hi[a]));
            const double top = e::marginal_cdf(r->marginal[a], 1e9);
            worst = std::fmax(worst, std::fmax(std::fabs(lo - 0.025), std::fmax(std::fabs(hi - 0.975), std::fabs(top - 1.0))));
            VCAL_CHECK(e::marginal_cdf(r->marginal[a], -1e9) == 0.0);
        }
    }
    std::printf("marginal_cdf round trip: worst %.3g\n", worst);
    VCAL_CHECK(worst <= tol::TOL_POSTERIOR_CDF_ROUNDTRIP_ABS);
}

// A posterior narrower than three local grids can resolve (sigma 1e-7 in u): refused, not reported.
VCAL_TEST(posterior_refuses_after_three_refinements) {
    NormalInU obj{{logit(0.01), logit(0.12)}, {1e-7, 0.2}};
    const auto r = posterior_of(obj);
    VCAL_CHECK(r.flags & e::kPosteriorRefused);
    VCAL_CHECK_EQ(r.refinements, e::kPosteriorMaxRefinements);
}

// A posterior piled against the box's lower rho bound (a normal in u centred near it, truncated by the
// box): the local grid reaches the bound, and the fit is not refused.
VCAL_TEST(posterior_near_a_bound_is_not_refused) {
    NormalInU obj{{logit(0.01), logit(0.002)}, {0.05, 1.5}};
    const auto r = posterior_of(obj);
    VCAL_CHECK((r.flags & e::kPosteriorRefused) == 0u);
    VCAL_CHECK(std::fabs(r.grid.axis[1].lo / 1e-3 - 1.0) < 1e-12);  // the box's floor, up to the logit round trip
}

// Results do not depend on the thread count (surfaces are per point; sums are in a fixed order).
VCAL_TEST(posterior_is_identical_for_any_thread_count) {
    NormalInU obj{{logit(0.02), logit(0.2)}, {0.05, 0.1}};
    const auto a = posterior_of(obj, Backend{1});
    const auto b = posterior_of(obj, Backend{8});
    for (int k = 0; k < 2; ++k) {
        VCAL_CHECK(a.et_lo[k] == b.et_lo[k] && a.et_hi[k] == b.et_hi[k]);
        VCAL_CHECK(a.hpd_lo[k] == b.hpd_lo[k] && a.hpd_hi[k] == b.hpd_hi[k]);
        VCAL_CHECK(a.mean[k] == b.mean[k] && a.sd[k] == b.sd[k]);
    }
    VCAL_CHECK(a.flags == b.flags && a.refinements == b.refinements);
}

// The Jeffreys table for binomial-mixture periods: finite everywhere on the box and the same for any
// thread count. Its values are cross-checked independently by validation/scipy/grid_posterior.py.
VCAL_TEST(jeffreys_table_is_finite_and_deterministic) {
    using Objective = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
    const auto integrator = vcal::quadrature::parity_rule();
    const vcal::Grid<2> g{{{1e-4, 0.2, 13, vcal::AxisScale::Logit}, {1e-3, 0.5, 9, vcal::AxisScale::Logit}}};
    const auto a = e::jeffreys_table(Backend{1}, Objective{}, integrator, 100, g);
    const auto b = e::jeffreys_table(Backend{8}, Objective{}, integrator, 100, g);
    int nonfinite = 0;
    for (std::size_t k = 0; k < a.log_prior.size(); ++k) {
        nonfinite += std::isfinite(a.log_prior[k]) ? 0 : 1;
        VCAL_CHECK(a.log_prior[k] == b.log_prior[k]);
    }
    VCAL_CHECK_EQ(nonfinite, 0);
}

// The committed reference results (posterior_reference --write; replicated by
// validation/scipy/grid_posterior.py): recomputed, the rule's decisions and grids exactly, doubles to
// TOL_POSTERIOR_REFERENCE_REL.
VCAL_TEST(posterior_reference_results_reproduce) {
    int bad = 0;
    double worst = 0.0;
    for (const auto& [file, text] : {std::pair<std::string, std::string>{"posterior/jeffreys.csv", vcal::postref::jeffreys_csv()},
                                     {"posterior/posterior.csv", vcal::postref::posterior_csv()}}) {
        const auto want = vcal::test::read_golden_csv(file);
        std::vector<std::vector<std::string>> got;
        std::size_t start = 0;
        bool header = true;
        while (start < text.size()) {
            const std::size_t end = text.find('\n', start);
            const std::string line = text.substr(start, end - start);
            start = end + 1;
            if (line.empty() || line[0] == '#') continue;
            if (header) {
                header = false;
                continue;
            }
            std::vector<std::string> cells;
            std::size_t a = 0;
            for (std::size_t i = 0; i <= line.size(); ++i) {
                if (i == line.size() || line[i] == ',') {
                    cells.push_back(line.substr(a, i - a));
                    a = i + 1;
                }
            }
            got.push_back(cells);
        }
        VCAL_REQUIRE(got.size() == want.rows.size());
        for (std::size_t i = 0; i < got.size(); ++i) {
            for (std::size_t j = 0; j < want.header.size(); ++j) {
                const std::string& h = want.header[j];
                const std::string& x = want.rows[i][j];
                const std::string& y = got[i][j];
                bool ok = x == y;
                if (h.size() > 4 && h.compare(h.size() - 4, 4, "_hex") == 0) {
                    const double u = vcal::test::parse_double(x), v = vcal::test::parse_double(y);
                    if (h == "log_prior_hex") {
                        // log sqrt det I_u can be near 0, so it is compared absolutely: its central
                        // differences divide a platform's last-bit differences by 2e-4.
                        ok = std::fabs(v - u) <= tol::TOL_POSTERIOR_JEFFREYS_ABS;
                    } else {
                        const double rel = u == v ? 0.0 : std::fabs(v / u - 1.0);
                        worst = std::fmax(worst, rel);
                        ok = rel <= tol::TOL_POSTERIOR_REFERENCE_REL;
                    }
                }
                if (!ok) {
                    ++bad;
                    vcal::test::note(file + " row " + std::to_string(i) + " " + h + ": " + y + ", committed " + x);
                }
            }
        }
    }
    vcal::test::note("worst relative difference: " + sci(worst));
    VCAL_CHECK_EQ(bad, 0);
}
