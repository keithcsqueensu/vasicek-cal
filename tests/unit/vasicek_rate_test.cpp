// SPDX-License-Identifier: Apache-2.0
// core/objectives/vasicek_rate.hpp (M3): the Vasicek-rate log-likelihood against mpmath goldens,
// its zero-rate treatments (D-044), and the engine's grid fit against the closed-form maximum.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/vasicek_rate.hpp"
#include "core/quadrature/parity.hpp"
#include "dgp/dgp.hpp"
#include "engine/calibrate.hpp"
#include "engine/profile.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace o = vcal::objectives;
namespace e = vcal::engine;
namespace tol = vcal::tol;
using vcal::PrecisionF64;
using Refuse = o::VasicekRate<PrecisionF64, o::ZeroRates::Refuse>;
using Censor = o::VasicekRate<PrecisionF64, o::ZeroRates::Censor>;
using Substitute = o::VasicekRate<PrecisionF64, o::ZeroRates::Substitute>;
constexpr double kEps = 2.220446049250313e-16;

double hex(const std::string& s) { return std::strtod(s.c_str(), nullptr); }

std::string sci(double v) {
    char b[32];
    std::snprintf(b, sizeof b, "%.3g", v);
    return b;
}

// The M1.7 box on logit axes (D-115), as the recovery harness uses.
vcal::Grid<2> box_grid() {
    return {{{1e-4, 0.2, 61, vcal::AxisScale::Logit}, {1e-3, 0.5, 41, vcal::AxisScale::Logit}}};
}

// Rates in the n -> infinity limit of a recovery panel: the DGP's factor draws z_t (the same draws
// as the recovery panel with that scenario id and replicate), mapped through the model exactly.
std::vector<o::RateObs> rate_panel(std::uint32_t scenario, std::uint32_t replicate, double pd, double rho,
                                   std::int64_t periods) {
    std::vector<std::int64_t> n(static_cast<std::size_t>(periods), 1), d(n.size());
    std::vector<double> z(n.size());
    const vcal::dgp::PanelSpec spec{0x4D31385245434F56ull, scenario, replicate, pd, rho, n.data(), periods};
    if (vcal::dgp::simulate_panel(spec, d.data(), z.data()) != vcal::dgp::Status::Ok) std::abort();
    const double c = vcal::special::probit(pd);
    std::vector<o::RateObs> obs(n.size());
    for (std::size_t t = 0; t < obs.size(); ++t) {
        const double x = (c - std::sqrt(rho) * z[t]) / std::sqrt(1.0 - rho);
        obs[t] = {0.5 * std::erfc(-x / std::sqrt(2.0)), 0.5e-6};
    }
    return obs;
}

}  // namespace

VCAL_TEST(vasicek_rate_log_density_matches_mpmath) {
    const auto g = vcal::test::read_golden_csv("vasicek_rate/log_density.csv");
    VCAL_REQUIRE(g.rows.size() == 250u);
    const auto integrator = vcal::quadrature::parity_rule();
    double worst = 0.0;
    for (const auto& r : g.rows) {
        const Refuse::Theta th{hex(r[g.column("pd_hex")]), hex(r[g.column("rho_hex")])};
        const double got = Refuse{}.log_contrib({hex(r[g.column("rate_hex")]), 0.25}, th, integrator);
        const double want = hex(r[g.column("log_density_hex")]);
        const double err = std::fabs(got - want) / (kEps * hex(r[g.column("scale_hex")]));
        worst = std::fmax(worst, err);
        VCAL_CHECK(err <= tol::TOL_VASICEK_RATE_LOGLIK_EPS);
    }
    vcal::test::note("worst |error| / (eps x term scale): " + sci(worst));
}

VCAL_TEST(vasicek_rate_boundary_treatments_match_mpmath) {
    const auto g = vcal::test::read_golden_csv("vasicek_rate/boundary.csv");
    VCAL_REQUIRE(g.rows.size() == 200u);
    const auto integrator = vcal::quadrature::parity_rule();
    double worst = 0.0;
    for (const auto& r : g.rows) {
        const Censor::Theta th{hex(r[g.column("pd_hex")]), hex(r[g.column("rho_hex")])};
        const o::RateObs y{r[g.column("rate")] == "0" ? 0.0 : 1.0, hex(r[g.column("detect_hex")])};
        const bool censor = r[g.column("treatment")] == "censor";
        const double got = censor ? Censor{}.log_contrib(y, th, integrator) : Substitute{}.log_contrib(y, th, integrator);
        const double want = hex(r[g.column("value_hex")]);
        const double err = std::fabs(got - want) / (kEps * hex(r[g.column("scale_hex")]));
        worst = std::fmax(worst, err);
        VCAL_CHECK(err <= tol::TOL_VASICEK_RATE_LOGLIK_EPS);
    }
    vcal::test::note("worst |error| / (eps x term scale): " + sci(worst));
}

VCAL_TEST(vasicek_rate_parity_refuses_boundary_rates_and_names_them) {
    const std::vector<o::RateObs> obs = {{0.01, 5e-3}, {0.0, 5e-3}, {0.02, 5e-3}, {1.0, 5e-3}};
    VCAL_CHECK(Refuse::panel_error(obs.data(), 4) != nullptr);
    VCAL_CHECK(Censor::panel_error(obs.data(), 4) == nullptr);
    VCAL_CHECK(Substitute::panel_error(obs.data(), 4) == nullptr);
    const auto named = o::zero_rate_periods(obs.data(), 4);
    VCAL_REQUIRE(named.size() == 2u);
    VCAL_CHECK(named[0] == 1 && named[1] == 3);
    const auto kept = o::drop_zero_rate_periods(obs.data(), 4);
    VCAL_REQUIRE(kept.size() == 2u);
    VCAL_CHECK(kept[0] == obs[0] && kept[1] == obs[2]);
    // Count data: rate d/n and a detection limit of half a default.
    const auto y = o::rate_obs(200, 3);
    VCAL_CHECK(y.rate == 3.0 / 200.0 && y.detect == 0.5 / 200.0);
    // The engine refuses the panel as a whole (D-044).
    const auto primary = vcal::quadrature::parity_rule();
    const auto check = vcal::quadrature::parity_rule(true);
    std::vector<double> L;
    e::Estimate2 est{};
    VCAL_CHECK(e::calibrate(vcal::backends::CpuBackend{1}, Refuse{}, primary, check, obs.data(), 4, box_grid(), L,
                            est) == e::Status::InvalidPanel);
    // Out-of-range inputs are panel errors under every treatment.
    const std::vector<o::RateObs> bad = {{-0.1, 5e-3}, {0.01, 0.0}};
    VCAL_CHECK(Censor::panel_error(bad.data(), 1) != nullptr);
    VCAL_CHECK(Censor::panel_error(bad.data() + 1, 1) != nullptr);
}

// The grid path (surface, argmax, refinement; the profile code's polished maximum) against the
// closed-form maximum, on rate panels from the recovery DGP's factor draws.
VCAL_TEST(vasicek_rate_grid_fit_matches_the_closed_form) {
    struct Case {
        std::uint32_t scenario;
        double pd, rho;
        std::int64_t periods;
    };
    const Case cases[] = {{38, 0.01, 0.12, 20}, {41, 0.01, 0.12, 40}, {44, 0.01, 0.12, 100}, {50, 0.01, 0.24, 40},
                          {59, 0.05, 0.02, 40}, {26, 0.001, 0.24, 100}, {71, 0.05, 0.12, 100}};
    const auto primary = vcal::quadrature::parity_rule();
    const auto check = vcal::quadrature::parity_rule(true);
    const auto g = box_grid();
    double worst_u = 0.0, worst_ll = 0.0;
    for (const auto& c : cases) {
        for (std::uint32_t r = 0; r < 5; ++r) {
            const auto obs = rate_panel(c.scenario, r, c.pd, c.rho, c.periods);
            o::RateClosedForm cf{};
            VCAL_REQUIRE(o::vasicek_rate_closed_form(obs.data(), c.periods, cf));
            const bool inside = cf.pd > g.axis[0].lo && cf.pd < g.axis[0].hi && cf.rho > g.axis[1].lo && cf.rho < g.axis[1].hi;
            if (!inside) continue;  // the grid path is confined to the box; the closed form is not
            std::vector<double> L;
            e::Estimate2 est{};
            VCAL_REQUIRE(e::calibrate(vcal::backends::CpuBackend{1}, Refuse{}, primary, check, obs.data(), c.periods, g,
                                      L, est) == e::Status::Ok);
            VCAL_CHECK(est.quad_check_flagged == 0);
            const auto prof = e::profile_intervals(Refuse{}, primary, obs.data(), c.periods, g, L, est);
            const double du[2] = {std::fabs(vcal::grid::to_scaled(vcal::AxisScale::Logit, prof.max_at[0]) -
                                            vcal::grid::to_scaled(vcal::AxisScale::Logit, cf.pd)),
                                  std::fabs(vcal::grid::to_scaled(vcal::AxisScale::Logit, prof.max_at[1]) -
                                            vcal::grid::to_scaled(vcal::AxisScale::Logit, cf.rho))};
            double ll_cf = 0.0;
            for (const auto& y : obs) ll_cf += Refuse{}.log_contrib(y, {cf.pd, cf.rho}, primary);
            worst_u = std::fmax(worst_u, std::fmax(du[0], du[1]));
            worst_ll = std::fmax(worst_ll, ll_cf - prof.loglik_max);
            VCAL_CHECK(du[0] <= tol::TOL_VASICEK_RATE_POLISH_U && du[1] <= tol::TOL_VASICEK_RATE_POLISH_U);
            VCAL_CHECK(ll_cf - prof.loglik_max <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);
            VCAL_CHECK(prof.loglik_max - ll_cf <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);  // cannot beat the exact max
        }
    }
    vcal::test::note("worst polished-vs-closed-form distance in logit units: " + sci(worst_u) +
                     "; worst log-likelihood shortfall: " + sci(worst_ll));
}
