// SPDX-License-Identifier: Apache-2.0
// core/objectives/vasicek_rate.hpp (M3): the Vasicek-rate log-likelihood against mpmath goldens,
// its zero-rate treatments (D-044), and the engine's grid fit against the closed-form maximum.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

// The committed recovery summary (rate_recovery_harness, R = 1,000; docs/methodology/vasicek_rate_mle.md
// section 5): every profile coverage in the band except the reviewed findings, which share their
// factor draws with the binomial recovery's small-T findings (29/PD, 68/PD; 74's q, S-23), and the
// closed form within tolerance everywhere. A verdict that moves needs a new review.
VCAL_TEST(vasicek_rate_recovery_summary_verdicts) {
    const auto g = vcal::test::read_golden_csv("vasicek_rate/recovery_summary.csv");
    VCAL_REQUIRE(g.rows.size() == 27u);
    struct Reviewed {
        const char* scenario;
        const char* column;
    };
    // The same factor draws as binomial scenarios 29, 68 and 74 (the n -> infinity limit of those
    // panels), whose profile verdicts are pinned small-T findings; S-13 put them at about 0.94 at
    // R = 10,000 (D-158). The rate model reproduces the dips on the same draws.
    const Reviewed reviewed[] = {{"29", "pd_class"}, {"68", "pd_class"}, {"74", "q_class"}};
    std::string unreviewed;
    int out_of_band = 0;
    for (const auto& r : g.rows) {
        for (const char* col : {"pd_class", "rho_class", "q_class"}) {
            const std::string c = r[g.column(col)];
            bool listed = false;
            for (const auto& v : reviewed) listed = listed || (r[g.column("z_source_scenario")] == v.scenario && std::string(col) == v.column);
            if (c != "PASS") {
                ++out_of_band;
                if (!(listed && c == "BELOW")) unreviewed += " " + r[g.column("z_source_scenario")] + "/" + col + " " + c;
            } else if (listed) {
                unreviewed += " " + r[g.column("z_source_scenario")] + "/" + col + " now PASS";
            }
        }
        VCAL_CHECK(hex(r[g.column("closed_form_gap_max_hex")]) <= tol::TOL_VASICEK_RATE_POLISH_U);
        VCAL_CHECK(hex(r[g.column("residual_max_hex")]) <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);
        VCAL_CHECK_EQ(vcal::test::parse_int(r[g.column("not_computed")]), 0);
        // Every Wald flag is near-bound (the rho cap at T = 20, rho = 0.24): no edge or flat fits.
        VCAL_CHECK_EQ(vcal::test::parse_int(r[g.column("flagged_edge")]), 0);
        VCAL_CHECK_EQ(vcal::test::parse_int(r[g.column("flagged_flat")]), 0);
    }
    if (!unreviewed.empty()) vcal::test::note("unreviewed:" + unreviewed);
    VCAL_CHECK(unreviewed.empty());
    VCAL_CHECK_EQ(out_of_band, 3);
}

namespace {

// The derivatives chained to logit coordinates: g_u = g J, H_u = J_i J_j H + delta_ij g_i J_i (1 - 2 theta_i).
template <class O>
struct LogitRate {
    double l, g[2], h[3];
};

template <class O>
LogitRate<O> logit_rate(const O& obj, const o::RateObs& y, double u, double w) {
    const auto primary = vcal::quadrature::parity_rule();
    const double pd = vcal::grid::from_scaled(vcal::AxisScale::Logit, u);
    const double rho = vcal::grid::from_scaled(vcal::AxisScale::Logit, w);
    double g[2], h[3];
    LogitRate<O> r{};
    r.l = obj.log_contrib_derivs(y, {pd, rho}, primary, g, h);
    const double J[2] = {pd * (1.0 - pd), rho * (1.0 - rho)};
    r.g[0] = g[0] * J[0];
    r.g[1] = g[1] * J[1];
    r.h[0] = h[0] * J[0] * J[0] + g[0] * J[0] * (1.0 - 2.0 * pd);
    r.h[1] = h[1] * J[1] * J[1] + g[1] * J[1] * (1.0 - 2.0 * rho);
    r.h[2] = h[2] * J[0] * J[1];
    return r;
}

template <class O>
void check_rate_derivs(const O& obj, const std::vector<o::RateObs>& ys, const char* name, double& worst_g,
                       double& worst_h, double& bits) {
    const auto primary = vcal::quadrature::parity_rule();
    const double h = 0.02;
    for (const auto& y : ys) {
        for (const double pd : {1.0001e-4, 0.001, 0.01, 0.05, 0.19999}) {
            for (const double rho : {1.0001e-3, 0.02, 0.12, 0.24, 0.49999}) {
                const double u = vcal::grid::to_scaled(vcal::AxisScale::Logit, pd);
                const double w = vcal::grid::to_scaled(vcal::AxisScale::Logit, rho);
                const auto a = logit_rate(obj, y, u, w);
                const double plain = obj.log_contrib(y, {vcal::grid::from_scaled(vcal::AxisScale::Logit, u),
                                                         vcal::grid::from_scaled(vcal::AxisScale::Logit, w)},
                                                     primary);
                if (std::memcmp(&plain, &a.l, sizeof plain) != 0) bits = 1.0;
                const double s0 = h / std::sqrt(std::fmax(1.0, std::fabs(a.h[0])));
                const double s1 = h / std::sqrt(std::fmax(1.0, std::fabs(a.h[1])));
                const auto diff = [&](double k, double (&dg)[2], double (&dh)[3]) {
                    const auto up = logit_rate(obj, y, u + k * s0, w), um = logit_rate(obj, y, u - k * s0, w);
                    const auto wp = logit_rate(obj, y, u, w + k * s1), wm = logit_rate(obj, y, u, w - k * s1);
                    dg[0] = (up.l - um.l) / (2 * k * s0);
                    dg[1] = (wp.l - wm.l) / (2 * k * s1);
                    dh[0] = (up.g[0] - um.g[0]) / (2 * k * s0);
                    dh[1] = (wp.g[1] - wm.g[1]) / (2 * k * s1);
                    dh[2] = (wp.g[0] - wm.g[0]) / (2 * k * s1);
                };
                double g1[2], h1[3], g2[2], h2[3];
                diff(1.0, g1, h1);
                diff(0.5, g2, h2);
                for (int i = 0; i < 2; ++i) {
                    const double fd = (4.0 * g2[i] - g1[i]) / 3.0;
                    worst_g = std::fmax(worst_g, std::fabs(a.g[i] - fd) / std::fmax(1.0, std::fabs(a.g[i])));
                }
                for (int i = 0; i < 3; ++i) {
                    const double fd = (4.0 * h2[i] - h1[i]) / 3.0;
                    worst_h = std::fmax(worst_h, std::fabs(a.h[i] - fd) / std::fmax(1.0, std::fabs(a.h[i])));
                }
            }
        }
    }
    vcal::test::note(std::string(name) + ": worst relative score error " + sci(worst_g) + ", Hessian " + sci(worst_h));
}

}  // namespace

// The analytic score and Hessian (D-417) against Richardson-extrapolated central differences in logit coordinates,
// each axis's step 0.02 / sqrt(max(1, |l_aa|)), relative to max(1, |analytic|), as the binomial objective's test:
// interior rates from 1e-8 to 1 - 1e-6 under every treatment, and rates of 0 and 1 (detect = 1/(2n), n = 100 and
// 1e6) under the censored and substituted treatments, over the recovery box's corners and interior. The value is
// log_contrib's bit for bit. Tolerance TOL_SCORE_HESSIAN_FD_REL, the binomial objective's.
VCAL_TEST(vasicek_rate_score_hessian_match_finite_differences) {
    std::vector<o::RateObs> interior;
    for (const double r : {1e-8, 1e-4, 0.003, 0.01, 0.05, 0.3, 0.7, 0.99, 1.0 - 1e-6}) interior.push_back({r, 0.5e-6});
    std::vector<o::RateObs> boundary;
    for (const double n : {100.0, 1e6}) {
        boundary.push_back({0.0, 0.5 / n});
        boundary.push_back({1.0, 0.5 / n});
    }
    double worst_g = 0.0, worst_h = 0.0, bits = 0.0;
    check_rate_derivs(Refuse{}, interior, "density", worst_g, worst_h, bits);
    double cg = 0.0, ch = 0.0;
    check_rate_derivs(Censor{}, boundary, "censored boundary terms", cg, ch, bits);
    double sg = 0.0, sh = 0.0;
    check_rate_derivs(Substitute{}, boundary, "substituted boundary terms", sg, sh, bits);
    VCAL_CHECK(bits == 0.0);
    for (const double e : {worst_g, worst_h, cg, ch, sg, sh}) VCAL_CHECK(e <= tol::TOL_SCORE_HESSIAN_FD_REL);
}

// At the closed-form MLE of a panel without boundary rates, the engine's analytic SE is the closed-form SE
// (D-417). x_t = Phi^-1(r_t) is a normal sample, so at its MLE the observed information in (mu, v) is
// diag(T/v, T/(2 v^2)) exactly. With logit rho = log v, SE(logit rho) = sqrt(2/T); and with
// a = mu/sqrt(1 + v), PD = Phi(a), Var(logit PD) = k^2 (v/((1 + v) T) + mu^2 v^2 / (2 (1 + v)^3 T)),
// k = phi(a)/(Phi(a)(1 - Phi(a))). Within 1e-10 relative, stated before the run.
VCAL_TEST(vasicek_rate_analytic_se_is_the_closed_form) {
    const auto primary = vcal::quadrature::parity_rule();
    const auto g = box_grid();
    double worst = 0.0;
    int compared = 0;
    for (const std::uint32_t scenario : {38u, 41u, 44u, 50u, 59u, 26u, 71u}) {
        const double pd = scenario == 59 || scenario == 71 ? 0.05 : scenario == 26 ? 0.001 : 0.01;
        const double rho = scenario == 50 || scenario == 26 ? 0.24 : scenario == 59 ? 0.02 : 0.12;
        const std::int64_t T = scenario == 38 ? 20 : scenario == 41 || scenario == 50 || scenario == 59 ? 40 : 100;
        for (std::uint32_t r = 0; r < 5; ++r) {
            const auto obs = rate_panel(scenario, r, pd, rho, T);
            o::RateClosedForm cf{};
            VCAL_REQUIRE(o::vasicek_rate_closed_form(obs.data(), T, cf));
            if (!(cf.pd > g.axis[0].lo && cf.pd < g.axis[0].hi && cf.rho > g.axis[1].lo && cf.rho < g.axis[1].hi)) continue;
            double mu = 0.0;
            for (const auto& y : obs) mu += o::detail::rate_probit(y.rate);
            mu /= static_cast<double>(T);
            const double v = cf.rho / (1.0 - cf.rho);
            const double a = mu / std::sqrt(1.0 + v);
            const double Pa = 0.5 * std::erfc(-a / std::sqrt(2.0));
            const double k = std::exp(-0.5 * a * a) * vcal::special::constants::kInvSqrt2Pi / (Pa * (1.0 - Pa));
            const double Td = static_cast<double>(T);
            const double se_pd = k * std::sqrt(v / ((1.0 + v) * Td) + mu * mu * v * v / (2.0 * std::pow(1.0 + v, 3) * Td));
            const double se_rho = std::sqrt(2.0 / Td);
            const double us[2] = {vcal::grid::to_scaled(vcal::AxisScale::Logit, cf.pd),
                                  vcal::grid::to_scaled(vcal::AxisScale::Logit, cf.rho)};
            const auto an = e::analytic_se_at_scaled(Refuse{}, primary, obs.data(), T, g, us);
            VCAL_REQUIRE(an.ok);
            const double e0 = std::fabs(an.se_u[0] / se_pd - 1.0), e1 = std::fabs(an.se_u[1] / se_rho - 1.0);
            worst = std::fmax(worst, std::fmax(e0, e1));
            VCAL_CHECK(e0 <= 1e-10 && e1 <= 1e-10);
            ++compared;
        }
    }
    VCAL_CHECK(compared >= 20);
    vcal::test::note(std::to_string(compared) + " panels: worst relative SE gap " + sci(worst));
}
