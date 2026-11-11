// SPDX-License-Identifier: Apache-2.0
// The multi-grade shared-rho MLE (M6; D-411, D-413; benchmark_models/multigrade/fit.hpp), with tolerances stated
// before the tests ran (D-413):
//   - at K = 1 it is calibrate's estimate, SEs and rho profile interval, within solver tolerance;
//   - two identical grades give one PD, and the pooled single-grade fit (pooling by symmetry);
//   - a grade with no defaults lands on the PD box's lower end, flagged, with rho barely moved;
//   - an empty grade changes nothing, bit for bit; the fit is deterministic;
//   - the estimate is a maximum of ref's independent K-grade likelihood;
//   - recovery on shared-factor panels, K in {2, 5, 17}, T in {20, 100}, rho in {0.24, 0.12, 0.02} (slow).
// Panels: until PR 3's DGP layer (D-151), the test draws Z_t with dgp::factor_draw and each grade's defaults given
// Z_t with dgp::count_defaults, grade k on scenario (S << 6) + k, so the grades' Bernoulli streams are disjoint.
// Seed 0x4D474649540A ("MGFIT"): fresh, disjoint from the other tests' and studies' seeds.
#include "benchmark_models/multigrade/fit.hpp"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "core/special/lbinom.hpp"
#include "dgp/dgp.hpp"
#include "engine/calibrate.hpp"
#include "engine/profile.hpp"
#include "ref/ref.hpp"
#include "tests/harness/golden.hpp"
#include "tests/recovery/grade_matrix.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace mg = vcal::benchmark_models::multigrade;
using mg::GradeCount;
using vcal::test::describe;

constexpr std::uint64_t kSeed = 0x4D474649540Aull;
const mg::Box kBox{};

struct Panel {
    std::vector<GradeCount> y;  // T x K
    int K = 0;
    std::string name;
    std::int64_t T() const { return static_cast<std::int64_t>(y.size()) / K; }
};

// Shared-factor panel: per period one Z_t, then each grade's defaults given Z_t.
Panel shared_panel(std::uint32_t scenario, std::uint32_t replicate, const std::vector<double>& pd,
                   const std::vector<std::int64_t>& n, double rho, std::int64_t T) {
    Panel p;
    p.K = static_cast<int>(pd.size());
    const auto kernel = vcal::dgp::bernoulli_kernel_selected();
    for (std::int64_t t = 0; t < T; ++t) {
        const auto period = static_cast<std::uint32_t>(t);
        double z = 0.0;
        VCAL_REQUIRE(vcal::dgp::factor_draw(kSeed, scenario << 6, replicate, period, vcal::dgp::kFactorBlocks, &z) ==
                     vcal::dgp::Status::Ok);
        for (std::size_t k = 0; k < pd.size(); ++k) {
            const double c = vcal::dgp::det_probit(pd[k]);
            const double pk = vcal::dgp::det_ncdf((c - std::sqrt(rho) * z) / std::sqrt(1.0 - rho));
            const std::int64_t d = vcal::dgp::count_defaults(kSeed, (scenario << 6) + static_cast<std::uint32_t>(k),
                                                             replicate, period, n[k], pk, kernel);
            VCAL_REQUIRE(d >= 0);
            p.y.push_back({n[k], d});
        }
    }
    p.name = "scenario " + std::to_string(scenario) + "/" + std::to_string(replicate) + " (K " +
             std::to_string(p.K) + ", rho " + describe(rho) + ", T " + std::to_string(T) + ")";
    return p;
}

mg::Fit fit(const Panel& p, const mg::Box& box = kBox) {
    mg::Fit f;
    VCAL_REQUIRE(mg::fit(p.y.data(), p.T(), p.K, box, f) == mg::Status::Ok);
    return f;
}

double logit(double x) { return mg::detail::logit(x); }
bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

}  // namespace

// K = 1: the estimate is calibrate's (a 61 x 41 grid, refined) as polished by profile_intervals, within 1e-5 in
// each logit coordinate; the log-likelihoods within 1e-9; the rho profile interval's ends within 1e-6 in logit rho.
// Different searches on the same likelihood (K = 1 is the single-grade one bit for bit, D-412). The SEs within 1e-8
// relative of the engine's analytic SE at the same point (D-185); calibrate's own SE, central differences of the
// surface (D-119), differs by up to 2% and is reported only (D-413: the stated 1e-3 against it was the wrong
// comparand).
VCAL_TEST(multigrade_fit_at_one_grade_is_calibrate) {
    using Objective = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
    const auto primary = vcal::quadrature::parity_rule();
    const auto check = vcal::quadrature::parity_rule(true);
    const vcal::Grid<2> g{{{1e-4, 0.2, 61, vcal::AxisScale::Logit}, {1e-3, 0.5, 41, vcal::AxisScale::Logit}}};
    std::uint32_t scenario = 1;
    for (const auto& spec : {std::vector<double>{0.01, 0.12, 20, 1000}, std::vector<double>{0.05, 0.24, 40, 100},
                             std::vector<double>{0.001, 0.12, 20, 10000}, std::vector<double>{0.02, 0.03, 30, 500}}) {
        const auto p = shared_panel(scenario++, 0, {spec[0]}, {static_cast<std::int64_t>(spec[3])}, spec[1],
                                    static_cast<std::int64_t>(spec[2]));
        std::vector<Objective::Obs> obs;
        for (const auto& o : p.y) obs.push_back({o.n, o.d});
        std::vector<double> L;
        vcal::engine::Estimate2 est{};
        VCAL_REQUIRE(vcal::engine::calibrate(vcal::backends::CpuBackend{1}, Objective{}, primary, check, obs.data(),
                                             p.T(), g, L, est) == vcal::engine::Status::Ok);
        VCAL_REQUIRE(est.flags == 0u);
        const auto prof = vcal::engine::profile_intervals(Objective{}, primary, obs.data(), p.T(), g, L, est);
        const auto f = fit(p);
        VCAL_CHECK_EQ(f.flags, 0u);
        const double gap_pd = std::fabs(logit(f.pd[0]) - logit(prof.max_at[0]));
        const double gap_rho = std::fabs(logit(f.rho) - logit(prof.max_at[1]));
        const double gap_ll = std::fabs(f.loglik - prof.loglik_max);
        const double gap_lo = std::fabs(logit(f.rho_lo) - logit(prof.lo[1]));
        const double gap_hi = std::fabs(logit(f.rho_hi) - logit(prof.hi[1]));
        // The SEs against the engine's analytic SE at this estimate (the same computation in another code path);
        // calibrate's own SE is central differences (D-119), reported only.
        const double us[2] = {logit(f.pd[0]), logit(f.rho)};
        const auto an = vcal::engine::analytic_se_at_scaled(Objective{}, primary, obs.data(), p.T(), g, us);
        VCAL_REQUIRE(an.ok);
        const double se_pd = std::fabs(f.pd_se[0] / (f.pd[0] * (1.0 - f.pd[0])) / an.se_u[0] - 1.0);
        const double se_rho = std::fabs(f.rho_se / (f.rho * (1.0 - f.rho)) / an.se_u[1] - 1.0);
        const double cd_pd = f.pd_se[0] / est.se[0] - 1.0, cd_rho = f.rho_se / est.se[1] - 1.0;
        VCAL_CHECK(gap_pd <= 1e-5);
        VCAL_CHECK(gap_rho <= 1e-5);
        VCAL_CHECK(gap_ll <= 1e-9);
        VCAL_CHECK(gap_lo <= 1e-6);
        VCAL_CHECK(gap_hi <= 1e-6);
        VCAL_CHECK(se_pd <= 1e-8);
        VCAL_CHECK(se_rho <= 1e-8);
        VCAL_CHECK_EQ(f.rho_interval_flags, prof.flags[1]);
        vcal::test::note(p.name + ": logit gaps " + describe(gap_pd) + ", " + describe(gap_rho) + "; loglik " +
                         describe(gap_ll) + "; interval ends " + describe(gap_lo) + ", " + describe(gap_hi) +
                         "; SEs " + describe(se_pd) + ", " + describe(se_rho) + " relative (calibrate's differenced " +
                         describe(cd_pd) + ", " + describe(cd_rho) + "); " +
                         std::to_string(f.evaluations) + " + " + std::to_string(f.derivative_evaluations) +
                         " evaluations");
    }
}

// Pooling by symmetry: a grade duplicated (each period's counts twice) gives one PD, |u_1 - u_2| <= 1e-8, and the
// K = 1 fit on the doubled counts: u and v within 1e-6, the interval's ends within 1e-6, and the log-likelihoods equal
// once their binomial coefficients are removed, within 1e-9. (D-411's plan asked for grades generated with one PD and
// PD_1 ~ PD_2; that holds only to sampling error, so the exact identity replaces it, D-413.)
VCAL_TEST(multigrade_fit_pools_a_duplicated_grade) {
    std::uint32_t scenario = 10;
    for (const auto& spec : {std::vector<double>{0.01, 0.12, 20, 1000}, std::vector<double>{0.03, 0.2, 30, 300}}) {
        const auto base = shared_panel(scenario++, 0, {spec[0]}, {static_cast<std::int64_t>(spec[3])}, spec[1],
                                       static_cast<std::int64_t>(spec[2]));
        Panel two, pooled;
        two.K = 2;
        pooled.K = 1;
        double lc_two = 0.0, lc_pooled = 0.0;
        for (const auto& o : base.y) {
            two.y.push_back(o);
            two.y.push_back(o);
            pooled.y.push_back({2 * o.n, 2 * o.d});
            lc_two += 2.0 * vcal::special::lbinom(o.n, o.d);
            lc_pooled += vcal::special::lbinom(2 * o.n, 2 * o.d);
        }
        const auto f2 = fit(two);
        const auto f1 = fit(pooled);
        VCAL_CHECK_EQ(f2.flags, 0u);
        VCAL_CHECK_EQ(f1.flags, 0u);
        const double split = std::fabs(logit(f2.pd[0]) - logit(f2.pd[1]));
        const double gap_pd = std::fmax(std::fabs(logit(f2.pd[0]) - logit(f1.pd[0])),
                                        std::fabs(logit(f2.pd[1]) - logit(f1.pd[0])));
        const double gap_rho = std::fabs(logit(f2.rho) - logit(f1.rho));
        const double gap_ll = std::fabs((f2.loglik - lc_two) - (f1.loglik - lc_pooled));
        const double gap_end = std::fmax(std::fabs(logit(f2.rho_lo) - logit(f1.rho_lo)),
                                         std::fabs(logit(f2.rho_hi) - logit(f1.rho_hi)));
        VCAL_CHECK(split <= 1e-8);
        VCAL_CHECK(gap_pd <= 1e-6);
        VCAL_CHECK(gap_rho <= 1e-6);
        VCAL_CHECK(gap_ll <= 1e-9);
        VCAL_CHECK(gap_end <= 1e-6);
        vcal::test::note(base.name + ": split " + describe(split) + "; gaps " + describe(gap_pd) + ", " +
                         describe(gap_rho) + "; loglik " + describe(gap_ll) + "; ends " + describe(gap_end));
    }
}

// A grade with no defaults in any period: its PD is exactly the box's lower end, kGradeAtLower and kFlagGradeEdge
// are set, its SE and Wald ends are NaN, and nothing else is flagged. rho is estimated from every grade; its move
// from the fit without the grade is within 0.25 of rho's SE in v.
VCAL_TEST(multigrade_fit_puts_a_zero_default_grade_on_the_edge) {
    std::uint32_t scenario = 20;
    for (const std::int64_t n0 : {std::int64_t{100}, std::int64_t{2000}}) {
        auto p = shared_panel(scenario++, 0, {0.01, 0.04}, {1000, 500}, 0.12, 20);
        Panel with;
        with.K = 3;
        with.name = p.name + " + a zero-default grade of n " + std::to_string(n0);
        for (std::int64_t t = 0; t < p.T(); ++t) {
            with.y.push_back(p.y[static_cast<std::size_t>(2 * t)]);
            with.y.push_back({n0, 0});
            with.y.push_back(p.y[static_cast<std::size_t>(2 * t + 1)]);
        }
        const auto f = fit(with);
        const auto f0 = fit(p);
        VCAL_CHECK(bits_equal(f.pd[1], kBox.pd_lo));
        VCAL_CHECK_EQ(f.grade_flags[1], mg::kGradeAtLower);
        VCAL_CHECK_EQ(f.grade_flags[0], 0u);
        VCAL_CHECK_EQ(f.grade_flags[2], 0u);
        VCAL_CHECK_EQ(f.flags, mg::kFlagGradeEdge);
        VCAL_CHECK(std::isnan(f.pd_se[1]) && std::isnan(f.pd_lo[1]) && std::isnan(f.pd_hi[1]));
        VCAL_CHECK(std::isfinite(f.pd_se[0]) && std::isfinite(f.pd_se[2]) && std::isfinite(f.rho_se));
        const double se_v = f.rho_se / (f.rho * (1.0 - f.rho));
        const double move = std::fabs(logit(f.rho) - logit(f0.rho)) / se_v;
        VCAL_CHECK(move <= 0.25);
        vcal::test::note(with.name + ": rho " + describe(f.rho) + " against " + describe(f0.rho) + " without it, " +
                         describe(move) + " SE");
    }
}

// An empty grade (n = 0 in every period), inserted first, in the middle and last, changes nothing bit for bit; the
// fit is deterministic (two runs, bit for bit).
VCAL_TEST(multigrade_fit_ignores_an_empty_grade_and_is_deterministic) {
    const auto p = shared_panel(30, 0, {0.005, 0.02, 0.08}, {2000, 1000, 300}, 0.15, 25);
    const auto f = fit(p);
    const auto again = fit(p);
    VCAL_CHECK(bits_equal(f.loglik, again.loglik) && bits_equal(f.rho, again.rho) &&
               bits_equal(f.rho_lo, again.rho_lo) && bits_equal(f.rho_hi, again.rho_hi));
    for (int at = 0; at <= p.K; ++at) {
        Panel q;
        q.K = p.K + 1;
        for (std::int64_t t = 0; t < p.T(); ++t) {
            for (int k = 0; k <= p.K; ++k) {
                if (k == at) {
                    q.y.push_back({0, 0});
                } else {
                    q.y.push_back(p.y[static_cast<std::size_t>(t * p.K + (k < at ? k : k - 1))]);
                }
            }
        }
        const auto e = fit(q);
        VCAL_CHECK_EQ(e.grade_flags[static_cast<std::size_t>(at)], mg::kGradeEmpty);
        VCAL_CHECK(std::isnan(e.pd[static_cast<std::size_t>(at)]));
        bool same = bits_equal(e.loglik, f.loglik) && bits_equal(e.rho, f.rho) && bits_equal(e.rho_se, f.rho_se) &&
                    bits_equal(e.rho_lo, f.rho_lo) && bits_equal(e.rho_hi, f.rho_hi) && e.flags == f.flags;
        for (int k = 0; k < p.K; ++k) {
            const auto ke = static_cast<std::size_t>(k < at ? k : k + 1);
            const auto kf = static_cast<std::size_t>(k);
            same = same && bits_equal(e.pd[ke], f.pd[kf]) && bits_equal(e.pd_se[ke], f.pd_se[kf]) &&
                   bits_equal(e.pd_lo[ke], f.pd_lo[kf]) && bits_equal(e.pd_hi[ke], f.pd_hi[kf]);
        }
        VCAL_CHECK(same);
    }
}

// The estimate is a maximum of ref's independent likelihood (no quadrature code in common, D-031): every step of
// 1e-3 in one logit coordinate, and along each pair's two diagonals, lowers ref's panel log-likelihood. Interior
// estimates only.
VCAL_TEST(multigrade_fit_is_a_maximum_of_the_reference) {
    const auto p = shared_panel(40, 0, {0.004, 0.015, 0.06}, {3000, 1500, 400}, 0.1, 20);
    const auto f = fit(p);
    VCAL_REQUIRE(f.flags == 0u);
    const int n = p.K + 1;
    std::vector<double> w(static_cast<std::size_t>(n));
    for (int k = 0; k < p.K; ++k) w[static_cast<std::size_t>(k)] = logit(f.pd[static_cast<std::size_t>(k)]);
    w[static_cast<std::size_t>(p.K)] = logit(f.rho);
    const auto ref_ll = [&](const std::vector<double>& x) {
        std::vector<double> pd(static_cast<std::size_t>(p.K));
        for (int k = 0; k < p.K; ++k) pd[static_cast<std::size_t>(k)] = mg::detail::inv_logit(x[static_cast<std::size_t>(k)]);
        const double rho = mg::detail::inv_logit(x[static_cast<std::size_t>(p.K)]);
        double l = 0.0;
        for (std::int64_t t = 0; t < p.T(); ++t) {
            std::vector<vcalref::Period> g;
            for (int k = 0; k < p.K; ++k) {
                const auto& o = p.y[static_cast<std::size_t>(t * p.K + k)];
                g.push_back({o.n, o.d});
            }
            l += vcalref::period_loglik_grades(pd, rho, g);
        }
        return l;
    };
    const double at = ref_ll(w);
    double least = HUGE_VAL;
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            for (const double si : {-1.0, 1.0}) {
                for (const double sj : {-1.0, 1.0}) {
                    if (i == j && sj < 0.0) continue;
                    auto x = w;
                    x[static_cast<std::size_t>(i)] += si * 1e-3;
                    if (j != i) x[static_cast<std::size_t>(j)] += sj * 1e-3;
                    const double drop = at - ref_ll(x);
                    least = std::fmin(least, drop);
                    VCAL_CHECK(drop > 0.0);
                }
            }
        }
    }
    vcal::test::note(p.name + ": least drop of ref " + describe(least) + "; ref minus core at the estimate " +
                     describe(at - f.loglik));
}

// The two committed panels of tests/golden/multigrade/estimates.csv (tests/studies/multigrade_goldens.cpp; the
// panels validation/scipy/multigrade_mle.py maximises on its own): the refit's PDs, rho and rho interval ends
// within TOL_MULTIGRADE_FIT_REGRESSION_LOGIT of the committed values in logit units, its log-likelihood within the
// same relative to max(1, |l|), and the same flags.
VCAL_TEST(multigrade_fit_matches_the_committed_estimates) {
    const auto table = vcal::test::read_golden_csv("multigrade/estimates.csv");
    VCAL_REQUIRE(table.rows.size() == 2u);
    const auto col = [&](const char* c) { return table.column(c); };
    for (const auto& row : table.rows) {
        Panel p;
        p.K = std::stoi(row[col("K")]);
        p.name = row[col("name")];
        const std::string& counts = row[col("counts")];
        std::size_t at = 0;
        while (at < counts.size()) {
            const std::size_t colon = counts.find(':', at);
            std::size_t end = counts.find_first_of(";|", colon);
            if (end == std::string::npos) end = counts.size();
            p.y.push_back({std::stoll(counts.substr(at, colon - at)), std::stoll(counts.substr(colon + 1, end - colon - 1))});
            at = end + 1;
        }
        VCAL_REQUIRE(p.T() == std::stoll(row[col("T")]));
        const auto f = fit(p);
        std::vector<double> pd;
        const std::string& ph = row[col("pd_hex")];
        for (std::size_t a = 0; a <= ph.size();) {
            std::size_t b = ph.find(';', a);
            if (b == std::string::npos) b = ph.size();
            pd.push_back(vcal::test::parse_double(ph.substr(a, b - a)));
            a = b + 1;
        }
        VCAL_REQUIRE(static_cast<int>(pd.size()) == p.K);
        double worst = 0.0;
        for (int k = 0; k < p.K; ++k) {
            worst = std::fmax(worst, std::fabs(logit(f.pd[static_cast<std::size_t>(k)]) - logit(pd[static_cast<std::size_t>(k)])));
        }
        for (const auto& [name, got] : {std::pair<const char*, double>{"rho_hex", f.rho}, {"rho_lo_hex", f.rho_lo},
                                        {"rho_hi_hex", f.rho_hi}}) {
            worst = std::fmax(worst, std::fabs(logit(got) - logit(vcal::test::parse_double(row[col(name)]))));
        }
        const double ll = vcal::test::parse_double(row[col("loglik_hex")]);
        const double ll_gap = std::fabs(f.loglik - ll) / std::fmax(1.0, std::fabs(ll));
        VCAL_CHECK(worst <= vcal::tol::TOL_MULTIGRADE_FIT_REGRESSION_LOGIT);
        VCAL_CHECK(ll_gap <= vcal::tol::TOL_MULTIGRADE_FIT_REGRESSION_LOGIT);
        VCAL_CHECK_EQ(static_cast<long long>(f.flags), std::stoll(row[col("flags")]));
        vcal::test::note(p.name + ": worst logit gap " + describe(worst) + ", loglik " + describe(ll_gap));
    }
}

// Invalid input is refused before any work.
VCAL_TEST(multigrade_fit_refuses_invalid_input) {
    const std::vector<GradeCount> y = {{100, 2}, {50, 1}};
    mg::Fit f;
    VCAL_CHECK(mg::fit(y.data(), 1, 2, mg::Box{0.0, 0.5, 1e-3, 0.5}, f) == mg::Status::InvalidBox);
    VCAL_CHECK(mg::fit(y.data(), 1, 2, mg::Box{1e-5, 0.5, 0.6, 0.5}, f) == mg::Status::InvalidBox);
    VCAL_CHECK(mg::fit(y.data(), 0, 2, kBox, f) == mg::Status::InvalidPanel);
    VCAL_CHECK(mg::fit(y.data(), 1, 0, kBox, f) == mg::Status::InvalidPanel);
    const std::vector<GradeCount> bad = {{100, 101}, {50, 1}};
    VCAL_CHECK(mg::fit(bad.data(), 1, 2, kBox, f) == mg::Status::InvalidPanel);
}

namespace {

struct Scale {
    std::vector<double> pd;
    std::vector<std::int64_t> n;
    std::string name;
};

// The recovery panels' grade scales: two and five grades of 2,000 to 1,000 obligors, and D-151's 17-grade master
// scale, PD_k = 0.03% x 1.5^(k-1), with 10,000 obligors spread by a discretised normal of standard deviation 4
// grades centred on grade 9.
std::vector<Scale> recovery_scales() {
    Scale s17;
    s17.name = "17-grade master scale";
    double wsum = 0.0;
    std::vector<double> w;
    for (int k = 0; k < 17; ++k) {
        s17.pd.push_back(0.0003 * std::pow(1.5, k));
        w.push_back(std::exp(-0.5 * (k - 8) * (k - 8) / 16.0));
        wsum += w.back();
    }
    for (int k = 0; k < 17; ++k) s17.n.push_back(static_cast<std::int64_t>(std::llround(10000.0 * w[static_cast<std::size_t>(k)] / wsum)));
    return {{{0.005, 0.03}, {2000, 1000}, "2 grades"},
            {{0.001, 0.003, 0.01, 0.03, 0.1}, {1000, 1000, 1000, 1000, 1000}, "5 grades"},
            s17};
}

// R, 100 unless VCAL_MG_RECOVERY_R sets it; VCAL_MG_RECOVERY_SCALE (0, 1 or 2) and VCAL_MG_RECOVERY_RHO (0, 1 or 2)
// run one scale or one rho only, so the matrix
// can be split across processes. About 70 CPU-minutes at R = 100, nine tenths of it the 17-grade cells.
int recovery_replicates() {
    const std::string e = vcal::test::env_or_empty("VCAL_MG_RECOVERY_R");
    return !e.empty() && std::atoi(e.c_str()) > 0 ? std::atoi(e.c_str()) : 100;
}
int recovery_rho_only() {  // VCAL_MG_RECOVERY_RHO (0, 1 or 2): one rho only
    const std::string e = vcal::test::env_or_empty("VCAL_MG_RECOVERY_RHO");
    return e.empty() ? -1 : std::atoi(e.c_str());
}
int recovery_scale_only() {
    const std::string e = vcal::test::env_or_empty("VCAL_MG_RECOVERY_SCALE");
    return e.empty() ? -1 : std::atoi(e.c_str());
}

}  // namespace

namespace {

// One recovery cell (D-413's assertions): R panels from make(r), fitted, and rho's profile coverage, the bias
// of v, the pooled PD Wald coverage (grades with n PD T >= 20) and the flags checked and reported.
template <class Make>
void recovery_cell(const std::string& cell, const std::vector<double>& pd, const std::vector<std::int64_t>& n,
                   double rho, std::int64_t T, int R, Make make) {
    const double band = 0.95 - vcal::tol::TOL_RECOVERY_COVERAGE_BAND_Z * std::sqrt(0.95 * 0.05 / R);
    int covered = 0, flat = 0, bad = 0, edge = 0, grade_edge = 0, wald_n = 0, wald_covered = 0;
    double sum = 0.0, sum2 = 0.0;
    std::int64_t evals = 0;
    for (int r = 0; r < R; ++r) {
        const Panel p = make(static_cast<std::uint32_t>(r));
        const auto f = fit(p);
        if (f.rho_lo <= rho && rho <= f.rho_hi) ++covered;
        if (f.flags & vcal::engine::kFlagFlatSurface) ++flat;
        if (f.flags & (mg::kFlagNotConverged | vcal::engine::kFlagNumeric)) ++bad;
        if (f.flags & vcal::engine::kFlagGridEdge) ++edge;
        if (f.flags & mg::kFlagGradeEdge) ++grade_edge;
        const double e = logit(f.rho) - logit(rho);
        sum += e;
        sum2 += e * e;
        evals += f.evaluations + f.derivative_evaluations;
        for (std::size_t k = 0; k < pd.size(); ++k) {
            if (static_cast<double>(n[k]) * pd[k] * static_cast<double>(T) < 20.0) continue;
            ++wald_n;
            if (f.pd_lo[k] <= pd[k] && pd[k] <= f.pd_hi[k]) ++wald_covered;
        }
    }
    const double cov = static_cast<double>(covered) / R;
    const double mean = sum / R;
    const double sd = std::sqrt(std::fmax(0.0, sum2 / R - mean * mean));
    const double wald = wald_n > 0 ? static_cast<double>(wald_covered) / wald_n : 1.0;
    if (!(cov >= band)) vcal::test::report_failure(__FILE__, __LINE__, cell + ": profile coverage " + describe(cov));
    if (!(std::fabs(mean) <= 0.5 * sd + vcal::tol::TOL_RECOVERY_COVERAGE_BAND_Z * sd / std::sqrt(R))) {
        vcal::test::report_failure(__FILE__, __LINE__, cell + ": bias " + describe(mean) + ", sd " + describe(sd));
    }
    if (!(wald >= band)) vcal::test::report_failure(__FILE__, __LINE__, cell + ": Wald coverage " + describe(wald));
    VCAL_CHECK_EQ(bad, 0);
    VCAL_CHECK(flat <= R / 20);
    vcal::test::note(cell + ": profile coverage " + describe(cov) + ", bias " + describe(mean) + " (sd " +
                     describe(sd) + "), Wald " + describe(wald) + " over " + std::to_string(wald_n) +
                     ", rho on the box " + std::to_string(edge) + ", a grade on the box " +
                     std::to_string(grade_edge) + ", flat " + std::to_string(flat) + ", " +
                     std::to_string(evals / R) + " evaluations a fit");
}

// Replicate r of grade scenario id as a Panel.
Panel grade_panel(int id, std::uint32_t r) {
    namespace gm = vcal::recovery::grades;
    std::vector<std::int64_t> n, d;
    VCAL_REQUIRE(gm::panel(id, r, n, d));
    Panel p;
    p.K = gm::kGrades;
    for (std::size_t i = 0; i < n.size(); ++i) p.y.push_back({n[i], d[i]});
    const auto sc = gm::scenario(id);
    p.name = std::string("grade scenario ") + std::to_string(id) + " (" + sc.scale->name + ", T " +
             std::to_string(sc.periods) + ")/" + std::to_string(r);
    return p;
}

std::vector<double> master_scale() {
    std::vector<double> pd;
    for (int k = 1; k <= vcal::recovery::grades::kGrades; ++k) pd.push_back(vcal::recovery::grades::master_pd(k));
    return pd;
}

}  // namespace

// Recovery (slow), tolerances stated before the run (D-413): 3 scales x T in {20, 100} x rho in {0.24, 0.12, 0.02},
// R replicates each. Per cell: rho's profile coverage at least 0.95 - z sqrt(0.0475 / R) (z =
// TOL_RECOVERY_COVERAGE_BAND_Z); the bias of v = logit rho under half its spread plus its Monte Carlo error; the
// PDs' Wald coverage, pooled over grades with n PD T >= 20 expected defaults, in the same band; never
// kFlagNotConverged or kFlagNumeric, and kFlagFlatSurface on at most 5% of fits.
VCAL_TEST(multigrade_fit_recovers_parameters) {
    const int R = recovery_replicates();
    const auto scales = recovery_scales();
    for (int si = 0; si < static_cast<int>(scales.size()); ++si) {
        const auto& scale = scales[static_cast<std::size_t>(si)];
        std::uint32_t scenario = 100 + 6 * static_cast<std::uint32_t>(si);
        if (recovery_scale_only() >= 0 && recovery_scale_only() != si) continue;
        for (const std::int64_t T : {std::int64_t{20}, std::int64_t{100}}) {
            int ri = -1;
            for (const double rho : {0.24, 0.12, 0.02}) {
                const std::uint32_t sc = scenario++;
                ++ri;
                if (recovery_rho_only() >= 0 && recovery_rho_only() != ri) continue;
                recovery_cell(scale.name + ", T " + std::to_string(T) + ", rho " + describe(rho), scale.pd, scale.n,
                              rho, T, R, [&](std::uint32_t r) { return shared_panel(sc, r, scale.pd, scale.n, rho, T); });
            }
        }
    }
}

// The grade scenario matrix's obligor counts are the ones listed in studies/README.md (D-151), summing to
// 10,000, with the listed portfolio PDs (to two decimals, in percent); the master scale runs 0.03% to 19.7%.
VCAL_TEST(grade_matrix_counts_are_the_studies_index_counts) {
    namespace gm = vcal::recovery::grades;
    const std::vector<std::vector<std::int64_t>> listed = {
        {696, 866, 1012, 1112, 1147, 1112, 1012, 866, 696, 525, 373, 248, 155, 91, 50, 26, 13},
        {224, 337, 475, 629, 783, 916, 1006, 1037, 1006, 916, 783, 629, 475, 337, 224, 140, 83},
        {46, 84, 143, 229, 343, 484, 641, 797, 932, 1024, 1056, 1024, 932, 797, 641, 484, 343}};
    const double portfolio[3] = {0.61, 1.54, 3.42};
    for (std::size_t si = 0; si < 3; ++si) {
        const auto n = gm::obligors(gm::kScales[si]);
        VCAL_CHECK(n == listed[si]);
        std::int64_t total = 0;
        double defaults = 0.0;
        for (int k = 1; k <= gm::kGrades; ++k) {
            total += n[static_cast<std::size_t>(k - 1)];
            defaults += static_cast<double>(n[static_cast<std::size_t>(k - 1)]) * gm::master_pd(k);
        }
        VCAL_CHECK_EQ(total, gm::kObligors);
        VCAL_CHECK(std::fabs(100.0 * defaults / static_cast<double>(total) - portfolio[si]) < 0.005);
    }
    VCAL_CHECK(bits_equal(gm::master_pd(1), 0.0003));
    VCAL_CHECK(std::fabs(gm::master_pd(17) - 0.197) < 0.0005);
    for (int id = 0; id < gm::kScenarios; ++id) {
        const auto sc = gm::scenario(id);
        VCAL_CHECK_EQ(sc.periods, gm::kPeriods[id % 3]);
        VCAL_CHECK(sc.scale == &gm::kScales[id / 3]);
    }
}

// Replicate 0 of the six grade scenarios with T in {20, 40}: the panel is the same twice, and its fit
// converges, is finite and not flat. The CI's look at the matrix; recovery on it is the slow test below.
VCAL_TEST(grade_matrix_replicate_zero_fits) {
    namespace gm = vcal::recovery::grades;
    const auto pd = master_scale();
    for (int id = 0; id < gm::kScenarios; ++id) {
        if (id % 3 == 2) continue;  // T = 100: in the slow test
        const auto p = grade_panel(id, 0);
        const auto again = grade_panel(id, 0);
        bool same = p.y.size() == again.y.size();
        for (std::size_t i = 0; same && i < p.y.size(); ++i) same = p.y[i].d == again.y[i].d;
        VCAL_CHECK(same);
        const auto f = fit(p);
        VCAL_CHECK_EQ(f.flags & (mg::kFlagNotConverged | vcal::engine::kFlagNumeric | vcal::engine::kFlagFlatSurface),
                      0u);
        VCAL_CHECK(std::isfinite(f.loglik) && std::isfinite(f.rho_lo) && std::isfinite(f.rho_hi));
        int edge = 0;
        for (const auto g : f.grade_flags) edge += (g & (mg::kGradeAtLower | mg::kGradeAtUpper)) != 0 ? 1 : 0;
        vcal::test::note(p.name + ": rho " + describe(f.rho) + " [" + describe(f.rho_lo) + ", " + describe(f.rho_hi) +
                         "] against " + describe(gm::scenario(id).scale->rho) + ", " + std::to_string(edge) +
                         " grades on the box, flags " + std::to_string(f.flags));
    }
}

// Recovery on the grade scenario matrix (slow; D-414), with D-413's assertions, R = 100 unless
// VCAL_MG_RECOVERY_R sets it. VCAL_MG_GRADE_SCENARIO runs one scenario (0 .. 8). About 70 CPU-minutes at
// R = 100, two thirds of it the three T = 100 scenarios.
VCAL_TEST(multigrade_fit_recovers_on_the_grade_matrix) {
    namespace gm = vcal::recovery::grades;
    const int R = recovery_replicates();
    const std::string only = vcal::test::env_or_empty("VCAL_MG_GRADE_SCENARIO");
    const auto pd = master_scale();
    for (int id = 0; id < gm::kScenarios; ++id) {
        if (!only.empty() && std::atoi(only.c_str()) != id) continue;
        const auto sc = gm::scenario(id);
        recovery_cell(std::string("grade scenario ") + std::to_string(id) + " (" + sc.scale->name + ", T " +
                          std::to_string(sc.periods) + ", rho " + describe(sc.scale->rho) + ")",
                      pd, gm::obligors(*sc.scale), sc.scale->rho, sc.periods, R,
                      [&](std::uint32_t r) { return grade_panel(id, r); });
    }
}
