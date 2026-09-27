// SPDX-License-Identifier: Apache-2.0
//
// M1.7 cross-reference: the engine (core/ + engine/) against the independent reference (ref/),
// and the engine against its own recorded estimates across platforms. Three comparisons,
// three separately named tolerances (D-107), never interchangeable:
//
//   1. per-period l_t, core (the parity rule, D-118) vs ref, on a fixed grid of (PD, rho, n, d)
//      points, all periods, in term-scaled eps.
//   2. estimates, core vs ref, on a fixed set of DGP panels, in units of core's SE. Both
//      engines search the identical box (D-115). Fits flagged edge/flat/rejected are excluded
//      from the SE comparison: both engines must instead agree on the bound.
//   3. estimates, core on this platform vs core's checked-in regression values, 1e-12 relative
//      (D-102). This is a determinism check (libm last-digit differences only, D-099), not a
//      correctness check; correctness is (1) and (2).
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "dgp/dgp.hpp"
#include "engine/calibrate.hpp"
#include "engine/profile.hpp"
#include "ref/ref.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace e = vcal::engine;
namespace tol = vcal::tol;
using vcal::test::describe;
using Objective = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;



#define VCAL_TEST_STR2(x) #x
#define VCAL_TEST_STR(x) VCAL_TEST_STR2(x)
#if defined(_MSC_VER)
#define VCAL_TEST_COMPILER "MSVC " VCAL_TEST_STR(_MSC_VER)
#elif defined(__clang__)
#define VCAL_TEST_COMPILER "Clang " __clang_version__
#elif defined(__GNUC__)
#define VCAL_TEST_COMPILER "GCC " __VERSION__
#else
#define VCAL_TEST_COMPILER "unknown compiler"
#endif

double logit(double v) { return std::log(v) - std::log1p(-v); }

struct Worst {
    double value = 0.0;
    std::string where;
    int count = 0;
    void add(double v, const std::string& at) {
        ++count;
        if (!(v <= value)) {
            value = v;
            where = at;
        }
    }
};

// --- the estimation setting shared by core and ref (identical box, D-115) ----------------------

constexpr double kPdLo = 1e-4, kPdHi = 0.2, kRhoLo = 1e-3, kRhoHi = vcal::kDefaultRhoUpper;

vcal::Grid<2> estimation_grid() {
    return {{{kPdLo, kPdHi, 61, vcal::AxisScale::Logit}, {kRhoLo, kRhoHi, 41, vcal::AxisScale::Logit}}};
}

// The fixed panel set (D-115): DGP seed, scenario list and replicate numbers are part of the test.
constexpr std::uint64_t kPanelSeed = 0x4D31375852454631ull;  // arbitrary fixed seed
struct PanelCase {
    std::uint32_t scenario;
    std::uint32_t replicate;
    double pd;
    double rho;
    std::int64_t periods;
    std::int64_t obligors;
};
const PanelCase kPanels[] = {
    {0, 0, 0.005, 0.05, 20, 1000}, {1, 0, 0.01, 0.12, 20, 1000}, {2, 0, 0.02, 0.20, 20, 1000},
    {3, 0, 0.001, 0.12, 20, 1000}, {4, 0, 0.05, 0.08, 20, 1000}, {5, 0, 0.01, 0.30, 20, 1000},
};

struct Panel {
    std::vector<Objective::Obs> core;
    std::vector<vcalref::Period> ref;
};

Panel simulate(const PanelCase& c) {
    std::vector<std::int64_t> n(static_cast<std::size_t>(c.periods), c.obligors);
    std::vector<std::int64_t> d(n.size());
    const vcal::dgp::PanelSpec spec{kPanelSeed, c.scenario, c.replicate, c.pd, c.rho, n.data(), c.periods};
    if (vcal::dgp::simulate_panel(spec, d.data(), nullptr) != vcal::dgp::Status::Ok) std::abort();
    Panel p;
    for (std::size_t t = 0; t < n.size(); ++t) {
        p.core.push_back({n[t], d[t]});
        p.ref.push_back({n[t], d[t]});
    }
    return p;
}

struct CoreResult {
    e::Estimate2 est;
    e::ProfileIntervals2 prof;
};

CoreResult core_fit_with_profile(const Panel& p) {
    const auto primary = vcal::quadrature::parity_rule();
    const auto check = vcal::quadrature::parity_rule(true);
    std::vector<double> L;
    CoreResult r{};
    const auto T = static_cast<std::int64_t>(p.core.size());
    const auto st = e::calibrate(vcal::backends::CpuBackend{}, Objective{}, primary, check, p.core.data(), T,
                                 estimation_grid(), L, r.est);
    if (st != e::Status::Ok) std::abort();
    r.prof = e::profile_intervals(Objective{}, primary, p.core.data(), T, estimation_grid(), L, r.est);
    return r;
}

e::Estimate2 core_fit(const Panel& p) { return core_fit_with_profile(p).est; }

std::string panel_label(const PanelCase& c) {
    return "scenario " + std::to_string(c.scenario) + " (PD " + describe(c.pd) + ", rho " + describe(c.rho) + ")";
}

}  // namespace

// --- 1. per-period log-likelihood --------------------------------------------------------------

// l_t = log C(n, d) + log I can be a small difference of two large terms (n = 1e6, d = 350000:
// about +-6.5e5 each, sum -15), so the error scale is the size of the terms, not of the result:
// eps * max(1, |log C| + |log I|).
double term_scaled_eps(double core, double ref, double log_choose) {
    const double scale = std::fmax(1.0, std::fabs(log_choose) + std::fabs(ref - log_choose));
    return std::fabs(core - ref) / (DBL_EPSILON * scale);
}

// Every period, core (the parity rule, D-118) against ref, all rho, in term-scaled eps: adaptive
// GH for 0 < d < n and composite Gauss-Legendre for d in {0, n}.
VCAL_TEST(xref_period_loglik_core_vs_ref) {
    const double pds[] = {1e-6, 1e-5, 1e-4, 1e-3, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5};
    const double rhos[] = {1e-4, 1e-3, 0.01, 0.03, 0.06, 0.12, 0.18, 0.24, 0.35, 0.5, 0.7, 0.9};
    const std::int64_t ns[] = {1, 2, 10, 50, 200, 1000, 10000, 100000, 1000000};
    const auto parity = vcal::quadrature::parity_rule();
    const auto doubled = vcal::quadrature::parity_rule(true);
    Worst interior, one_sided_low, one_sided_high, worst_abs, worst_check;
    int flagged = 0;
    for (const double pd : pds) {
        for (const double rho : rhos) {
            for (const std::int64_t n : ns) {
                const auto expected = static_cast<std::int64_t>(std::floor(static_cast<double>(n) * pd));
                std::vector<std::int64_t> ds = {0, n, std::min(n, std::max<std::int64_t>(1, expected)),
                                                std::min(n, 3 * expected + 1)};
                std::sort(ds.begin(), ds.end());
                ds.erase(std::unique(ds.begin(), ds.end()), ds.end());
                for (const std::int64_t d : ds) {
                    const double v[2] = {pd, rho};
                    const double core = Objective{}.log_contrib({n, d}, Objective::theta(v), parity);
                    const double ref = vcalref::period_loglik(pd, rho, {n, d});
                    const std::string at = "PD " + describe(pd) + " rho " + describe(rho) + " n " + std::to_string(n) +
                                           " d " + std::to_string(d);
                    const double e = term_scaled_eps(core, ref, vcalref::lchoose(n, d));
                    if (d > 0 && d < n) {
                        interior.add(e, at);
                    } else {
                        (rho <= tol::TOL_AGH_REGIME_B_MAX_RHO ? one_sided_low : one_sided_high).add(e, at);
                    }
                    worst_abs.add(std::fabs(core - ref), at);
                    // The per-run check compares the rule with its doubled variant. Where core is at
                    // rounding level, that difference is rounding too; record how large it gets, since
                    // for n ~ 1e6 the terms of l_t are ~1e6 and rounding alone is ~1e-10.
                    const double check = std::fabs(core - Objective{}.log_contrib({n, d}, Objective::theta(v), doubled));
                    worst_check.add(check, at);
                    const double threshold =
                        std::fmax(e::kQuadratureCheckFlagAbs, e::kQuadratureCheckRoundingUlps * DBL_EPSILON *
                                                                  Objective::rounding_scale({n, d}, core));
                    if (check > threshold) ++flagged;
                }
            }
        }
    }
    const auto report = [](const char* name, const Worst& w) {
        vcal::test::note(std::string(name) + ": " + std::to_string(w.count) + " points, worst " + describe(w.value) +
                         " eps (term-scaled) at " + w.where);
    };
    report("0 < d < n", interior);
    report("d in {0, n}, rho <= 0.5", one_sided_low);
    report("d in {0, n}, rho > 0.5", one_sided_high);
    vcal::test::note("worst absolute error vs ref " + describe(worst_abs.value) + " at " + worst_abs.where);
    vcal::test::note("worst |rule - doubled rule| " + describe(worst_check.value) + " at " + worst_check.where + "; " +
                     std::to_string(flagged) + " periods the engine would flag (D-120 threshold)");
    VCAL_CHECK_EQ(flagged, 0);  // every period is at rounding level, so nothing may be flagged
    VCAL_CHECK(interior.value <= tol::TOL_XREF_PERIOD_LL_EPS);
    VCAL_CHECK(one_sided_low.value <= tol::TOL_XREF_PERIOD_LL_EPS);
    VCAL_CHECK(one_sided_high.value <= tol::TOL_XREF_PERIOD_LL_EPS);
}

// --- 2. estimates, core vs ref (Release only: ref takes seconds per fit) -----------------------

VCAL_TEST(xref_estimates_core_vs_ref) {
    Worst pd_se, rho_se, se_ratio, se_near_bound;
    for (const auto& c : kPanels) {
        const Panel p = simulate(c);
        const auto core = core_fit(p);
        const auto ref = vcalref::fit(p.ref, kPdLo, kPdHi, kRhoLo, kRhoHi);
        const std::string at = panel_label(c);
        vcal::test::note(at + ": core PD " + describe(core.value[0]) + " rho " + describe(core.value[1]) + " flags " +
                         std::to_string(core.flags) + " | ref PD " + describe(ref.pd) + " rho " + describe(ref.rho) +
                         (ref.on_boundary ? " (boundary)" : ""));
        const std::uint32_t curvature_flags = e::kFlagGridEdge | e::kFlagFlatSurface | e::kFlagRefinementRejected;
        if (core.flags & curvature_flags) {
            // No meaningful SE: both engines must put the estimate on the same bound, or both flag it.
            const bool core_edge = (core.flags & e::kFlagGridEdge) != 0;
            VCAL_CHECK(core_edge == ref.on_boundary);
            continue;
        }
        VCAL_CHECK(!ref.on_boundary);
        pd_se.add(std::fabs(core.value[0] - ref.pd) / core.se[0], at);
        rho_se.add(std::fabs(core.value[1] - ref.rho) / core.se[1], at);
        const double se_diff = std::fmax(std::fabs(core.se[0] / ref.se_pd - 1), std::fabs(core.se[1] / ref.se_rho - 1));
        // Near a bound Wald SEs are unreliable in any engine (D-119): reported, not compared.
        ((core.flags & e::kFlagNearBound) ? se_near_bound : se_ratio).add(se_diff, at);
    }
    vcal::test::note("PD: worst |core - ref| = " + describe(pd_se.value) + " SE over " + std::to_string(pd_se.count) +
                     " panels, at " + pd_se.where);
    vcal::test::note("rho: worst |core - ref| = " + describe(rho_se.value) + " SE, at " + rho_se.where);
    vcal::test::note("SEs: worst relative difference " + describe(se_ratio.value) + " over " +
                     std::to_string(se_ratio.count) + " panels, at " + se_ratio.where);
    vcal::test::note("SEs near a bound (reported, not compared): worst relative difference " +
                     describe(se_near_bound.value) + " over " + std::to_string(se_near_bound.count) + " panels, at " +
                     se_near_bound.where);
    VCAL_CHECK(pd_se.count >= 4);
    VCAL_CHECK(pd_se.value <= tol::TOL_XREF_ESTIMATE_SE);
    VCAL_CHECK(rho_se.value <= tol::TOL_XREF_ESTIMATE_SE);
    VCAL_CHECK(se_ratio.value <= tol::TOL_XREF_SE_REL);
}

// Both engines on a panel with no defaults: PD's MLE is below the box, so both must report the
// lower bound (core: edge flag and value exactly at the bound; ref: boundary).
VCAL_TEST(xref_estimates_agree_on_the_bound) {
    Panel p;
    for (int t = 0; t < 10; ++t) {
        p.core.push_back({1000, 0});
        p.ref.push_back({1000, 0});
    }
    const auto core = core_fit(p);
    const auto ref = vcalref::fit(p.ref, kPdLo, kPdHi, kRhoLo, kRhoHi);
    VCAL_CHECK(core.flags & e::kFlagGridEdge);
    VCAL_CHECK_EQ(core.value[0], kPdLo);
    VCAL_CHECK(ref.on_boundary);
    VCAL_CHECK_REL(ref.pd, kPdLo, tol::TOL_XREF_REF_BOUND_REL);
}

// --- 2b. profile-likelihood intervals, core vs ref (M2, D-128) -----------------------------------

// Three panels (typical, near the lower rho bound, and near the upper one): both engines solve
// P_a(u) = l_max - 1.92 by unrelated methods (core: Brent on a grid-bracketed profile; ref:
// golden-section over the whole range and bisection). Endpoints compared in logit units;
// truncation must agree.
VCAL_TEST(xref_profile_intervals_core_vs_ref) {
    Worst worst;
    for (const std::uint32_t s : {1u, 3u, 5u}) {
        const auto& c = kPanels[s];
        const Panel p = simulate(c);
        const auto core = core_fit_with_profile(p);
        const auto rf = vcalref::fit(p.ref, kPdLo, kPdHi, kRhoLo, kRhoHi);
        for (int a = 0; a < 2; ++a) {
            const auto ri = vcalref::profile_interval(p.ref, rf, a, kPdLo, kPdHi, kRhoLo, kRhoHi, e::kProfileThreshold95);
            const std::string at = panel_label(c) + (a == 0 ? ", PD" : ", rho");
            vcal::test::note(at + ": core [" + describe(core.prof.lo[a]) + ", " + describe(core.prof.hi[a]) +
                             "] flags " + std::to_string(core.prof.flags[a]) + " | ref [" + describe(ri.lo) + ", " +
                             describe(ri.hi) + "]" + (ri.lo_at_bound ? " lower at bound" : "") +
                             (ri.hi_at_bound ? " upper at bound" : ""));
            VCAL_CHECK_EQ((core.prof.flags[a] & e::kIntervalLowerTruncated) != 0, ri.lo_at_bound);
            VCAL_CHECK_EQ((core.prof.flags[a] & e::kIntervalUpperTruncated) != 0, ri.hi_at_bound);
            worst.add(std::fabs(logit(core.prof.lo[a]) - logit(ri.lo)), at + " lower");
            worst.add(std::fabs(logit(core.prof.hi[a]) - logit(ri.hi)), at + " upper");
        }
    }
    vcal::test::note("worst |core - ref| endpoint, logit units: " + describe(worst.value) + " at " + worst.where);
    VCAL_CHECK(worst.value <= tol::TOL_XREF_PROFILE_ENDPOINT_U);
}

// --- 3. core across platforms: regression values at 1e-12 relative ------------------------------

// Set VCAL_WRITE_REGRESSION_GOLDENS=1 to rewrite tests/golden/xref/core_estimates.csv from this
// platform (then commit it with the platform recorded in its header).
VCAL_TEST(xref_core_estimates_match_regression_values) {
    const std::string path = vcal::test::golden_path("xref/core_estimates.csv");
    if (vcal::test::env_or_empty("VCAL_WRITE_REGRESSION_GOLDENS") == "1") {
        std::ofstream out(path, std::ios::binary);
        out << "# Regression values: core's own estimates on the fixed DGP panels (D-107 comparison 3).\n"
               "# NOT a correctness reference (see xref_estimates_core_vs_ref). Written by this test with\n"
               "# VCAL_WRITE_REGRESSION_GOLDENS=1; compared on every platform at TOL_XREF_CROSS_PLATFORM_REL.\n"
               "# Written by: " VCAL_TEST_COMPILER "\n"
               "# Profile-likelihood intervals (M2) are compared at TOL_PROFILE_CROSS_PLATFORM_REL.\n"
               "scenario,pd_hex,rho_hex,se_pd_hex,se_rho_hex,loglik_hex,flags,pd_lo_hex,pd_hi_hex,rho_lo_hex,"
               "rho_hi_hex,profile_flags\n";
        for (const auto& c : kPanels) {
            const auto r = core_fit_with_profile(simulate(c));
            const auto& est = r.est;
            out << c.scenario << ',' << vcal::test::to_hex(est.value[0]) << ',' << vcal::test::to_hex(est.value[1]) << ','
                << vcal::test::to_hex(est.se[0]) << ',' << vcal::test::to_hex(est.se[1]) << ','
                << vcal::test::to_hex(est.loglik) << ',' << est.flags << ',' << vcal::test::to_hex(r.prof.lo[0]) << ','
                << vcal::test::to_hex(r.prof.hi[0]) << ',' << vcal::test::to_hex(r.prof.lo[1]) << ','
                << vcal::test::to_hex(r.prof.hi[1]) << ',' << (r.prof.flags[0] | (r.prof.flags[1] << 4)) << '\n';
        }
        vcal::test::note("wrote " + path);
        return;
    }
    const auto t = vcal::test::read_golden_csv("xref/core_estimates.csv");
    double worst = 0.0;
    double worst_se = 0.0;
    double worst_profile = 0.0;
    std::string worst_at;
    std::string worst_se_at;
    std::size_t i = 0;
    for (const auto& c : kPanels) {
        VCAL_REQUIRE(i < t.rows.size());
        const auto& r = t.rows[i++];
        const auto res = core_fit_with_profile(simulate(c));
        const auto& est = res.est;
        VCAL_CHECK_EQ(static_cast<std::int64_t>(est.flags), vcal::test::parse_int(r[t.column("flags")]));
        VCAL_CHECK_EQ(static_cast<std::int64_t>(res.prof.flags[0] | (res.prof.flags[1] << 4)),
                      vcal::test::parse_int(r[t.column("profile_flags")]));
        const double ends[4] = {res.prof.lo[0], res.prof.hi[0], res.prof.lo[1], res.prof.hi[1]};
        const char* end_cols[4] = {"pd_lo_hex", "pd_hi_hex", "rho_lo_hex", "rho_hi_hex"};
        for (int k = 0; k < 4; ++k) {
            const double want = vcal::test::parse_double(r[t.column(end_cols[k])]);
            worst_profile = std::fmax(worst_profile, std::fabs(ends[k] / want - 1.0));
        }
        const double got[5] = {est.value[0], est.value[1], est.se[0], est.se[1], est.loglik};
        const char* cols[5] = {"pd_hex", "rho_hex", "se_pd_hex", "se_rho_hex", "loglik_hex"};
        for (int k = 0; k < 5; ++k) {
            const double want = vcal::test::parse_double(r[t.column(cols[k])]);
            if (std::isnan(want) && std::isnan(got[k])) continue;
            const double rel = std::fabs(got[k] / want - 1.0);
            const bool is_se = k == 2 || k == 3;
            double& w = is_se ? worst_se : worst;
            std::string& w_at = is_se ? worst_se_at : worst_at;
            if (rel > w) {
                w = rel;
                w_at = std::string(cols[k]) + " of scenario " + std::to_string(c.scenario);
            }
        }
    }
    vcal::test::note("core vs its regression values: estimates and loglik worst relative difference " +
                     describe(worst) + (worst_at.empty() ? std::string() : " (" + worst_at + ")") + "; SEs " +
                     describe(worst_se) + (worst_se_at.empty() ? std::string() : " (" + worst_se_at + ")") +
                     "; profile endpoints " + describe(worst_profile));
    VCAL_CHECK(worst <= tol::TOL_XREF_CROSS_PLATFORM_REL);
    VCAL_CHECK(worst_se <= tol::TOL_XREF_CROSS_PLATFORM_SE_REL);
    VCAL_CHECK(worst_profile <= tol::TOL_PROFILE_CROSS_PLATFORM_REL);
}

// --- Q16 condition: the parity surface is smooth for zero-default panels (D-118) ---------------

// Refinement and SEs fit a quadratic through the surface, so integration error must not jump from
// grid point to grid point. On zero-default panels (all periods d = 0, and one mixing a few
// defaults), core's second differences along each axis must match ref's; a kink would show up as
// an outlier.
VCAL_TEST(xref_surface_is_smooth_for_zero_default_panels) {
    const vcal::Grid<2> g{{{1e-4, 0.05, 21, vcal::AxisScale::Logit}, {0.01, kRhoHi, 21, vcal::AxisScale::Logit}}};
    const auto parity = vcal::quadrature::parity_rule();
    const std::vector<std::vector<vcalref::Period>> panels = {
        std::vector<vcalref::Period>(8, vcalref::Period{1000, 0}),
        {{200, 0}, {5000, 0}, {1000, 3}, {1000, 0}, {100000, 0}, {50, 0}},
    };
    double worst = 0.0;
    double scale = 0.0;
    for (const auto& panel : panels) {
        std::vector<double> core(static_cast<std::size_t>(g.size()));
        std::vector<double> ref(core.size());
        for (std::int64_t k = 0; k < g.size(); ++k) {
            double v[2];
            g.values(k, v);
            double c = 0.0;
            for (const auto& y : panel) c += Objective{}.log_contrib({y.n, y.d}, Objective::theta(v), parity);
            core[static_cast<std::size_t>(k)] = c;
            ref[static_cast<std::size_t>(k)] = vcalref::loglik(v[0], v[1], panel);
        }
        const auto at = [&](const std::vector<double>& s, std::int32_t i, std::int32_t j) {
            const std::int32_t idx[2] = {i, j};
            return s[static_cast<std::size_t>(g.flatten(idx))];
        };
        for (std::int32_t i = 1; i + 1 < g.axis[0].n; ++i) {
            for (std::int32_t j = 1; j + 1 < g.axis[1].n; ++j) {
                const double d2c[2] = {at(core, i + 1, j) - 2 * at(core, i, j) + at(core, i - 1, j),
                                       at(core, i, j + 1) - 2 * at(core, i, j) + at(core, i, j - 1)};
                const double d2r[2] = {at(ref, i + 1, j) - 2 * at(ref, i, j) + at(ref, i - 1, j),
                                       at(ref, i, j + 1) - 2 * at(ref, i, j) + at(ref, i, j - 1)};
                for (int a = 0; a < 2; ++a) {
                    worst = std::fmax(worst, std::fabs(d2c[a] - d2r[a]));
                    scale = std::fmax(scale, std::fabs(d2r[a]));
                }
            }
        }
    }
    vcal::test::note("zero-default surfaces: worst |second difference, core - ref| = " + describe(worst) +
                     " (largest second difference " + describe(scale) + ", ratio " + describe(worst / scale) + ")");
    VCAL_CHECK(worst / scale <= tol::TOL_XREF_SURFACE_SECOND_DIFF_REL);
}

// --- validation goldens for the scipy replication script (M1.9, D-126) ----------------------------

// The scipy script (validation/scipy/binomial_mixture_mle.py) is compared against core on the
// M1.7 panels: their observations, core's log-likelihood surface on every 5th point of the
// estimation grid (13 x 9 cells), and core's estimates (xref/core_estimates.csv). This test keeps
// the first two current: set VCAL_WRITE_VALIDATION_GOLDENS=1 to rewrite them; otherwise core must
// reproduce them at TOL_XREF_CROSS_PLATFORM_REL (the same comparison as D-107 no. 3, on the same
// panels).
VCAL_TEST(validation_goldens_are_current) {
    constexpr std::int32_t kStride = 5;
    const auto g = estimation_grid();
    const auto parity = vcal::quadrature::parity_rule();
    struct Cell {
        std::uint32_t scenario;
        std::int32_t i0, i1;
        double pd, rho, loglik;
    };
    std::vector<Cell> cells;
    std::string panels = "scenario,period,n,d\n";
    for (const auto& c : kPanels) {
        const Panel p = simulate(c);
        const auto T = static_cast<std::int64_t>(p.core.size());
        for (std::int64_t t = 0; t < T; ++t) {
            const auto& y = p.core[static_cast<std::size_t>(t)];
            panels += std::to_string(c.scenario) + "," + std::to_string(t) + "," + std::to_string(y.n) + "," +
                      std::to_string(y.d) + "\n";
        }
        std::vector<double> L(static_cast<std::size_t>(T * g.size()));
        e::evaluate_surface(vcal::backends::CpuBackend{}, Objective{}, parity, p.core.data(), T, g, L.data());
        const std::vector<double> ones(static_cast<std::size_t>(T), 1.0);
        for (std::int32_t i0 = 0; i0 < g.axis[0].n; i0 += kStride) {
            for (std::int32_t i1 = 0; i1 < g.axis[1].n; i1 += kStride) {
                const std::int32_t idx[2] = {i0, i1};
                const std::int64_t k = g.flatten(idx);
                double v[2];
                g.values(k, v);
                cells.push_back({c.scenario, i0, i1, v[0], v[1], e::weighted_sum(L.data(), T, g.size(), ones.data(), k)});
            }
        }
    }
    const std::string panels_path = vcal::test::golden_path("validation/panels.csv");
    const std::string surface_path = vcal::test::golden_path("validation/surface.csv");
    if (vcal::test::env_or_empty("VCAL_WRITE_VALIDATION_GOLDENS") == "1") {
        std::filesystem::create_directories(std::filesystem::path(panels_path).parent_path());
        std::ofstream pf(panels_path, std::ios::binary);
        pf << "# The six M1.7 DGP panels (D-115; seed 0x4D31375852454631, replicate 0), as observed data for\n"
              "# validation/scipy/binomial_mixture_mle.py. Written by unit_xref validation_goldens_are_current.\n"
           << panels;
        std::ofstream sf(surface_path, std::ios::binary);
        sf << "# Core's panel log-likelihood (the parity rule, D-118) on every 5th point of the estimation grid:\n"
              "# PD logit [1e-4, 0.2] x 61, rho logit [1e-3, 0.5] x 41. Written by " VCAL_TEST_COMPILER ".\n"
              "scenario,i_pd,i_rho,pd_hex,rho_hex,loglik_hex,pd,rho,loglik\n";
        for (const auto& x : cells) {
            sf << x.scenario << ',' << x.i0 << ',' << x.i1 << ',' << vcal::test::to_hex(x.pd) << ','
               << vcal::test::to_hex(x.rho) << ',' << vcal::test::to_hex(x.loglik) << ',' << describe(x.pd) << ','
               << describe(x.rho) << ',' << describe(x.loglik) << '\n';
        }
        pf.close();
        sf.close();
        VCAL_REQUIRE(pf && sf);
        vcal::test::note("wrote " + panels_path + " and " + surface_path);
        return;
    }
    const auto pt = vcal::test::read_golden_csv("validation/panels.csv");
    std::string recorded = "scenario,period,n,d\n";
    for (const auto& r : pt.rows) recorded += r[0] + "," + r[1] + "," + r[2] + "," + r[3] + "\n";
    VCAL_CHECK(recorded == panels);  // the DGP is bitwise reproducible (D-057)
    const auto st = vcal::test::read_golden_csv("validation/surface.csv");
    VCAL_REQUIRE(st.rows.size() == cells.size());
    double worst = 0.0;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const double want = vcal::test::parse_double(st.rows[i][st.column("loglik_hex")]);
        worst = std::fmax(worst, std::fabs(cells[i].loglik / want - 1.0));
    }
    vcal::test::note(std::to_string(cells.size()) + " surface cells: worst relative difference from the golden " +
                     describe(worst));
    VCAL_CHECK(worst <= tol::TOL_XREF_CROSS_PLATFORM_REL);
}
