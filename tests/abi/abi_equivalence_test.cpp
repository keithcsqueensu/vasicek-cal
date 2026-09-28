// SPDX-License-Identifier: Apache-2.0
//
// The C ABI against the engine (M2c, D-144): every vcal_* result must be bitwise what the engine
// gives when called directly with the same inputs. Together with the recovery replay (engine
// against goldens) this ties the ABI to the pinned results. The C test (abi_c_test.c) covers
// the conventions a C caller sees; this one covers the numbers.
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/quadrature/parity.hpp"
#include "dgp/dgp.hpp"
#include "core/objectives/vasicek_rate.hpp"
#include "engine/calibrate.hpp"
#include "engine/moments.hpp"
#include "engine/posterior.hpp"
#include "engine/profile.hpp"
#include "resample/bootstrap.hpp"
#include "resample/weights.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/recovery/recovery.hpp"
#include "vcal/vcal.h"

namespace {

namespace e = vcal::engine;
namespace rs = vcal::resample;
using Objective = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
using Obs = Objective::Obs;

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

bool bits_equal(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

template <class T>
T sized() {
    T v{};
    v.struct_size = static_cast<std::uint32_t>(sizeof(T));
    return v;
}

struct Context {
    explicit Context(std::int32_t threads = 0) {
        auto o = sized<vcal_context_options>();
        o.n_threads = threads;
        VCAL_REQUIRE(vcal_context_create(&o, &ctx) == VCAL_OK);
    }
    ~Context() { vcal_context_destroy(ctx); }
    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;
    vcal_context* ctx = nullptr;
};

// A recovery-study panel, simulated through the engine's DGP.
struct Panel {
    Panel(std::uint32_t scenario, std::uint32_t replicate) {
        const auto s = vcal::recovery::scenario(scenario);
        n.assign(static_cast<std::size_t>(s.periods), s.obligors);
        d.resize(n.size());
        const vcal::dgp::PanelSpec spec{vcal::recovery::kSeed, s.id, replicate, s.pd, s.rho, n.data(), s.periods};
        VCAL_REQUIRE(vcal::dgp::simulate_panel(spec, d.data(), nullptr) == vcal::dgp::Status::Ok);
        for (std::size_t t = 0; t < n.size(); ++t) obs.push_back({n[t], d[t]});
        abi = sized<vcal_panel>();
        abi.n_periods = static_cast<std::int64_t>(n.size());
        abi.n_obligors = n.data();
        abi.n_defaults = d.data();
    }
    std::int64_t periods() const { return static_cast<std::int64_t>(n.size()); }
    std::vector<std::int64_t> n, d;
    std::vector<Obs> obs;
    vcal_panel abi{};
};

vcal_grid default_grid() {
    auto g = sized<vcal_grid>();
    VCAL_REQUIRE(vcal_grid_default(&g) == VCAL_OK);
    return g;
}

// Scenarios: 40 (PD 1%, rho 0.12, T = 40, n = 1000), an ordinary fit; 0 (PD 0.1%, rho 0.02,
// T = 20, n = 100), mostly zero-default periods with the estimate at a bound.
constexpr std::uint32_t kScenarios[] = {40, 0};

}  // namespace

VCAL_TEST(abi_default_grid_is_the_recovery_grid) {
    const vcal_grid g = default_grid();
    const vcal::Grid<2> r = vcal::recovery::grid();
    VCAL_CHECK(g.pd_lo == r.axis[0].lo && g.pd_hi == r.axis[0].hi && g.pd_points == r.axis[0].n);
    VCAL_CHECK(g.rho_lo == r.axis[1].lo && g.rho_hi == r.axis[1].hi && g.rho_points == r.axis[1].n);
    VCAL_CHECK(g.pd_scale == VCAL_SCALE_LOGIT && g.rho_scale == VCAL_SCALE_LOGIT);
    for (int a = 0; a < 2; ++a) {
        std::vector<double> v(static_cast<std::size_t>(r.axis[a].n));
        VCAL_REQUIRE(vcal_grid_values(&g, a, v.data(), static_cast<std::int64_t>(v.size()), nullptr) == VCAL_OK);
        for (std::int32_t i = 0; i < r.axis[a].n; ++i) VCAL_CHECK(bits_equal(v[static_cast<std::size_t>(i)], r.axis[a].value_at(i)));
    }
}

VCAL_TEST(abi_calibrate_and_profile_equal_the_engine) {
    const Context c;
    const vcal_grid grid = default_grid();
    const vcal::Grid<2> g = vcal::recovery::grid();
    const auto primary = vcal::quadrature::parity_rule();
    const auto check = vcal::quadrature::parity_rule(true);
    for (const auto scenario : kScenarios) {
        const Panel p(scenario, 0);
        std::vector<double> L;
        e::Estimate2 est{};
        VCAL_REQUIRE(e::calibrate(vcal::backends::CpuBackend{}, Objective{}, primary, check, p.obs.data(), p.periods(), g,
                                  L, est) == e::Status::Ok);
        const e::ProfileIntervals2 prof = e::profile_intervals(Objective{}, primary, p.obs.data(), p.periods(), g, L, est);

        auto a = sized<vcal_estimate>();
        auto ap = sized<vcal_profile_intervals>();
        VCAL_REQUIRE(vcal_calibrate(c.ctx, &p.abi, &grid, &a, &ap) == VCAL_OK);
        VCAL_CHECK_EQ(a.flags, est.flags);
        VCAL_CHECK(bits_equal(a.pd, est.value[0]) && bits_equal(a.rho, est.value[1]));
        VCAL_CHECK(bits_equal(a.se_pd, est.se[0]) && bits_equal(a.se_rho, est.se[1]) && bits_equal(a.corr_pd_rho, est.corr));
        VCAL_CHECK(bits_equal(a.loglik, est.loglik));
        VCAL_CHECK(bits_equal(a.quad_check_max, est.quad_check_max) && bits_equal(a.quad_check_total, est.quad_check_total));
        VCAL_CHECK_EQ(a.quad_check_flagged, est.quad_check_flagged);
        VCAL_CHECK_EQ(a.grid_index, est.grid_index);
        VCAL_CHECK_EQ(a.nan_count, est.nan_count);
        VCAL_CHECK_EQ(ap.pd_flags, prof.flags[0]);
        VCAL_CHECK_EQ(ap.rho_flags, prof.flags[1]);
        VCAL_CHECK(bits_equal(ap.pd_lo, prof.lo[0]) && bits_equal(ap.pd_hi, prof.hi[0]));
        VCAL_CHECK(bits_equal(ap.rho_lo, prof.lo[1]) && bits_equal(ap.rho_hi, prof.hi[1]));
        VCAL_CHECK(bits_equal(ap.loglik_max, prof.loglik_max) && bits_equal(ap.residual_max, prof.residual_max));
        VCAL_CHECK(bits_equal(ap.pd_at_max, prof.max_at[0]) && bits_equal(ap.rho_at_max, prof.max_at[1]));
        VCAL_CHECK_EQ(ap.evaluations, prof.evaluations);
        vcal::test::note("scenario " + std::to_string(scenario) + ": flags " + std::to_string(a.flags) + ", PD " +
                         vcal::test::describe(a.pd) + " [" + vcal::test::describe(ap.pd_lo) + ", " +
                         vcal::test::describe(ap.pd_hi) + "], rho " + vcal::test::describe(a.rho) + " [" +
                         vcal::test::describe(ap.rho_lo) + ", " + vcal::test::describe(ap.rho_hi) + "]");

        // The surface is calibrate's L.
        std::int64_t required = 0;
        VCAL_REQUIRE(vcal_surface(c.ctx, &p.abi, &grid, nullptr, 0, &required) == VCAL_OK);
        std::vector<double> S(static_cast<std::size_t>(required));
        VCAL_REQUIRE(vcal_surface(c.ctx, &p.abi, &grid, S.data(), required, nullptr) == VCAL_OK);
        VCAL_CHECK(bits_equal(S, L));
    }
}

VCAL_TEST(abi_results_do_not_depend_on_the_thread_count) {
    const vcal_grid grid = default_grid();
    const Panel p(40, 1);
    auto spec = sized<vcal_resample_spec>();
    spec.scheme = VCAL_RESAMPLE_IID_BOOTSTRAP;
    spec.seed = 7;
    spec.replicates = 64;
    std::vector<vcal_estimate> est;
    std::vector<std::vector<double>> pd;
    for (const std::int32_t threads : {1, 3, 0}) {
        const Context c(threads);
        auto a = sized<vcal_estimate>();
        VCAL_REQUIRE(vcal_calibrate(c.ctx, &p.abi, &grid, &a, nullptr) == VCAL_OK);
        est.push_back(a);
        std::vector<double> v(64);
        auto reps = sized<vcal_replicates>();
        reps.capacity = 64;
        reps.pd = v.data();
        VCAL_REQUIRE(vcal_resample(c.ctx, &p.abi, &grid, &spec, &reps, nullptr, nullptr) == VCAL_OK);
        pd.push_back(v);
    }
    for (std::size_t i = 1; i < est.size(); ++i) {
        VCAL_CHECK(std::memcmp(&est[i], &est[0], sizeof(vcal_estimate)) == 0);
        VCAL_CHECK(bits_equal(pd[i], pd[0]));
    }
}

VCAL_TEST(abi_resample_equals_the_engine_for_every_scheme) {
    const Context c;
    const vcal_grid grid = default_grid();
    const vcal::Grid<2> g = vcal::recovery::grid();
    const Panel p(40, 0);
    const std::int64_t T = p.periods();
    std::vector<double> L(static_cast<std::size_t>(T * g.size()));
    e::evaluate_surface(vcal::backends::CpuBackend{}, Objective{}, vcal::quadrature::parity_rule(), p.obs.data(), T, g,
                        L.data());

    std::vector<std::int64_t> supplied;  // 3 replicates of 25 draws (an m-out-of-n bootstrap)
    for (std::int64_t i = 0; i < 75; ++i) supplied.push_back((i * 17 + 3) % T);
    std::vector<double> custom(static_cast<std::size_t>(2 * T), 1.0);
    for (std::int64_t t = 0; t < T; ++t) custom[static_cast<std::size_t>(T + t)] = 0.5 + static_cast<double>(t % 3);

    struct Case {
        const char* name;
        vcal_resample_spec spec;
        std::vector<double> W;  // what the engine's builders give
    };
    std::vector<Case> cases;
    auto make = [&](std::int32_t scheme) {
        auto s = sized<vcal_resample_spec>();
        s.scheme = scheme;
        return s;
    };
    {
        auto s = make(VCAL_RESAMPLE_IID_BOOTSTRAP);
        s.seed = 11;
        s.replicates = 99;
        s.level = 0.9;
        const auto idx = rs::bootstrap_indices(11, rs::Scheme::IidBootstrap, 99, T);
        cases.push_back({"iid", s, rs::weights_from_indices(idx.data(), 99, T, T)});
    }
    {
        auto s = make(VCAL_RESAMPLE_BLOCK_BOOTSTRAP);
        s.seed = 12;
        s.replicates = 50;  // block_length 0: ceil(40^(1/3)) = 4
        const auto idx = rs::bootstrap_indices(12, rs::Scheme::BlockBootstrap, 50, T, 4);
        cases.push_back({"block, default length", s, rs::weights_from_indices(idx.data(), 50, T, T)});
        s.block_length = 7;
        const auto idx7 = rs::bootstrap_indices(12, rs::Scheme::BlockBootstrap, 50, T, 7);
        cases.push_back({"block, length 7", s, rs::weights_from_indices(idx7.data(), 50, T, T)});
    }
    cases.push_back({"jackknife", make(VCAL_RESAMPLE_JACKKNIFE), rs::jackknife_weights(T)});
    {
        auto s = make(VCAL_RESAMPLE_WALK_FORWARD);
        s.window = 20;
        s.first_end = 25;
        s.step = 3;
        cases.push_back({"walk-forward", s, rs::walk_forward_weights(T, 20, 25, 3)});
    }
    {
        auto s = make(VCAL_RESAMPLE_INDICES);
        s.replicates = 3;
        s.draws = 25;
        s.indices = supplied.data();
        cases.push_back({"indices", s, rs::weights_from_indices(supplied.data(), 3, 25, T)});
    }
    {
        auto s = make(VCAL_RESAMPLE_WEIGHTS);
        s.replicates = 2;
        s.weights = custom.data();
        cases.push_back({"weights", s, custom});
    }

    for (const auto& k : cases) {
        const std::int64_t B = static_cast<std::int64_t>(k.W.size()) / T;
        std::int64_t required = 0;
        VCAL_REQUIRE(vcal_resample_weights(c.ctx, &k.spec, T, nullptr, 0, &required) == VCAL_OK);
        VCAL_CHECK_EQ(required, B * T);
        std::vector<double> W(static_cast<std::size_t>(required));
        VCAL_REQUIRE(vcal_resample_weights(c.ctx, &k.spec, T, W.data(), required, nullptr) == VCAL_OK);
        VCAL_CHECK(bits_equal(W, k.W));

        std::vector<rs::Replicate2> expect(static_cast<std::size_t>(B));
        rs::replicate_estimates(vcal::backends::CpuBackend{}, g, L.data(), T, k.W.data(), B, expect.data());
        const double level = k.spec.level == 0.0 ? 0.95 : k.spec.level;
        const auto i0 = rs::percentile_interval(expect.data(), B, 0, level);
        const auto i1 = rs::percentile_interval(expect.data(), B, 1, level);

        std::vector<double> pd(static_cast<std::size_t>(B)), rho(pd.size()), ll(pd.size());
        std::vector<std::int64_t> gi(pd.size());
        std::vector<std::uint32_t> flags(pd.size());
        auto reps = sized<vcal_replicates>();
        reps.capacity = B;
        reps.pd = pd.data();
        reps.rho = rho.data();
        reps.grid_loglik = ll.data();
        reps.grid_index = gi.data();
        reps.flags = flags.data();
        auto iv = sized<vcal_percentile_intervals>();
        VCAL_REQUIRE(vcal_resample(c.ctx, &p.abi, &grid, &k.spec, &reps, &iv, &required) == VCAL_OK);
        VCAL_CHECK_EQ(required, B);
        bool same = true;
        std::int64_t edge = 0;
        for (std::size_t b = 0; b < expect.size(); ++b) {
            const auto& x = expect[b];
            same = same && bits_equal(pd[b], x.value[0]) && bits_equal(rho[b], x.value[1]) &&
                   bits_equal(ll[b], x.surface_max) && gi[b] == x.grid_index && flags[b] == x.flags;
            edge += (x.flags & e::kFlagGridEdge) ? 1 : 0;
        }
        VCAL_CHECK(same);
        VCAL_CHECK(bits_equal(iv.level, level));
        VCAL_CHECK(bits_equal(iv.pd_lo, i0.lo) && bits_equal(iv.pd_hi, i0.hi));
        VCAL_CHECK(bits_equal(iv.rho_lo, i1.lo) && bits_equal(iv.rho_hi, i1.hi));
        VCAL_CHECK_EQ(iv.replicates, B);
        VCAL_CHECK_EQ(iv.pd_excluded, i0.excluded);
        VCAL_CHECK_EQ(iv.rho_excluded, i1.excluded);
        VCAL_CHECK_EQ(iv.grid_edge, edge);
        vcal::test::note(std::string(k.name) + ": B = " + std::to_string(B) + ", bitwise " + (same ? "equal" : "DIFFERENT"));
    }
}

VCAL_TEST(abi_default_block_length_is_the_integer_cube_root) {
    const Context c;
    // ceil(T^(1/3)) at and around perfect cubes, where a floating-point cube root could round.
    const std::int64_t periods[] = {1, 2, 8, 9, 26, 27, 28, 63, 64, 65, 124, 125, 126, 1000, 1001};
    const std::int64_t expect[] = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 10, 11};
    for (std::size_t i = 0; i < sizeof periods / sizeof periods[0]; ++i) {
        const std::int64_t T = periods[i];
        auto s = sized<vcal_resample_spec>();
        s.scheme = VCAL_RESAMPLE_BLOCK_BOOTSTRAP;
        s.seed = 5;
        s.replicates = 4;
        std::vector<double> dflt(static_cast<std::size_t>(4 * T)), explicit_l(dflt.size());
        VCAL_REQUIRE(vcal_resample_weights(c.ctx, &s, T, dflt.data(), 4 * T, nullptr) == VCAL_OK);
        s.block_length = expect[i];
        VCAL_REQUIRE(vcal_resample_weights(c.ctx, &s, T, explicit_l.data(), 4 * T, nullptr) == VCAL_OK);
        VCAL_CHECK(bits_equal(dflt, explicit_l));
    }
}

VCAL_TEST(abi_dgp_equals_the_engine) {
    const Context c;
    std::vector<std::int64_t> n = {10, 0, 1000, 250, 1, 5000};
    auto spec = sized<vcal_dgp_spec>();
    spec.seed = vcal::recovery::kSeed;
    spec.scenario = 3;
    spec.replicate = 17;
    spec.pd = 0.05;
    spec.rho = 0.24;
    spec.n_periods = static_cast<std::int64_t>(n.size());
    spec.n_obligors = n.data();
    std::vector<std::int64_t> d(n.size()), d_engine(n.size());
    std::vector<double> z(n.size()), z_engine(n.size());
    VCAL_REQUIRE(vcal_dgp_simulate(c.ctx, &spec, d.data(), z.data(), spec.n_periods, nullptr) == VCAL_OK);
    const vcal::dgp::PanelSpec ps{spec.seed, 3, 17, 0.05, 0.24, n.data(), spec.n_periods};
    VCAL_REQUIRE(vcal::dgp::simulate_panel(ps, d_engine.data(), z_engine.data()) == vcal::dgp::Status::Ok);
    VCAL_CHECK(d == d_engine);
    VCAL_CHECK(bits_equal(z, z_engine));
    VCAL_CHECK_EQ(d[1], std::int64_t{0});  // a period with no obligors has no defaults
}

// ---- ABI 0.3: the M3 estimators (D-169) ----------------------------------------------------------

namespace {

namespace ob = vcal::objectives;
using RateRefuse = ob::VasicekRate<vcal::PrecisionF64, ob::ZeroRates::Refuse>;
using RateCensor = ob::VasicekRate<vcal::PrecisionF64, ob::ZeroRates::Censor>;
using RateSubstitute = ob::VasicekRate<vcal::PrecisionF64, ob::ZeroRates::Substitute>;

struct NativeContext {
    NativeContext() {
        auto o = sized<vcal_context_options>();
        o.profile = VCAL_PROFILE_NATIVE;
        VCAL_REQUIRE(vcal_context_create(&o, &ctx) == VCAL_OK);
    }
    ~NativeContext() { vcal_context_destroy(ctx); }
    NativeContext(const NativeContext&) = delete;
    NativeContext& operator=(const NativeContext&) = delete;
    vcal_context* ctx = nullptr;
};

template <class O>
void check_rate_fit(const std::vector<ob::RateObs>& obs, const vcal_estimate& a, const vcal_profile_intervals& ap) {
    const vcal::Grid<2> g = vcal::recovery::grid();
    const auto primary = vcal::quadrature::parity_rule();
    const auto check = vcal::quadrature::parity_rule(true);
    const auto T = static_cast<std::int64_t>(obs.size());
    std::vector<double> L;
    e::Estimate2 est{};
    VCAL_REQUIRE(e::calibrate(vcal::backends::CpuBackend{}, O{}, primary, check, obs.data(), T, g, L, est) == e::Status::Ok);
    const auto prof = e::profile_intervals(O{}, primary, obs.data(), T, g, L, est);
    VCAL_CHECK_EQ(a.flags, est.flags);
    VCAL_CHECK(bits_equal(a.pd, est.value[0]) && bits_equal(a.rho, est.value[1]));
    VCAL_CHECK(bits_equal(a.se_pd, est.se[0]) && bits_equal(a.se_rho, est.se[1]) && bits_equal(a.loglik, est.loglik));
    VCAL_CHECK(bits_equal(ap.pd_lo, prof.lo[0]) && bits_equal(ap.pd_hi, prof.hi[0]));
    VCAL_CHECK(bits_equal(ap.rho_lo, prof.lo[1]) && bits_equal(ap.rho_hi, prof.hi[1]));
    VCAL_CHECK(bits_equal(ap.pd_at_max, prof.max_at[0]) && bits_equal(ap.rho_at_max, prof.max_at[1]));
}

}  // namespace

// The rate MLE on count data and on a rate series, each treatment, equals the engine bit for bit.
VCAL_TEST(abi_rate_mle_equals_the_engine) {
    const vcal_grid grid = default_grid();
    // Count data without zeros (scenario 41: PD 1%, rho 0.12, T = 40, n = 10^4): parity refuses nothing.
    {
        const Context c;
        const Panel p(41, 0);
        std::vector<ob::RateObs> obs;
        for (const auto& y : p.obs) obs.push_back(ob::rate_obs(y.n, y.d));
        VCAL_REQUIRE(ob::zero_rate_periods(obs.data(), static_cast<std::int64_t>(obs.size())).empty());
        auto a = sized<vcal_estimate>();
        auto ap = sized<vcal_profile_intervals>();
        VCAL_REQUIRE(vcal_calibrate_rate(c.ctx, &p.abi, nullptr, &grid, VCAL_ZERO_RATES_REFUSE, &a, &ap) == VCAL_OK);
        check_rate_fit<RateRefuse>(obs, a, ap);
        VCAL_CHECK(a.quad_check_max == 0.0 && a.quad_check_flagged == 0);
    }
    // Count data with zero-default periods (scenario 30: PD 1%, rho 0.02, T = 40, n = 100), native.
    const NativeContext c;
    const Panel p(30, 0);
    std::vector<ob::RateObs> obs;
    for (const auto& y : p.obs) obs.push_back(ob::rate_obs(y.n, y.d));
    VCAL_REQUIRE(!ob::zero_rate_periods(obs.data(), static_cast<std::int64_t>(obs.size())).empty());
    {
        auto a = sized<vcal_estimate>();
        auto ap = sized<vcal_profile_intervals>();
        VCAL_REQUIRE(vcal_calibrate_rate(c.ctx, &p.abi, nullptr, &grid, VCAL_ZERO_RATES_CENSOR, &a, &ap) == VCAL_OK);
        check_rate_fit<RateCensor>(obs, a, ap);
    }
    {
        auto a = sized<vcal_estimate>();
        auto ap = sized<vcal_profile_intervals>();
        VCAL_REQUIRE(vcal_calibrate_rate(c.ctx, &p.abi, nullptr, &grid, VCAL_ZERO_RATES_SUBSTITUTE, &a, &ap) == VCAL_OK);
        check_rate_fit<RateSubstitute>(obs, a, ap);
    }
    {
        auto a = sized<vcal_estimate>();
        auto ap = sized<vcal_profile_intervals>();
        VCAL_REQUIRE(vcal_calibrate_rate(c.ctx, &p.abi, nullptr, &grid, VCAL_ZERO_RATES_DROP, &a, &ap) == VCAL_OK);
        check_rate_fit<RateRefuse>(ob::drop_zero_rate_periods(obs.data(), static_cast<std::int64_t>(obs.size())), a, ap);
    }
    // The same data as a rate series with its detection limits.
    {
        std::vector<double> r, lim;
        for (const auto& y : obs) {
            r.push_back(y.rate);
            lim.push_back(y.detect);
        }
        auto rs = sized<vcal_rate_series>();
        rs.n_periods = static_cast<std::int64_t>(r.size());
        rs.rates = r.data();
        rs.detection_limits = lim.data();
        auto a = sized<vcal_estimate>();
        auto ap = sized<vcal_profile_intervals>();
        VCAL_REQUIRE(vcal_calibrate_rate(c.ctx, nullptr, &rs, &grid, VCAL_ZERO_RATES_CENSOR, &a, &ap) == VCAL_OK);
        check_rate_fit<RateCensor>(obs, a, ap);
    }
}

// Parity refuses zero rates and names the periods; native treatments need a native context; a
// series with zero rates needs its detection limits.
VCAL_TEST(abi_rate_mle_refusals) {
    const vcal_grid grid = default_grid();
    const Panel p(30, 0);
    const Context parity;
    auto a = sized<vcal_estimate>();
    VCAL_CHECK(vcal_calibrate_rate(parity.ctx, &p.abi, nullptr, &grid, VCAL_ZERO_RATES_REFUSE, &a, nullptr) ==
               VCAL_E_INVALID_ARGUMENT);
    char msg[512];
    VCAL_REQUIRE(vcal_last_error(msg, sizeof msg, nullptr) == VCAL_OK);
    VCAL_CHECK(std::strstr(msg, "refused (D-044): periods ") != nullptr);
    VCAL_CHECK(vcal_calibrate_rate(parity.ctx, &p.abi, nullptr, &grid, VCAL_ZERO_RATES_CENSOR, &a, nullptr) ==
               VCAL_E_INVALID_ARGUMENT);
    VCAL_REQUIRE(vcal_last_error(msg, sizeof msg, nullptr) == VCAL_OK);
    VCAL_CHECK(std::strstr(msg, "VCAL_PROFILE_NATIVE") != nullptr);
    const NativeContext native;
    VCAL_CHECK(vcal_calibrate_rate(native.ctx, &p.abi, nullptr, &grid, 7, &a, nullptr) == VCAL_E_INVALID_ARGUMENT);
    VCAL_CHECK(vcal_calibrate_rate(native.ctx, nullptr, nullptr, &grid, VCAL_ZERO_RATES_CENSOR, &a, nullptr) ==
               VCAL_E_INVALID_ARGUMENT);
    const double r[3] = {0.01, 0.0, 0.02};
    auto rs = sized<vcal_rate_series>();
    rs.n_periods = 3;
    rs.rates = r;
    VCAL_CHECK(vcal_calibrate_rate(native.ctx, nullptr, &rs, &grid, VCAL_ZERO_RATES_CENSOR, &a, nullptr) ==
               VCAL_E_INVALID_ARGUMENT);  // a zero rate without detection limits
}

// The method of moments, counts and rates, equals the engine bit for bit.
VCAL_TEST(abi_moments_equal_the_engine) {
    const Context c;
    const vcal_grid grid = default_grid();
    const vcal::Grid<2> g = vcal::recovery::grid();
    const auto integrator = vcal::quadrature::parity_rule();
    for (const auto scenario : {40u, 0u, 30u}) {
        const Panel p(scenario, 0);
        const auto m = e::mom_from_counts(integrator, p.obs.data(), p.periods(), g.axis[0].lo, g.axis[0].hi, g.axis[1].lo,
                                          g.axis[1].hi);
        auto a = sized<vcal_moments_estimate>();
        VCAL_REQUIRE(vcal_calibrate_moments(c.ctx, &p.abi, nullptr, &grid, &a) == VCAL_OK);
        VCAL_CHECK_EQ(a.flags, m.flags);
        VCAL_CHECK(bits_equal(a.pd, m.pd) && bits_equal(a.rho, m.rho) && bits_equal(a.pd2, m.pd2));
        std::vector<double> r;
        for (const auto& y : p.obs) r.push_back(static_cast<double>(y.d) / static_cast<double>(y.n));
        const auto mr = e::mom_from_rates(integrator, r.data(), p.periods(), g.axis[0].lo, g.axis[0].hi, g.axis[1].lo,
                                          g.axis[1].hi);
        auto rs = sized<vcal_rate_series>();
        rs.n_periods = p.periods();
        rs.rates = r.data();
        VCAL_REQUIRE(vcal_calibrate_moments(c.ctx, nullptr, &rs, &grid, &a) == VCAL_OK);
        VCAL_CHECK_EQ(a.flags, mr.flags);
        VCAL_CHECK(bits_equal(a.pd, mr.pd) && bits_equal(a.rho, mr.rho) && bits_equal(a.pd2, mr.pd2));
    }
}

// The grid posterior, both priors, equals the engine bit for bit; the Jeffreys table is reused
// within a context; unequal n with Jeffreys is unsupported.
VCAL_TEST(abi_posterior_equals_the_engine) {
    const Context c;
    const vcal_grid grid = default_grid();
    const vcal::Grid<2> g = vcal::recovery::grid();
    const auto primary = vcal::quadrature::parity_rule();
    const auto check = vcal::quadrature::parity_rule(true);
    const Panel p(37, 0);  // PD 1%, rho 0.12, T = 20, n = 1000
    std::vector<double> L;
    e::Estimate2 est{};
    VCAL_REQUIRE(e::calibrate(vcal::backends::CpuBackend{}, Objective{}, primary, check, p.obs.data(), p.periods(), g, L,
                              est) == e::Status::Ok);
    double se_u[2];
    for (int a = 0; a < 2; ++a) {
        se_u[a] = est.se[a] / vcal::grid::dvalue_dscaled(g.axis[a].scale, vcal::grid::to_scaled(g.axis[a].scale, est.value[a]));
    }
    const auto jt = e::jeffreys_table(vcal::backends::CpuBackend{}, Objective{}, primary, p.n[0], g);
    for (const auto prior : {VCAL_PRIOR_FLAT, VCAL_PRIOR_JEFFREYS}) {
        const auto r = e::grid_posterior(vcal::backends::CpuBackend{}, Objective{}, primary, p.obs.data(), p.periods(), g, L,
                                         prior == VCAL_PRIOR_FLAT ? e::Prior::Flat : e::Prior::Jeffreys, &jt, se_u);
        for (int call = 0; call < 2; ++call) {  // the second Jeffreys call uses the context's cached table
            auto a = sized<vcal_posterior>();
            VCAL_REQUIRE(vcal_calibrate_posterior(c.ctx, &p.abi, &grid, prior, &a) == VCAL_OK);
            VCAL_CHECK_EQ(a.flags, r.flags);
            VCAL_CHECK_EQ(a.refinements, r.refinements);
            VCAL_CHECK(bits_equal(a.pd_et_lo, r.et_lo[0]) && bits_equal(a.pd_et_hi, r.et_hi[0]));
            VCAL_CHECK(bits_equal(a.pd_hpd_lo, r.hpd_lo[0]) && bits_equal(a.pd_hpd_hi, r.hpd_hi[0]));
            VCAL_CHECK(bits_equal(a.rho_et_lo, r.et_lo[1]) && bits_equal(a.rho_et_hi, r.et_hi[1]));
            VCAL_CHECK(bits_equal(a.rho_hpd_lo, r.hpd_lo[1]) && bits_equal(a.rho_hpd_hi, r.hpd_hi[1]));
            VCAL_CHECK(bits_equal(a.pd_mean_logit, r.mean[0]) && bits_equal(a.rho_sd_logit, r.sd[1]));
        }
    }
    std::vector<std::int64_t> n = p.n, d = p.d;
    n[0] += 1;
    auto uneq = sized<vcal_panel>();
    uneq.n_periods = p.periods();
    uneq.n_obligors = n.data();
    uneq.n_defaults = d.data();
    auto a = sized<vcal_posterior>();
    VCAL_CHECK(vcal_calibrate_posterior(c.ctx, &uneq, &grid, VCAL_PRIOR_JEFFREYS, &a) == VCAL_E_UNSUPPORTED);
    VCAL_CHECK(vcal_calibrate_posterior(c.ctx, &uneq, &grid, VCAL_PRIOR_FLAT, &a) == VCAL_OK);
    VCAL_CHECK(vcal_calibrate_posterior(c.ctx, &p.abi, &grid, 5, &a) == VCAL_E_INVALID_ARGUMENT);
}
