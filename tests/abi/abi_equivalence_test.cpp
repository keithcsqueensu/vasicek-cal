// SPDX-License-Identifier: Apache-2.0
//
// The C ABI against the engine (M2c, D-144): every vcal_* result must be bitwise what the engine
// gives when called directly with the same inputs. Together with the recovery replay (engine
// against goldens) this ties the ABI to the pinned results. The C test (abi_c_test.c) covers
// the conventions a C caller sees; this one covers the numbers.
#include <cstdint>
#include <cstring>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/quadrature/parity.hpp"
#include "dgp/dgp.hpp"
#include "engine/calibrate.hpp"
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
