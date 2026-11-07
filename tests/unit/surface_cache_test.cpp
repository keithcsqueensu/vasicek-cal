// SPDX-License-Identifier: Apache-2.0
// engine/surface_cache.hpp (D-167): rows cached across panels are identical to evaluate_surface's; a
// change in any one setting that determines a row is a cache miss, never a stale hit; rows are filled
// on demand; and a recovery fit through the cache equals the fit without it, bit for bit.
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/objectives/vasicek_rate.hpp"
#include "core/quadrature/composite_legendre.hpp"
#include "core/quadrature/gauss_hermite.hpp"
#include "core/quadrature/parity.hpp"
#include "engine/calibrate.hpp"
#include "engine/surface.hpp"
#include "engine/surface_cache.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

namespace e = vcal::engine;
namespace q = vcal::quadrature;
using Backend = vcal::backends::CpuBackend;
using Binomial = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
using BinomialMixed = vcal::objectives::BinomialMixture<vcal::PrecisionMixed>;

vcal::Grid<2> small_grid() { return {{{1e-3, 0.2, 9, vcal::AxisScale::Logit}, {1e-2, 0.5, 7, vcal::AxisScale::Logit}}}; }

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

// The row evaluate_surface computes for one observation under a configuration.
template <class Objective, class Integrator>
std::vector<double> direct_row(const Objective& o, const Integrator& in, const vcal::Grid<2>& g,
                               const typename Objective::Obs& y) {
    std::vector<double> L(static_cast<std::size_t>(g.size()));
    e::evaluate_surface(Backend{1}, o, in, &y, 1, g, L.data());
    return L;
}

}  // namespace

VCAL_TEST(cached_surfaces_equal_evaluate_surface_bit_for_bit) {
    const auto in = q::parity_rule();
    const auto g = small_grid();
    const std::vector<std::vector<Binomial::Obs>> panels = {
        {{1000, 3}, {1000, 12}, {1000, 3}, {1000, 0}}, {{1000, 12}, {1000, 7}, {1000, 0}}, {{500, 12}, {1000, 12}}};
    e::SurfaceRowCache cache;
    for (const auto& p : panels) {
        const auto T = static_cast<std::int64_t>(p.size());
        std::vector<double> a(static_cast<std::size_t>(T * g.size())), b(a.size());
        e::evaluate_surface(Backend{1}, Binomial{}, in, p.data(), T, g, a.data());
        e::evaluate_surface_cached(Backend{1}, cache, Binomial{}, in, p.data(), T, g, b.data());
        VCAL_CHECK(same_bits(a, b));
    }
    // Rows filled on demand, once each: (1000, 3), (1000, 12), (1000, 0), (1000, 7), (500, 12).
    VCAL_CHECK_EQ(cache.size(), 5u);
    VCAL_CHECK_EQ(cache.misses(), 5);
    VCAL_CHECK_EQ(cache.hits(), 3);  // (1000, 12) and (1000, 0) in panel 2, (1000, 12) in panel 3
}

// P-26 (D-409): evaluate_surface_cached writes every element, so a reused surface needs no clearing: over a
// buffer filled with a NaN pattern, and over a previous panel's surface of the same size, the result equals the
// zeroed fill's bit for bit; and calibrate_cached on a reused vector equals the fit on a fresh one.
VCAL_TEST(cached_surface_overwrites_every_element) {
    const auto in = q::parity_rule();
    const auto g = small_grid();
    const std::vector<std::vector<Binomial::Obs>> panels = {
        {{1000, 3}, {1000, 12}, {1000, 3}, {1000, 0}}, {{1000, 12}, {1000, 7}, {1000, 0}, {1000, 7}}};
    e::SurfaceRowCache cache;
    std::vector<double> reused;
    for (const auto& p : panels) {
        const auto T = static_cast<std::int64_t>(p.size());
        std::vector<double> zeroed(static_cast<std::size_t>(T * g.size()), 0.0);
        std::vector<double> junk(zeroed.size());
        const std::uint64_t pattern = 0x7FF8DEADBEEF0001ull;  // a quiet NaN with a payload
        for (double& x : junk) std::memcpy(&x, &pattern, sizeof x);
        e::evaluate_surface_cached(Backend{1}, cache, Binomial{}, in, p.data(), T, g, zeroed.data());
        e::evaluate_surface_cached(Backend{1}, cache, Binomial{}, in, p.data(), T, g, junk.data());
        reused.resize(zeroed.size());  // as calibrate_cached: the previous panel's values stay
        e::evaluate_surface_cached(Backend{1}, cache, Binomial{}, in, p.data(), T, g, reused.data());
        VCAL_CHECK(same_bits(zeroed, junk));
        VCAL_CHECK(same_bits(zeroed, reused));
    }
}

// Every setting that determines a row is part of its key: changing any one of them is a miss, and the
// row returned is the one that configuration computes.
VCAL_TEST(changing_any_row_setting_is_a_cache_miss) {
    const auto base_in = q::parity_rule();
    const auto base_g = small_grid();
    const Binomial::Obs y{1000, 12};
    e::SurfaceRowCache cache;
    const Backend b{1};
    const auto request = [&](const auto& objective, const auto& integrator, const vcal::Grid<2>& g, const auto& obs) {
        const auto config = e::surface_configuration(objective, integrator, g);
        const auto before = cache.misses();
        const auto& row = cache.row(b, config, objective, integrator, g, obs);
        return std::make_pair(cache.misses() - before, same_bits(row, direct_row(objective, integrator, g, obs)));
    };
    VCAL_CHECK(request(Binomial{}, base_in, base_g, y) == std::make_pair(std::int64_t{1}, true));
    VCAL_CHECK(request(Binomial{}, base_in, base_g, y) == std::make_pair(std::int64_t{0}, true));  // the same: a hit

    std::vector<std::pair<std::string, vcal::Grid<2>>> grids;
    for (int a = 0; a < 2; ++a) {
        auto g = base_g;
        g.axis[a].lo *= 1.5;
        grids.push_back({"lower bound " + std::to_string(a), g});
        g = base_g;
        g.axis[a].hi *= 0.9;
        grids.push_back({"upper bound " + std::to_string(a), g});
        g = base_g;
        g.axis[a].n += 2;
        grids.push_back({"spacing " + std::to_string(a), g});
        g = base_g;
        g.axis[a].scale = vcal::AxisScale::Probit;
        grids.push_back({"scale " + std::to_string(a), g});
    }
    for (const auto& [what, g] : grids) {
        const auto r = request(Binomial{}, base_in, g, y);
        if (r.first != 1 || !r.second) vcal::test::note("grid change not a clean miss: " + what);
        VCAL_CHECK(r.first == 1 && r.second);
    }
    // The integrator: Gauss-Hermite N, the composite rule's panels and its Legendre points, and the
    // doubled check rule.
    auto gh = base_in;
    gh.interior.rule = q::gauss_hermite_rule(64);
    auto panels = base_in;
    panels.one_sided.panels = 32;
    auto points = base_in;
    points.one_sided.rule = q::gauss_legendre_rule(12);
    for (const auto& in : {gh, panels, points, q::parity_rule(true)}) {
        const auto r = request(Binomial{}, in, base_g, y);
        VCAL_CHECK(r.first == 1 && r.second);
    }
    // The precision policy (a different objective type).
    const auto rm = request(BinomialMixed{}, q::SplitRule<vcal::PrecisionMixed>{
                                                 q::GaussHermiteAdaptive<vcal::PrecisionMixed>{base_in.interior.rule},
                                                 q::CompositeLegendre<vcal::PrecisionMixed>{base_in.one_sided.rule,
                                                                                            base_in.one_sided.panels}},
                            base_g, BinomialMixed::Obs{1000, 12});
    VCAL_CHECK(rm.first == 1);
    // The zero/all-default rule: the rate objective's treatments, on the same observation.
    using namespace vcal::objectives;
    const RateObs z{0.0, 5e-3};
    VCAL_CHECK(request(VasicekRate<vcal::PrecisionF64, ZeroRates::Censor>{}, base_in, base_g, z).first == 1);
    VCAL_CHECK(request(VasicekRate<vcal::PrecisionF64, ZeroRates::Substitute>{}, base_in, base_g, z).first == 1);
    VCAL_CHECK(request(VasicekRate<vcal::PrecisionF64, ZeroRates::Censor>{}, base_in, base_g, z).first == 0);
    // ... and a different detection limit is a different observation.
    VCAL_CHECK(request(VasicekRate<vcal::PrecisionF64, ZeroRates::Censor>{}, base_in, base_g, RateObs{0.0, 5e-4}).first == 1);
}

// A recovery fit through the cache equals the fit without it, bit for bit, field by field.
VCAL_TEST(recovery_fits_through_the_cache_are_identical) {
    namespace rc = vcal::recovery;
    e::SurfaceRowCache cache;
    int differing = 0;
    for (const std::uint32_t id : {29u, 55u}) {  // T = 20, n = 10^4 and 10^3
        for (std::uint32_t r = 0; r < 2; ++r) {
            const auto a = rc::fit(rc::scenario(id), r);
            const auto c = rc::fit(rc::scenario(id), r, &cache);
            const auto eq = [](const auto& x, const auto& y) { return std::memcmp(&x, &y, sizeof x) == 0; };
            const bool same = eq(a.value, c.value) && eq(a.se, c.se) && eq(a.loglik, c.loglik) && a.flags == c.flags &&
                              eq(a.prof_lo, c.prof_lo) && eq(a.prof_hi, c.prof_hi) && eq(a.prof_flags, c.prof_flags) &&
                              eq(a.boot_lo, c.boot_lo) && eq(a.boot_hi, c.boot_hi) && a.boot_edge == c.boot_edge &&
                              eq(a.q_hat, c.q_hat) && eq(a.q_se_s, c.q_se_s) && eq(a.q_prof_lo, c.q_prof_lo) &&
                              eq(a.q_prof_hi, c.q_prof_hi) && a.q_prof_flags == c.q_prof_flags &&
                              eq(a.q_boot_lo, c.q_boot_lo) && eq(a.q_boot_hi, c.q_boot_hi);
            differing += same ? 0 : 1;
        }
    }
    VCAL_CHECK_EQ(differing, 0);
    VCAL_CHECK(cache.hits() > 0);  // replicates of one scenario share rows
}

// P-27 (D-410): calibrate_from_rows (the cache's rows in place, the argmax bounded by the rows' cached tile maxima)
// gives calibrate_from_surface's estimate on the copied rows, bit for bit, in every field:
//   - on recovery panels of every scenario (two replicates each) and on bootstrap-like panels drawn at their
//     estimates, through the cache;
//   - on adversarial finite surfaces (flat stretches, ties, huge and tiny values), where the bound visits many tiles;
//   - on surfaces with an infinite or NaN value, which take the full grid (D-172's guard) on the copied rows.
// The cache's row tile maxima equal tile_maxima on the row, and its finiteness flag is the row's.
VCAL_TEST(calibrate_from_rows_is_calibrate_from_surface_bit_for_bit) {
    namespace rc = vcal::recovery;
    using Obs = rc::Objective::Obs;
    const auto primary = q::parity_rule();
    const auto check = q::parity_rule(true);
    const auto g = rc::grid();
    const std::int64_t K = g.size();
    auto same_est = [](const e::Estimate2& a, const e::Estimate2& b) {
        auto eq = [](double x, double y) { return std::memcmp(&x, &y, sizeof x) == 0; };
        bool s = a.grid_index == b.grid_index && a.flags == b.flags && a.nan_count == b.nan_count &&
                 eq(a.loglik, b.loglik) && eq(a.corr, b.corr) && eq(a.quad_check_max, b.quad_check_max) &&
                 eq(a.quad_check_total, b.quad_check_total) && a.quad_check_flagged == b.quad_check_flagged;
        for (int i = 0; i < 2; ++i) s = s && eq(a.value[i], b.value[i]) && eq(a.se[i], b.se[i]) && eq(a.scaled[i], b.scaled[i]);
        for (int i = 0; i < 3; ++i) s = s && eq(a.hess_u[i], b.hess_u[i]);
        return s;
    };
    e::SurfaceRowCache cache;
    int fits = 0, bad = 0, tiles_bad = 0;
    auto check_panel = [&](const std::vector<Obs>& obs, bool qc) {
        const auto T = static_cast<std::int64_t>(obs.size());
        std::vector<double> L(static_cast<std::size_t>(T * K)), Lr;
        e::evaluate_surface_cached(Backend{1}, cache, rc::Objective{}, primary, obs.data(), T, g, L.data());
        e::SurfaceRows rows;
        e::surface_rows_cached(Backend{1}, cache, rc::Objective{}, primary, obs.data(), T, g, rows);
        for (std::int64_t t = 0; t < T; ++t) {
            std::vector<double> U(static_cast<std::size_t>(e::argmax_tiles(g)));
            e::tile_maxima(g, rows.row[static_cast<std::size_t>(t)], U.data());
            tiles_bad += std::memcmp(U.data(), rows.tile_max[static_cast<std::size_t>(t)], U.size() * sizeof(double)) == 0 ? 0 : 1;
        }
        e::Estimate2 a{}, b{};
        const auto sa = e::calibrate_from_surface(Backend{1}, rc::Objective{}, primary, check, obs.data(), T, g, L, a,
                                                  e::kHessianStepFraction, qc);
        const auto sb = e::calibrate_from_rows(Backend{1}, rc::Objective{}, primary, check, obs.data(), T, g, rows, Lr, b,
                                               e::kHessianStepFraction, qc);
        ++fits;
        bad += (sa == sb && (sa != e::Status::Ok || same_est(a, b))) ? 0 : 1;
        return a;
    };
    for (std::uint32_t id = 0; id < rc::kScenarios; ++id) {
        const auto s = rc::scenario(id);
        for (std::uint32_t r = 0; r < 2; ++r) {
            std::vector<Obs> obs;
            for (const auto d : rc::panel(s, r)) obs.push_back({s.obligors, d});
            const auto est = check_panel(obs, r == 0);
            // A bootstrap-like panel: drawn at the estimate, as S-10 draws its parametric bootstrap.
            if (std::isfinite(est.value[0]) && std::isfinite(est.value[1])) {
                std::vector<std::int64_t> n(obs.size(), s.obligors), d(obs.size());
                const vcal::dgp::PanelSpec spec{0x5031u, id, 7u + r, est.value[0], est.value[1], n.data(),
                                                static_cast<std::int64_t>(obs.size())};
                if (vcal::dgp::simulate_panel(spec, d.data(), nullptr) == vcal::dgp::Status::Ok) {
                    std::vector<Obs> ob;
                    for (const auto x : d) ob.push_back({s.obligors, x});
                    check_panel(ob, false);
                }
            }
        }
    }
    vcal::test::note(std::to_string(fits) + " cached panels, " + std::to_string(bad) + " estimates differ, " +
                     std::to_string(tiles_bad) + " tile-maxima rows differ");
    VCAL_CHECK(bad == 0 && tiles_bad == 0);

    // Hand-made surfaces through row pointers: adversarial finite ones, and ones with a non-finite value.
    std::uint64_t st = 0x27u;
    auto uni = [&] {
        st = st * 6364136223846793005ull + 1442695040888963407ull;
        return static_cast<double>(st >> 11) * 0x1p-53;
    };
    int made = 0, made_bad = 0;
    for (int trial = 0; trial < 60; ++trial) {
        const std::int64_t T = 1 + trial % 13;
        std::vector<double> L(static_cast<std::size_t>(T * K));
        const int kind = trial % 4;  // 0 smooth, 1 flat stretches and ties, 2 huge and tiny, 3 a non-finite value
        for (std::int64_t t = 0; t < T; ++t) {
            for (std::int64_t k = 0; k < K; ++k) {
                const double x = static_cast<double>(k % 61) - 30.0, y = static_cast<double>(k / 61) - 20.0;
                double v = -(x * x + 0.7 * y * y) * (1.0 + 0.1 * static_cast<double>(t)) - 3.0 * uni();
                if (kind == 1) v = -std::floor(std::fabs(x) / 7.0) - std::floor(std::fabs(y) / 5.0);
                if (kind == 2 && uni() < 0.2) v = uni() < 0.5 ? -1e300 * uni() : -1e-300 * uni();
                L[static_cast<std::size_t>(t * K + k)] = v;
            }
        }
        if (kind == 3) {
            const double bad_values[] = {-HUGE_VAL, HUGE_VAL, std::nan("")};
            L[static_cast<std::size_t>((trial % static_cast<int>(T)) * K + (trial * 37) % K)] = bad_values[trial % 3];
        }
        e::SurfaceRows rows;
        rows.periods = T;
        rows.K = K;
        std::vector<std::vector<double>> U(static_cast<std::size_t>(T));
        for (std::int64_t t = 0; t < T; ++t) {
            const double* row = L.data() + t * K;
            U[static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(e::argmax_tiles(g)));
            e::tile_maxima(g, row, U[static_cast<std::size_t>(t)].data());
            rows.row.push_back(row);
            rows.tile_max.push_back(U[static_cast<std::size_t>(t)].data());
            for (std::int64_t k = 0; k < K; ++k) rows.finite = rows.finite && std::isfinite(row[k]);
        }
        std::vector<Obs> obs(static_cast<std::size_t>(T), Obs{1000, 3});
        std::vector<double> Lr;
        e::Estimate2 a{}, b{};
        const auto sa = e::calibrate_from_surface(Backend{1}, rc::Objective{}, primary, check, obs.data(), T, g, L, a);
        const auto sb = e::calibrate_from_rows(Backend{1}, rc::Objective{}, primary, check, obs.data(), T, g, rows, Lr, b);
        ++made;
        made_bad += (sa == sb && (sa != e::Status::Ok || same_est(a, b))) ? 0 : 1;
        if (kind == 3) VCAL_CHECK(!rows.finite && Lr.size() == static_cast<std::size_t>(T * K));
    }
    vcal::test::note(std::to_string(made) + " hand-made surfaces (15 with a non-finite value), " +
                     std::to_string(made_bad) + " differ");
    VCAL_CHECK(made_bad == 0);
}
