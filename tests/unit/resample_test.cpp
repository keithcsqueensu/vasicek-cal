// SPDX-License-Identifier: Apache-2.0
// resample/ (M2b): weight matrices, replicate estimates from W x L, percentile intervals
// (D-132..D-135); jackknife bias correction, BCa and leave-two-out weights (S-3, S-5, S-21; D-155).
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/gauss_hermite.hpp"
#include "dgp/philox.hpp"
#include "engine/calibrate.hpp"
#include "resample/bootstrap.hpp"
#include "resample/jackknife.hpp"
#include "resample/weights.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace e = vcal::engine;
namespace rs = vcal::resample;
namespace tol = vcal::tol;
using vcal::test::describe;
using Backend = vcal::backends::CpuBackend;
using Objective = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
using Adaptive = vcal::quadrature::GaussHermiteAdaptive<vcal::PrecisionF64>;
using Obs = Objective::Obs;

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

bool same_replicate(const rs::Replicate2& a, const rs::Replicate2& b) {
    return bits_equal(a.value[0], b.value[0]) && bits_equal(a.value[1], b.value[1]) &&
           bits_equal(a.surface_max, b.surface_max) && a.grid_index == b.grid_index && a.flags == b.flags;
}

// The M1.4 test panel: 20 periods of 1000 obligors.
std::vector<Obs> panel() {
    const std::int64_t d[] = {2, 5, 1, 9, 3, 0, 4, 12, 6, 2, 1, 3, 7, 15, 2, 4, 0, 5, 8, 3};
    std::vector<Obs> p;
    for (const auto x : d) p.push_back({1000, x});
    return p;
}

vcal::Grid<2> grid() {
    return {{{2e-4, 0.1, 41, vcal::AxisScale::Logit}, {5e-3, vcal::kDefaultRhoUpper, 31, vcal::AxisScale::Logit}}};
}

struct Rules {
    Adaptive primary{vcal::quadrature::gauss_hermite_rule(64)};
    Adaptive check{vcal::quadrature::gauss_hermite_rule(128)};
};

e::Estimate2 calibrate(const std::vector<Obs>& p, std::vector<double>& L) {
    const Rules r;
    e::Estimate2 est{};
    VCAL_REQUIRE(e::calibrate(Backend{}, Objective{}, r.primary, r.check, p.data(), static_cast<std::int64_t>(p.size()),
                              grid(), L, est) == e::Status::Ok);
    return est;
}

constexpr std::uint64_t kSeed = 0x5641534943454B31ull;  // the mirror's REFERENCE_SEED

}  // namespace

// --- draws: bit for bit against the Python mirror (tools/gen_dgp_tables.py) ---------------------

VCAL_TEST(resample_draws_match_the_python_mirror) {
    const auto t = vcal::test::read_golden_csv("dgp/resample_reference.csv");
    constexpr std::int64_t T = 23, ell = 5;
    constexpr std::uint32_t B = 4;
    const auto iid = rs::bootstrap_indices(kSeed, rs::Scheme::IidBootstrap, B, T);
    const auto block = rs::bootstrap_indices(kSeed, rs::Scheme::BlockBootstrap, B, T, ell);
    int rows = 0;
    for (const auto& r : t.rows) {
        const bool is_iid = r[t.column("scheme")] == "iid";
        const auto scheme = is_iid ? rs::Scheme::IidBootstrap : rs::Scheme::BlockBootstrap;
        const auto b = static_cast<std::uint32_t>(vcal::test::parse_int(r[t.column("replicate")]));
        const auto j = static_cast<std::uint32_t>(vcal::test::parse_int(r[t.column("draw")]));
        const double u = rs::draw_uniform(kSeed, scheme, b, j);
        VCAL_CHECK(bits_equal(u, vcal::test::parse_double(r[t.column("u_hex")])));
        VCAL_CHECK_EQ(rs::draw_index(kSeed, scheme, b, j, is_iid ? T : T - ell + 1),
                      vcal::test::parse_int(r[t.column("value")]));
        const auto pos = vcal::test::parse_int(r[t.column("position")]);
        const auto& idx = is_iid ? iid : block;
        VCAL_CHECK_EQ(idx[static_cast<std::size_t>(b * T + pos)], vcal::test::parse_int(r[t.column("index")]));
        ++rows;
    }
    VCAL_CHECK_EQ(rows, static_cast<int>(2 * B * T));
}

// Resampling has a key domain of its own: for the same seed, no DGP counter and no resampling
// counter share a key, so resampling a simulated panel cannot reuse its stream (D-132).
VCAL_TEST(resample_key_domain_differs_from_the_dgp) {
    const vcal::dgp::PhiloxKey dgp_key{static_cast<std::uint32_t>(kSeed), static_cast<std::uint32_t>(kSeed >> 32)};
    const auto rk = rs::resample_key(kSeed);
    VCAL_CHECK(rk[0] != dgp_key[0] && rk[1] != dgp_key[1]);
    const auto a = vcal::dgp::philox4x32({1u, 0u, 0u, 0u}, dgp_key);
    const auto b = vcal::dgp::philox4x32({1u, 0u, 0u, 0u}, rk);
    VCAL_CHECK(a != b);
}

// --- weight matrices ---------------------------------------------------------------------------

VCAL_TEST(resample_weight_matrices_have_the_right_shape) {
    constexpr std::int64_t T = 23;
    constexpr std::uint32_t B = 50;
    for (const auto scheme : {rs::Scheme::IidBootstrap, rs::Scheme::BlockBootstrap}) {
        const auto idx = rs::bootstrap_indices(kSeed, scheme, B, T, 5);
        const auto w = rs::weights_from_indices(idx.data(), B, T, T);
        for (std::uint32_t b = 0; b < B; ++b) {
            double total = 0.0;
            for (std::int64_t t = 0; t < T; ++t) total += w[static_cast<std::size_t>(b * T + t)];
            VCAL_CHECK_EQ(total, static_cast<double>(T));  // every replicate has T periods
            if (scheme == rs::Scheme::BlockBootstrap) {
                // Blocks of 5 consecutive periods, the last one cut at T.
                for (std::int64_t j = 0; j < T; ++j) {
                    if (j % 5 != 0) {
                        VCAL_CHECK_EQ(idx[static_cast<std::size_t>(b * T + j)], idx[static_cast<std::size_t>(b * T + j - 1)] + 1);
                    }
                }
            }
        }
    }
    const auto jk = rs::jackknife_weights(4);
    const std::vector<double> want_jk = {0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0};
    VCAL_CHECK(jk == want_jk);
    // Rolling windows of 3 ending at 3, 5 (T = 6, step 2); expanding windows ending at 4, 6.
    const std::vector<double> rolling = {1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0};
    VCAL_CHECK(rs::walk_forward_weights(6, 3, 3, 2) == rolling);
    const std::vector<double> expanding = {1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1};
    VCAL_CHECK(rs::walk_forward_weights(6, 0, 4, 2) == expanding);
    bool threw = false;
    try {
        const std::int64_t bad[2] = {0, 6};
        (void)rs::weights_from_indices(bad, 1, 2, 6);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    VCAL_CHECK(threw);
}

// --- replicate estimates ----------------------------------------------------------------------

// The fused, tiled reduction gives exactly what materialising each replicate's surface would,
// for any thread count.
VCAL_TEST(replicate_estimates_match_a_materialised_reduction) {
    const auto p = panel();
    std::vector<double> L;
    (void)calibrate(p, L);
    const auto g = grid();
    const std::int64_t T = static_cast<std::int64_t>(p.size()), K = g.size();
    constexpr std::uint32_t B = 64;
    const auto idx = rs::bootstrap_indices(kSeed, rs::Scheme::IidBootstrap, B, T);
    const auto W = rs::weights_from_indices(idx.data(), B, T, T);
    std::vector<rs::Replicate2> naive(B);
    for (std::uint32_t b = 0; b < B; ++b) {
        std::vector<double> s(static_cast<std::size_t>(K));
        std::int64_t best = 0;
        for (std::int64_t k = 0; k < K; ++k) {
            s[static_cast<std::size_t>(k)] = e::weighted_sum(L.data(), T, K, W.data() + b * T, k);
            if (s[static_cast<std::size_t>(k)] > s[static_cast<std::size_t>(best)]) best = k;
        }
        const auto f = e::refine_2d(g, best, [&](std::int32_t i0, std::int32_t i1) {
            const std::int32_t ii[2] = {i0, i1};
            return s[static_cast<std::size_t>(g.flatten(ii))];
        });
        naive[b] = {{f.value[0], f.value[1]}, s[static_cast<std::size_t>(best)], best, f.flags};
    }
    for (const int threads : {1, 3, 0}) {
        std::vector<rs::Replicate2> fused(B);
        rs::replicate_estimates(Backend{threads}, g, L.data(), T, W.data(), B, fused.data());
        bool same = true;
        for (std::uint32_t b = 0; b < B; ++b) same = same && same_replicate(fused[b], naive[b]);
        VCAL_CHECK(same);
    }
}

// A replicate's estimate is what calibrate gives on the panel the replicate stands for (its
// periods repeated as drawn). Only summation order differs: count x l_t against l_t repeated.
VCAL_TEST(replicate_estimate_equals_calibrate_on_the_resampled_panel) {
    const auto p = panel();
    std::vector<double> L;
    (void)calibrate(p, L);
    const auto g = grid();
    const std::int64_t T = static_cast<std::int64_t>(p.size());
    constexpr std::uint32_t B = 8;
    for (const auto scheme : {rs::Scheme::IidBootstrap, rs::Scheme::BlockBootstrap}) {
        const auto idx = rs::bootstrap_indices(kSeed, scheme, B, T, 4);
        const auto W = rs::weights_from_indices(idx.data(), B, T, T);
        std::vector<rs::Replicate2> reps(B);
        rs::replicate_estimates(Backend{}, g, L.data(), T, W.data(), B, reps.data());
        double worst = 0.0;
        for (std::uint32_t b = 0; b < B; ++b) {
            std::vector<Obs> resampled;
            for (std::int64_t j = 0; j < T; ++j) resampled.push_back(p[static_cast<std::size_t>(idx[static_cast<std::size_t>(b * T + j)])]);
            std::vector<double> L2;
            const auto est = calibrate(resampled, L2);
            VCAL_CHECK_EQ(reps[b].grid_index, est.grid_index);
            const std::uint32_t replicate_flags = e::kFlagGridEdge | e::kFlagFlatSurface | e::kFlagRefinementRejected;
            VCAL_CHECK_EQ(reps[b].flags & replicate_flags, est.flags & replicate_flags);
            for (int a = 0; a < 2; ++a) worst = std::fmax(worst, std::fabs(reps[b].value[a] / est.value[a] - 1.0));
        }
        vcal::test::note(std::string(scheme == rs::Scheme::IidBootstrap ? "iid" : "moving block") +
                         ": worst relative difference, replicate vs refit " + describe(worst));
        VCAL_CHECK(worst <= tol::TOL_RESAMPLE_VS_REFIT_REL);
    }
}

// Externally supplied indices reproduce generated replicates exactly (D-134).
VCAL_TEST(external_indices_reproduce_generated_replicates) {
    const auto p = panel();
    std::vector<double> L;
    (void)calibrate(p, L);
    const std::int64_t T = static_cast<std::int64_t>(p.size());
    constexpr std::uint32_t B = 16;
    const auto idx = rs::bootstrap_indices(kSeed, rs::Scheme::BlockBootstrap, B, T, 3);
    // "Another tool" hands over the same indices; the weights, and so the replicates, are identical.
    const std::vector<std::int64_t> external(idx.begin(), idx.end());
    const auto W1 = rs::weights_from_indices(idx.data(), B, T, T);
    const auto W2 = rs::weights_from_indices(external.data(), B, T, T);
    std::vector<rs::Replicate2> a(B), b(B);
    rs::replicate_estimates(Backend{}, grid(), L.data(), T, W1.data(), B, a.data());
    rs::replicate_estimates(Backend{}, grid(), L.data(), T, W2.data(), B, b.data());
    for (std::uint32_t r = 0; r < B; ++r) VCAL_CHECK(same_replicate(a[r], b[r]));
}

// --- percentile intervals ---------------------------------------------------------------------

VCAL_TEST(percentile_interval_is_type_7) {
    std::vector<rs::Replicate2> reps(11);
    for (int i = 0; i < 10; ++i) reps[static_cast<std::size_t>(i)] = {{static_cast<double>(10 - i), 0.0}, 0.0, 0, 0u};
    reps[10] = {{std::nan(""), 0.0}, 0.0, -1, e::kFlagNumeric};
    // Values 1..10 (unsorted), one NaN: at level 0.8, h = 9 x 0.1 = 0.9 and 9 x 0.9 = 8.1.
    const auto ci = rs::percentile_interval(reps.data(), 11, 0, 0.8);
    VCAL_CHECK_REL(ci.lo, 1.9, 1e-15);
    VCAL_CHECK_REL(ci.hi, 9.1, 1e-15);
    VCAL_CHECK_EQ(ci.used, 10);
    VCAL_CHECK_EQ(ci.excluded, 1);
}

// --- the shared jackknife run (S-3, S-5, S-21; D-155) -----------------------------------------------

// Leave-two-out rows: every pair once, in lexicographic order, with both periods deleted.
VCAL_TEST(delete_two_weights_cover_every_pair_once) {
    const std::int64_t T = 5;
    const auto w = rs::delete_two_weights(T);
    VCAL_REQUIRE(static_cast<std::int64_t>(w.size()) == T * (T - 1) / 2 * T);
    std::int64_t row = 0;
    for (std::int64_t s = 0; s < T; ++s) {
        for (std::int64_t t = s + 1; t < T; ++t, ++row) {
            for (std::int64_t u = 0; u < T; ++u) {
                VCAL_CHECK_EQ(w[static_cast<std::size_t>(row * T + u)], (u == s || u == t) ? 0.0 : 1.0);
            }
        }
    }
}

// Quenouille's correction turns the plug-in variance (divisor T) into the unbiased one (divisor
// T - 1) exactly: the textbook identity, up to rounding.
VCAL_TEST(jackknife_bias_correction_of_the_plug_in_variance) {
    const std::vector<double> x = {0.3, 1.7, 2.2, 2.9, 3.4, 3.8, 4.1, 5.0, 6.6, 9.3};
    const auto T = static_cast<std::int64_t>(x.size());
    const auto plug_in = [&](std::int64_t skip) {
        double m = 0.0, n = 0.0;
        for (std::int64_t i = 0; i < T; ++i) {
            if (i != skip) { m += x[static_cast<std::size_t>(i)]; n += 1.0; }
        }
        m /= n;
        double v = 0.0;
        for (std::int64_t i = 0; i < T; ++i) {
            if (i != skip) v += (x[static_cast<std::size_t>(i)] - m) * (x[static_cast<std::size_t>(i)] - m);
        }
        return v / n;
    };
    std::vector<double> minus(static_cast<std::size_t>(T));
    for (std::int64_t t = 0; t < T; ++t) minus[static_cast<std::size_t>(t)] = plug_in(t);
    const double corrected = rs::jackknife_bias_corrected(plug_in(-1), minus.data(), T);
    const double unbiased = plug_in(-1) * static_cast<double>(T) / static_cast<double>(T - 1);
    vcal::test::note("corrected " + describe(corrected) + ", unbiased " + describe(unbiased));
    VCAL_CHECK_REL(corrected, unbiased, tol::TOL_RESAMPLE_JACKKNIFE_REL);
}

// BCa against scipy/numpy on a small case (norm.ppf, norm.cdf, numpy.quantile's default type 7):
// z0 = -0.2533471031357997, a = -0.059709802247422554, ends 0.3318271240904089 and 6.928724371446333.
VCAL_TEST(bca_interval_matches_scipy) {
    const std::vector<double> x = {9.3, 0.3, 1.7, 2.2, 2.9, std::nan(""), 3.4, 3.8, 4.1, 5.0, 6.6};
    const std::vector<double> jack = {3.1, 3.3, 3.45, 3.6, 4.2};
    const auto ci = rs::bca_interval(x.data(), static_cast<std::int64_t>(x.size()), 3.0, jack.data(), 5);
    VCAL_REQUIRE(ci.computed);
    const double want[4] = {-0.2533471031357997, -0.059709802247422554, 0.3318271240904089, 6.928724371446333};
    const double got[4] = {ci.z0, ci.acceleration, ci.lo, ci.hi};
    double worst = 0.0;
    for (int k = 0; k < 4; ++k) worst = std::fmax(worst, std::fabs(got[k] / want[k] - 1.0));
    vcal::test::note("BCa against scipy: worst relative difference " + describe(worst));
    VCAL_CHECK(worst <= tol::TOL_RESAMPLE_JACKKNIFE_REL);
    // With z0 = 0 and a = 0 BCa is the percentile interval.
    const std::vector<double> sym = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    const std::vector<double> flat = {1.0, 2.0, 3.0};  // symmetric: a = 0
    const auto p = rs::bca_interval(sym.data(), 9, 5.0, flat.data(), 3);
    VCAL_REQUIRE(p.computed);
    VCAL_CHECK_EQ(p.z0, 0.0);
    VCAL_CHECK_EQ(p.acceleration, 0.0);
    std::vector<double> sorted = sym;
    VCAL_CHECK_REL(p.lo, rs::quantile_type7(sorted, 0.025), tol::TOL_RESAMPLE_JACKKNIFE_REL);
    VCAL_CHECK_REL(p.hi, rs::quantile_type7(sorted, 0.975), tol::TOL_RESAMPLE_JACKKNIFE_REL);
}

// Not computed: every replicate above the estimate (z0 infinite), or no jackknife spread.
VCAL_TEST(bca_interval_not_computed_at_a_boundary) {
    const std::vector<double> above = {2.0, 3.0, 4.0};
    const std::vector<double> jack = {1.0, 2.0, 3.0};
    const auto a = rs::bca_interval(above.data(), 3, 1.0, jack.data(), 3);
    VCAL_CHECK(!a.computed && std::isnan(a.lo) && std::isnan(a.hi));
    const std::vector<double> same = {2.0, 2.0, 2.0};
    const auto b = rs::bca_interval(above.data(), 3, 3.0, same.data(), 3);
    VCAL_CHECK(!b.computed);
}

// D-172: the reduction over distinct rows in (n, d) order. Its order is defined by the
// observations, not by the panel's layout: reversing the periods (and W's columns with them)
// gives the same replicates bit for bit, for any thread count.
VCAL_TEST(compact_replicates_do_not_depend_on_period_order_or_threads) {
    const auto p = panel();
    std::vector<double> L;
    (void)calibrate(p, L);
    const auto g = grid();
    const std::int64_t T = static_cast<std::int64_t>(p.size()), K = g.size();
    constexpr std::uint32_t B = 64;
    const auto idx = rs::bootstrap_indices(kSeed, rs::Scheme::IidBootstrap, B, T);
    const auto W = rs::weights_from_indices(idx.data(), B, T, T);
    std::vector<rs::Replicate2> ref(B);
    rs::replicate_estimates_compact(Backend{1}, g, p.data(), L.data(), T, W.data(), B, ref.data());
    std::vector<Obs> rp(p.rbegin(), p.rend());
    std::vector<double> rL(L.size()), rW(W.size());
    for (std::int64_t t = 0; t < T; ++t) {
        std::copy(L.begin() + (T - 1 - t) * K, L.begin() + (T - t) * K, rL.begin() + t * K);
        for (std::uint32_t b = 0; b < B; ++b) rW[static_cast<std::size_t>(b * T + t)] = W[static_cast<std::size_t>(b * T + T - 1 - t)];
    }
    for (const int threads : {1, 3, 0}) {
        std::vector<rs::Replicate2> out(B), rev(B);
        rs::replicate_estimates_compact(Backend{threads}, g, p.data(), L.data(), T, W.data(), B, out.data());
        rs::replicate_estimates_compact(Backend{threads}, g, rp.data(), rL.data(), T, rW.data(), B, rev.data());
        bool same = true;
        for (std::uint32_t b = 0; b < B; ++b) same = same && same_replicate(out[b], ref[b]) && same_replicate(rev[b], ref[b]);
        VCAL_CHECK(same);
    }
}

// The distinct-row reduction changes the sums only by rounding: fewer and reordered terms. Against
// the reduction over periods, the same argmax and estimates within TOL_RESAMPLE_VS_REFIT_REL.
VCAL_TEST(compact_replicates_match_the_period_reduction_to_rounding) {
    const auto p = panel();
    std::vector<double> L;
    (void)calibrate(p, L);
    const auto g = grid();
    const std::int64_t T = static_cast<std::int64_t>(p.size());
    constexpr std::uint32_t B = 256;
    const auto idx = rs::bootstrap_indices(kSeed, rs::Scheme::IidBootstrap, B, T);
    const auto W = rs::weights_from_indices(idx.data(), B, T, T);
    std::vector<rs::Replicate2> plain(B), compact(B);
    rs::replicate_estimates(Backend{}, g, L.data(), T, W.data(), B, plain.data());
    rs::replicate_estimates_compact(Backend{}, g, p.data(), L.data(), T, W.data(), B, compact.data());
    double worst = 0.0;
    std::int64_t bitwise = 0;
    for (std::uint32_t b = 0; b < B; ++b) {
        VCAL_CHECK_EQ(compact[b].grid_index, plain[b].grid_index);
        VCAL_CHECK_EQ(compact[b].flags, plain[b].flags);
        for (int a = 0; a < 2; ++a) worst = std::fmax(worst, std::fabs(compact[b].value[a] / plain[b].value[a] - 1.0));
        bitwise += same_replicate(compact[b], plain[b]) ? 1 : 0;
    }
    vcal::test::note("worst relative difference " + describe(worst) + "; bit-identical replicates " +
                     std::to_string(bitwise) + " of " + std::to_string(B));
    VCAL_CHECK(worst <= tol::TOL_RESAMPLE_VS_REFIT_REL);
}

// D-172: the bounded argmax equals the full grid's, bit for bit, on the surfaces where a coarse
// search is most likely to go wrong: rho not identified (n = 1 in every period, S-1's panel 26),
// estimates on a bound of the box (no defaults: PD at its lower end; extreme dispersion: rho at
// its upper end), a single distinct row, and the M1.4 panel; bootstrap and jackknife weights; any
// thread count. The at-bound panels must actually put replicates on the grid edge.
VCAL_TEST(bounded_argmax_equals_the_full_grid_on_hard_panels) {
    struct Case {
        const char* name;
        std::vector<Obs> p;
        bool expect_edge;
    };
    std::vector<Case> cases;
    cases.push_back({"M1.4 panel", panel(), false});
    {
        std::vector<Obs> p;
        for (int t = 0; t < 40; ++t) p.push_back({1, t % 7 == 0 ? 1 : 0});
        cases.push_back({"n = 1 in every period (rho not identified)", p, true});
    }
    cases.push_back({"no defaults (PD on its lower bound)", std::vector<Obs>(20, Obs{1000, 0}), true});
    {
        std::vector<Obs> p;
        const std::int64_t d[] = {0, 0, 0, 300, 0, 0, 250, 0, 1, 0, 0, 400, 0, 0, 0, 200, 0, 0, 0, 350};
        for (const auto x : d) p.push_back({1000, x});
        cases.push_back({"extreme dispersion (rho on its upper bound)", p, true});
    }
    cases.push_back({"one distinct row", std::vector<Obs>(25, Obs{500, 5}), false});
    const auto g = grid();
    for (const auto& c : cases) {
        std::vector<double> L;
        (void)calibrate(c.p, L);
        const std::int64_t T = static_cast<std::int64_t>(c.p.size());
        constexpr std::uint32_t B = 256;
        const auto idx = rs::bootstrap_indices(kSeed, rs::Scheme::IidBootstrap, B, T);
        const auto boot = rs::weights_from_indices(idx.data(), B, T, T);
        const auto jack = rs::jackknife_weights(T);
        std::int64_t compared = 0, edge = 0;
        bool same = true;
        for (const auto* W : {&boot, &jack}) {
            const std::int64_t R = static_cast<std::int64_t>(W->size()) / T;
            std::vector<rs::Replicate2> full(static_cast<std::size_t>(R));
            rs::replicate_estimates_compact(Backend{1}, g, c.p.data(), L.data(), T, W->data(), R, full.data(),
                                            rs::ArgmaxSearch::Full);
            for (const int threads : {1, 3, 0}) {
                std::vector<rs::Replicate2> bounded(static_cast<std::size_t>(R));
                rs::replicate_estimates_compact(Backend{threads}, g, c.p.data(), L.data(), T, W->data(), R, bounded.data(),
                                                rs::ArgmaxSearch::Bounded);
                for (std::int64_t b = 0; b < R; ++b) same = same && same_replicate(bounded[static_cast<std::size_t>(b)], full[static_cast<std::size_t>(b)]);
            }
            for (const auto& r : full) edge += (r.flags & e::kFlagGridEdge) ? 1 : 0;
            compared += R;
        }
        vcal::test::note(std::string(c.name) + ": " + std::to_string(compared) + " replicates, " + std::to_string(edge) +
                         " on the grid edge");
        VCAL_CHECK(same);
        if (c.expect_edge) VCAL_CHECK(edge > 0);
    }
}

// The bounded search's guard: a non-finite surface value or a negative weight sends the call to
// the full grid (argmax_bounded declines), so NaN counts and every other result stay the full grid's.
VCAL_TEST(bounded_argmax_guard_falls_back_to_the_full_grid) {
    const auto p = panel();
    std::vector<double> L;
    (void)calibrate(p, L);
    const auto g = grid();
    const std::int64_t T = static_cast<std::int64_t>(p.size());
    constexpr std::uint32_t B = 16;
    const auto idx = rs::bootstrap_indices(kSeed, rs::Scheme::IidBootstrap, B, T);
    auto W = rs::weights_from_indices(idx.data(), B, T, T);
    std::vector<vcal::reducers::ArgMax::State> best(B);
    VCAL_CHECK(rs::detail::argmax_bounded(Backend{}, g, L.data(), T, W.data(), B, best.data()));
    auto L_nan = L;
    L_nan[123] = std::nan("");
    VCAL_CHECK(!rs::detail::argmax_bounded(Backend{}, g, L_nan.data(), T, W.data(), B, best.data()));
    std::vector<rs::Replicate2> a(B), b(B);
    rs::replicate_estimates(Backend{}, g, L_nan.data(), T, W.data(), B, a.data(), rs::ArgmaxSearch::Full);
    rs::replicate_estimates(Backend{}, g, L_nan.data(), T, W.data(), B, b.data(), rs::ArgmaxSearch::Bounded);
    bool same = true;
    for (std::uint32_t i = 0; i < B; ++i) same = same && same_replicate(a[i], b[i]);
    VCAL_CHECK(same);
    W[5] = -1.0;
    VCAL_CHECK(!rs::detail::argmax_bounded(Backend{}, g, L.data(), T, W.data(), B, best.data()));
}
