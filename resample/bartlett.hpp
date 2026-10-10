// SPDX-License-Identifier: Apache-2.0
#pragma once

// The Bartlett-corrected profile-likelihood interval (S-4b; a native option, D-182, D-184).
//
// The profile interval's threshold c/2 (c = chi2_1(0.95)) assumes the likelihood-ratio statistic
// W_a = 2 (l_max - P_a(theta_a)) is chi2_1 at the true theta, so E[W_a] = 1. At short histories it
// is larger (S-4: about 1.1 at T = 20). The Bartlett correction rescales the threshold by an estimate
// k_a of E[W_a]: the interval is { theta_a : P_a(theta_a) >= l_max - k_a c / 2 }.
//
// k_a is estimated by parametric bootstrap at the estimate theta-hat:
//   1. panel b = 0 .. panels-1: defaults d*_t simulated from the one-factor DGP (docs/methodology/dgp.md)
//      at (PD-hat, rho-hat) with the panel's own n_t, seed `seed`, scenario `stream`, replicate
//      `first_panel + b`;
//   2. each panel fitted as the original (calibrate, the same grid and rules), and
//      W*_{a,b} = 2 (l*_max - P*_a(theta-hat_a)): the profile statistic of panel b at the value that
//      generated it;
//   3. k_a = the mean of the finite W*_{a,b}, in panel order.
// k_a falls back to 1 (the plain profile interval) when theta-hat_a is on the grid box's bound, where the
// bootstrap distribution is truncated, or when more than 10% of the W*_{a,b} are not finite.
//
// With seed 0x53313050424F4F54 ("S10PBOOT"), stream = the recovery scenario, first_panel = 1000 r and
// 199 panels this is exactly study S-4's S-4b for replicate r (tests: bartlett_reproduces_the_s4_run).
// S-4 found it the most accurate interval tested at T = 20 (rho 0.947, PD 0.945 against Jeffreys'
// 0.943, 0.942), at about 4-5% more width than the profile interval; the preferred interval stays
// Jeffreys equal-tailed (D-182).
//
// Needs the binomial objective (Obs {n, d}) with analytic derivatives (profile_lr_statistic), and
// vcal::dgp at link time.

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "core/grid.hpp"
#include "dgp/dgp.hpp"
#include "engine/calibrate.hpp"
#include "engine/profile.hpp"
#include "engine/surface_cache.hpp"

namespace vcal::resample {

inline constexpr std::uint32_t kBartlettPanels = 199;  // S-4b's count (D-176)
inline constexpr double kBartlettMaxLeftOut = 0.10;    // more non-finite W* than this: k = 1

// Per-parameter flags on the factor.
enum : std::uint32_t {
    kBartlettOnBound = 1u << 0,      // theta-hat_a on the grid bound: k_a = 1
    kBartlettTooFewFinite = 1u << 1, // more than 10% of W*_a not finite: k_a = 1
    kBartlettNotComputed = 1u << 2,  // no usable estimate (flat or numeric, or outside (0, 1)): no bootstrap
};

struct BartlettSpec {
    std::uint64_t seed;
    std::uint32_t stream;       // the DGP's scenario field
    std::uint32_t first_panel;  // panel b is the DGP's replicate first_panel + b
    std::uint32_t panels;       // >= 1; kBartlettPanels in S-4
};

struct BartlettIntervals2 {
    engine::ProfileIntervals2 intervals;  // per parameter at its own threshold k_a c / 2; residual_max over
                                          // both solves, evaluations summed over them (the bootstrap's not)
    double factor[2];                     // k_a (1 on a fallback; NaN if not computed)
    std::uint32_t factor_flags[2];        // kBartlett*
    std::int64_t finite[2];               // finite W*_a
    std::int64_t panels;                  // panels simulated and fitted
    std::int64_t panels_quad_flagged;     // bootstrap fits whose quadrature check flagged a period
};

enum class BartlettStatus : std::int32_t {
    Ok = 0,
    InvalidSpec = 1,    // panels == 0, an n_t above the DGP's limit, or first_panel + panels - 1 >= 2^32
    FpEnvironment = 2,  // flush-to-zero or denormals-are-zero is on (the DGP refuses)
    FactorStreamExhausted = 3,  // the DGP's factor draw ran out of its Philox block range
    KernelUnavailable = 4,      // VCAL_DGP_KERNEL names a Bernoulli kernel this CPU or build lacks (D-194)
};

// est and L are the panel's own fit (calibrate's outputs on the same grid and rules). The bootstrap panels'
// surfaces come from a row cache (D-167; identical results): the panels share the n_t, so their rows repeat.
// A caller fitting many panels on one configuration may pass its own cache; otherwise one is made per call.
template <class Backend, class Objective, class Integrator>
BartlettStatus bartlett_profile_intervals(const Backend& backend, const Objective& objective, const Integrator& primary,
                                          const Integrator& check, const typename Objective::Obs* obs,
                                          std::int64_t periods, const Grid<2>& grid, const std::vector<double>& L,
                                          const engine::Estimate2& est, const BartlettSpec& spec,
                                          BartlettIntervals2& out, engine::SurfaceRowCache* cache = nullptr) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    out = BartlettIntervals2{};
    out.factor[0] = out.factor[1] = nan;
    if (spec.panels == 0 || std::uint64_t{spec.first_panel} + spec.panels - 1 > 0xFFFFFFFFull) {
        return BartlettStatus::InvalidSpec;
    }
    std::vector<std::int64_t> n(static_cast<std::size_t>(periods)), d(n.size());
    for (std::int64_t t = 0; t < periods; ++t) {
        n[static_cast<std::size_t>(t)] = obs[t].n;
        if (obs[t].n > dgp::kMaxObligors) return BartlettStatus::InvalidSpec;
    }
    const double v_hat[2] = {est.value[0], est.value[1]};
    const bool usable = !(est.flags & (engine::kFlagFlatSurface | engine::kFlagNumeric)) && v_hat[0] > 0.0 &&
                        v_hat[0] < 1.0 && v_hat[1] > 0.0 && v_hat[1] < 1.0;
    if (!usable) {
        out.intervals = engine::profile_intervals(objective, primary, obs, periods, grid, L, est);
        out.factor_flags[0] = out.factor_flags[1] = kBartlettNotComputed;
        return BartlettStatus::Ok;
    }

    engine::SurfaceRowCache local;
    engine::SurfaceRowCache& rows = cache != nullptr ? *cache : local;
    double sum[2] = {0.0, 0.0};
    std::vector<typename Objective::Obs> obs_b(n.size());
    std::vector<double> Lb;
    for (std::uint32_t b = 0; b < spec.panels; ++b) {
        const dgp::PanelSpec ps{spec.seed, spec.stream, spec.first_panel + b, v_hat[0], v_hat[1], n.data(), periods};
        switch (dgp::simulate_panel(ps, d.data(), nullptr)) {
            case dgp::Status::Ok: break;
            case dgp::Status::FpEnvironment: return BartlettStatus::FpEnvironment;
            case dgp::Status::FactorStreamExhausted: return BartlettStatus::FactorStreamExhausted;
            case dgp::Status::KernelUnavailable: return BartlettStatus::KernelUnavailable;
            default: return BartlettStatus::InvalidSpec;
        }
        for (std::size_t t = 0; t < n.size(); ++t) obs_b[t] = {n[t], d[t]};
        engine::Estimate2 eb{};
        ++out.panels;
        if (engine::calibrate_cached(backend, rows, objective, primary, check, obs_b.data(), periods, grid, Lb, eb) !=
            engine::Status::Ok) {
            continue;  // no finite surface: W* not finite for this panel
        }
        if (eb.quad_check_flagged > 0) ++out.panels_quad_flagged;
        const auto lr = engine::profile_lr_statistic(objective, primary, obs_b.data(), periods, grid, Lb, eb, v_hat);
        for (int a = 0; a < 2; ++a) {
            if (std::isfinite(lr.w[a])) {
                sum[a] += lr.w[a];
                ++out.finite[a];
            }
        }
    }

    double threshold[2];
    for (int a = 0; a < 2; ++a) {
        const bool on_bound = !(v_hat[a] > grid.axis[a].lo && v_hat[a] < grid.axis[a].hi);
        const bool too_few = static_cast<double>(spec.panels) - static_cast<double>(out.finite[a]) >
                             kBartlettMaxLeftOut * spec.panels;
        out.factor_flags[a] = (on_bound ? kBartlettOnBound : 0u) | (too_few || out.finite[a] == 0 ? kBartlettTooFewFinite : 0u);
        out.factor[a] = out.factor_flags[a] != 0u ? 1.0 : sum[a] / static_cast<double>(out.finite[a]);
        threshold[a] = engine::kProfileThreshold95 * out.factor[a];
    }
    for (int a = 0; a < 2; ++a) {
        const auto p = engine::profile_intervals(objective, primary, obs, periods, grid, L, est, threshold[a]);
        if (a == 0) {
            out.intervals = p;
        } else {
            out.intervals.residual_max = std::fmax(out.intervals.residual_max, p.residual_max);
            out.intervals.evaluations += p.evaluations;
        }
        out.intervals.lo[a] = p.lo[a];
        out.intervals.hi[a] = p.hi[a];
        out.intervals.flags[a] = p.flags[a];
    }
    return BartlettStatus::Ok;
}

}  // namespace vcal::resample
