// SPDX-License-Identifier: Apache-2.0
//
// Panel simulation (see dgp/dgp.hpp and docs/methodology/dgp.md). Mirrors Mirror.factor, Mirror.panel and
// Mirror.panel_segments in tools/gen_dgp_tables.py operation for operation, the Bernoulli count through
// count_defaults (dgp/bernoulli.cpp; its vector kernels equal the mirror's loop bit for bit, P-17).
// Compiled with strict FP semantics and FP contraction off (D-056).
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "dgp/dgp.hpp"
#include "dgp/philox.hpp"

namespace vcal::dgp {
namespace {

PhiloxKey key_of(std::uint64_t seed) {
    return {static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)};
}

// The base DGP (ar1 false) and the AR(1)-factor variant (ar1 true), which differ only in Z_t.
Status simulate(const PanelSpec& spec, bool ar1, double phi, std::int64_t* defaults, double* factors) {
    if (!fp_environment_ok()) return Status::FpEnvironment;
    if (!(spec.pd > 0.0 && spec.pd < 1.0 && spec.rho > 0.0 && spec.rho < 1.0)) return Status::InvalidSpec;
    if (ar1 && !(phi > -1.0 && phi < 1.0)) return Status::InvalidSpec;
    if (!(spec.periods >= 1 && spec.periods <= (std::int64_t{1} << 32)) || spec.obligors == nullptr ||
        defaults == nullptr) {
        return Status::InvalidSpec;
    }
    for (std::int64_t t = 0; t < spec.periods; ++t) {
        if (!(spec.obligors[t] >= 0 && spec.obligors[t] <= kMaxObligors)) return Status::InvalidSpec;
    }

    const BernoulliKernel kernel = bernoulli_kernel_selected();
    if (!bernoulli_kernel_available(kernel)) return Status::KernelUnavailable;
    const double c = det_probit(spec.pd);
    const double sr = std::sqrt(spec.rho);
    const double s1 = std::sqrt(1.0 - spec.rho);
    const double se = ar1 ? std::sqrt(1.0 - phi * phi) : 1.0;  // the innovation's scale
    double z_prev = 0.0;
    for (std::int64_t t = 0; t < spec.periods; ++t) {
        const auto period = static_cast<std::uint32_t>(t);
        double e = 0.0;
        const Status st = factor_draw(spec.seed, spec.scenario, spec.replicate, period, kFactorBlocks, &e);
        if (st != Status::Ok) return st;
        const double z = ar1 && t > 0 ? phi * z_prev + se * e : e;
        z_prev = z;
        const double p = det_ncdf((c - sr * z) / s1);
        defaults[t] = count_defaults(spec.seed, spec.scenario, spec.replicate, period, spec.obligors[t], p, kernel);
        if (factors != nullptr) factors[t] = z;
    }
    return Status::Ok;
}

}  // namespace

Status simulate_panel_segments(const SegmentPanelSpec& spec, std::int64_t* defaults, double* factors) {
    if (!fp_environment_ok()) return Status::FpEnvironment;
    if (!(spec.segments >= 1 && spec.segments <= kMaxSegments) || spec.pd == nullptr) return Status::InvalidSpec;
    if (!(spec.rho > 0.0 && spec.rho < 1.0)) return Status::InvalidSpec;
    const auto K = static_cast<std::size_t>(spec.segments);
    for (std::size_t s = 0; s < K; ++s) {
        if (!(spec.pd[s] > 0.0 && spec.pd[s] < 1.0)) return Status::InvalidSpec;
    }
    if (!(spec.periods >= 1 && spec.periods <= (std::int64_t{1} << 32)) || spec.obligors == nullptr ||
        defaults == nullptr) {
        return Status::InvalidSpec;
    }
    const auto cells = static_cast<std::size_t>(spec.periods) * K;
    for (std::size_t i = 0; i < cells; ++i) {
        if (!(spec.obligors[i] >= 0 && spec.obligors[i] <= kMaxObligors)) return Status::InvalidSpec;
    }

    const BernoulliKernel kernel = bernoulli_kernel_selected();
    if (!bernoulli_kernel_available(kernel)) return Status::KernelUnavailable;
    double c[kMaxSegments];
    std::uint64_t seeds[kMaxSegments];
    for (std::size_t s = 0; s < K; ++s) {
        c[s] = det_probit(spec.pd[s]);
        seeds[s] = segment_seed(spec.seed, static_cast<int>(s));
    }
    const double sr = std::sqrt(spec.rho);
    const double s1 = std::sqrt(1.0 - spec.rho);
    for (std::int64_t t = 0; t < spec.periods; ++t) {
        const auto period = static_cast<std::uint32_t>(t);
        double z = 0.0;
        const Status st = factor_draw(spec.seed, spec.scenario, spec.replicate, period, kFactorBlocks, &z);
        if (st != Status::Ok) return st;
        const std::size_t row = static_cast<std::size_t>(t) * K;
        for (std::size_t s = 0; s < K; ++s) {
            const double p = det_ncdf((c[s] - sr * z) / s1);
            defaults[row + s] =
                count_defaults(seeds[s], spec.scenario, spec.replicate, period, spec.obligors[row + s], p, kernel);
        }
        if (factors != nullptr) factors[t] = z;
    }
    return Status::Ok;
}

Status factor_draw(std::uint64_t seed, std::uint32_t scenario, std::uint32_t replicate, std::uint32_t period,
                   std::uint32_t max_blocks, double* z) {
    const PhiloxKey key = key_of(seed);
    for (std::uint32_t block = 0; block < max_blocks; ++block) {
        const PhiloxCounter w = philox4x32({scenario, replicate, period, block}, key);
        const double v1 = 2.0 * uniform_from_words(w[0], w[1]) - 1.0;
        const double v2 = 2.0 * uniform_from_words(w[2], w[3]) - 1.0;
        const double s = v1 * v1 + v2 * v2;
        if (s > 0.0 && s < 1.0) {
            *z = v1 * std::sqrt((-2.0 * det_log(s)) / s);
            return Status::Ok;
        }
    }
    return Status::FactorStreamExhausted;
}

Status simulate_panel(const PanelSpec& spec, std::int64_t* defaults, double* factors) {
    return simulate(spec, false, 0.0, defaults, factors);
}

Status simulate_panel_ar1(const PanelSpec& spec, double phi, std::int64_t* defaults, double* factors) {
    return simulate(spec, true, phi, defaults, factors);
}

}  // namespace vcal::dgp
