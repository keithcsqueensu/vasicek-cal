// SPDX-License-Identifier: Apache-2.0
//
// Panel simulation (see dgp/dgp.hpp and docs/methodology/dgp.md). Mirrors Mirror.factor and
// Mirror.panel in tools/gen_dgp_tables.py operation for operation. Compiled with strict FP
// semantics and FP contraction off (D-056).
#include <cmath>
#include <cstdint>

#include "dgp/dgp.hpp"
#include "dgp/philox.hpp"

namespace vcal::dgp {
namespace {

PhiloxKey key_of(std::uint64_t seed) {
    return {static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)};
}

}  // namespace

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
    if (!fp_environment_ok()) return Status::FpEnvironment;
    if (!(spec.pd > 0.0 && spec.pd < 1.0 && spec.rho > 0.0 && spec.rho < 1.0)) return Status::InvalidSpec;
    if (!(spec.periods >= 1 && spec.periods <= (std::int64_t{1} << 32)) || spec.obligors == nullptr ||
        defaults == nullptr) {
        return Status::InvalidSpec;
    }
    for (std::int64_t t = 0; t < spec.periods; ++t) {
        if (!(spec.obligors[t] >= 0 && spec.obligors[t] <= kMaxObligors)) return Status::InvalidSpec;
    }

    const PhiloxKey key = key_of(spec.seed);
    const double c = det_probit(spec.pd);
    const double sr = std::sqrt(spec.rho);
    const double s1 = std::sqrt(1.0 - spec.rho);
    for (std::int64_t t = 0; t < spec.periods; ++t) {
        const auto period = static_cast<std::uint32_t>(t);
        double z = 0.0;
        const Status st = factor_draw(spec.seed, spec.scenario, spec.replicate, period, kFactorBlocks, &z);
        if (st != Status::Ok) return st;
        const double p = det_ncdf((c - sr * z) / s1);
        std::int64_t d = 0;
        PhiloxCounter w{};
        for (std::int64_t i = 0; i < spec.obligors[t]; ++i) {
            if (i % 2 == 0) {
                const auto block = static_cast<std::uint32_t>(kFactorBlocks + static_cast<std::uint64_t>(i / 2));
                w = philox4x32({spec.scenario, spec.replicate, period, block}, key);
            }
            const double u = i % 2 == 0 ? uniform_from_words(w[0], w[1]) : uniform_from_words(w[2], w[3]);
            if (u < p) ++d;
        }
        defaults[t] = d;
        if (factors != nullptr) factors[t] = z;
    }
    return Status::Ok;
}

}  // namespace vcal::dgp
