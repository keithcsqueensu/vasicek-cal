// SPDX-License-Identifier: Apache-2.0
//
// The grade scenario matrix (D-151; studies/README.md, "The grade scenario matrix"; D-414): nine scenarios,
// three grade scales x T in {20, 40, 100}, id = i_scale * 3 + i_T, scales in the order L, M, H.
//   - The master scale: 17 grades, PD_k = 0.03% x 1.5^(k-1), k = 1 .. 17.
//   - Obligors: 10,000 per period, the same in every period; grade k's share proportional to
//     exp(-(k - c)^2 / 32), c the scale's centre; n_k = 10,000 x share rounded down, the remainder one each
//     to the largest fractional parts (ties to the lower k).
//   - Replicate r of scenario id: dgp::simulate_panel_segments on seed 0x475241444553434E ("GRADESCN"),
//     scenario id, replicate r, grade k as segment k - 1 (grade 1 on the base stream), with the scale's rho.
// The shares use dgp::det_exp, so the counts are the same on every platform; unit_multigrade_fit checks them
// against the counts listed in studies/README.md.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "dgp/dgp.hpp"

namespace vcal::recovery::grades {

inline constexpr std::uint64_t kSeed = 0x475241444553434Eull;  // "GRADESCN"
inline constexpr int kGrades = 17;
inline constexpr int kScenarios = 9;
inline constexpr std::int64_t kObligors = 10000;

struct Scale {
    const char* name;
    int centre;  // c, a grade number 1 .. 17
    double rho;
};
inline constexpr Scale kScales[3] = {{"L", 5, 0.24}, {"M", 8, 0.12}, {"H", 11, 0.02}};
inline constexpr std::int64_t kPeriods[3] = {20, 40, 100};

// PD_k for k = 1 .. 17: 1.5^(k-1) is exact in binary64 (3^16 < 2^53), so the product is correctly rounded.
inline double master_pd(int k) {
    double f = 1.0;
    for (int j = 1; j < k; ++j) f *= 1.5;
    return 0.0003 * f;
}

// n_k for k = 1 .. 17 (index k - 1).
inline std::vector<std::int64_t> obligors(const Scale& scale) {
    std::vector<double> w(kGrades);
    double sum = 0.0;
    for (int k = 1; k <= kGrades; ++k) {
        const double x = static_cast<double>(k - scale.centre);
        w[static_cast<std::size_t>(k - 1)] = dgp::det_exp(-(x * x) / 32.0);
        sum += w[static_cast<std::size_t>(k - 1)];
    }
    std::vector<std::int64_t> n(kGrades);
    std::vector<double> frac(kGrades);
    std::int64_t given = 0;
    for (std::size_t k = 0; k < static_cast<std::size_t>(kGrades); ++k) {
        const double exact = static_cast<double>(kObligors) * w[k] / sum;
        const double whole = static_cast<double>(static_cast<std::int64_t>(exact));
        n[k] = static_cast<std::int64_t>(whole);
        frac[k] = exact - whole;
        given += n[k];
    }
    for (std::int64_t left = kObligors - given; left > 0; --left) {
        std::size_t best = 0;
        for (std::size_t k = 1; k < static_cast<std::size_t>(kGrades); ++k) {
            if (frac[k] > frac[best]) best = k;  // strict: ties go to the lower k
        }
        ++n[best];
        frac[best] = -1.0;
    }
    return n;
}

struct Scenario {
    int id;
    const Scale* scale;
    std::int64_t periods;
};

inline Scenario scenario(int id) {
    return {id, &kScales[id / 3], kPeriods[id % 3]};
}

// Replicate r of scenario id as a periods x 17 panel (row-major, grade k at index k - 1). False if the DGP
// refuses.
inline bool panel(int id, std::uint32_t replicate, std::vector<std::int64_t>& n, std::vector<std::int64_t>& d,
                  std::vector<double>* z = nullptr) {
    const Scenario sc = scenario(id);
    const auto per_period = obligors(*sc.scale);
    const auto T = static_cast<std::size_t>(sc.periods);
    n.resize(T * kGrades);
    d.assign(n.size(), 0);
    for (std::size_t t = 0; t < T; ++t) {
        for (std::size_t k = 0; k < static_cast<std::size_t>(kGrades); ++k) n[t * kGrades + k] = per_period[k];
    }
    double pd[kGrades];
    for (int k = 1; k <= kGrades; ++k) pd[k - 1] = master_pd(k);
    if (z != nullptr) z->assign(T, 0.0);
    const dgp::SegmentPanelSpec spec{kSeed, static_cast<std::uint32_t>(id), replicate, kGrades, pd,
                                     sc.scale->rho, n.data(), sc.periods};
    return dgp::simulate_panel_segments(spec, d.data(), z != nullptr ? z->data() : nullptr) == dgp::Status::Ok;
}

}  // namespace vcal::recovery::grades
