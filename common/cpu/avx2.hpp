// SPDX-License-Identifier: Apache-2.0
//
// Whether this CPU and its OS support AVX2 (CPUID leaf 7 EBX bit 5, with OSXSAVE and the YMM state in XCR0), for
// the kernels selected at run time: the DGP's Bernoulli count (P-17, D-194) and the AR(1) filter's transition
// (P-18, D-199). One definition shared by both (D-199). Include it only from files compiled for the baseline
// target, never from an AVX2-compiled one (dgp/philox_lanes.hpp explains the hazard).
#pragma once

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_AMD64))
#include <intrin.h>
#endif

namespace vcal::cpu {

inline bool has_avx2() {
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_AMD64))
    int r[4];
    __cpuid(r, 0);
    if (r[0] < 7) return false;
    __cpuid(r, 1);
    const bool osxsave = (r[2] & (1 << 27)) != 0;
    const bool avx = (r[2] & (1 << 28)) != 0;
    if (!osxsave || !avx || (_xgetbv(0) & 6u) != 6u) return false;
    __cpuidex(r, 7, 0);
    return (r[1] & (1 << 5)) != 0;
#elif defined(__x86_64__)
    // libgcc and compiler-rt check the OS's XCR0 for AVX features as well.
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") != 0;
#else
    return false;
#endif
}

}  // namespace vcal::cpu
