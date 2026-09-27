// SPDX-License-Identifier: Apache-2.0
//
// CPU backend (D-009, D-010, D-028): OpenMP `parallel for` over independent items, with
// no reduction clauses, so it is portable to MSVC's OpenMP 2.0 and deterministic by
// construction. Without OpenMP it runs serially and gives identical results.
//
// n_threads == 0 uses the OpenMP default. The callable must not throw: an exception cannot
// leave an OpenMP parallel region.
#pragma once

#include <cstdint>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace vcal::backends {

struct CpuBackend {
    int n_threads = 0;

    template <class F>
    void parallel_for(std::int64_t count, const F& f) const {
#if defined(_OPENMP)
        const int threads = n_threads > 0 ? n_threads : omp_get_max_threads();
#pragma omp parallel for schedule(static) num_threads(threads)
        for (std::int64_t i = 0; i < count; ++i) f(i);
#else
        for (std::int64_t i = 0; i < count; ++i) f(i);
#endif
    }

    static bool has_openmp() {
#if defined(_OPENMP)
        return true;
#else
        return false;
#endif
    }
};

}  // namespace vcal::backends
