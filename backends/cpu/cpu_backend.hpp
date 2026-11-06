// SPDX-License-Identifier: Apache-2.0
//
// CPU backend (D-009, D-010, D-028): OpenMP `parallel for` over independent items, with
// no reduction clauses, so it is portable to MSVC's OpenMP 2.0 and deterministic by
// construction. Without OpenMP it runs serially and gives identical results.
//
// n_threads == 0 uses the OpenMP default.
//
// Exceptions (D-218): an exception cannot leave an OpenMP parallel region (it would call
// std::terminate), so each item runs inside a try block. The first exception caught is kept, the
// items not yet started are skipped, and it is rethrown after the region, so a callable that
// throws (std::bad_alloc, say) reaches the caller exactly as it would in a serial build. Which
// exception is "first" can depend on scheduling when several items throw; results are only
// returned when none does, so determinism of results is unaffected.
#pragma once

#include <atomic>
#include <cstdint>
#include <exception>

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
        std::exception_ptr error;
        std::atomic<bool> failed{false};
#pragma omp parallel for schedule(static) num_threads(threads)
        for (std::int64_t i = 0; i < count; ++i) {
            if (failed.load(std::memory_order_relaxed)) continue;
            try {
                f(i);
            } catch (...) {
#pragma omp critical(vcal_cpu_backend_error)
                {
                    if (!error) error = std::current_exception();
                }
                failed.store(true, std::memory_order_relaxed);
            }
        }
        if (error) std::rethrow_exception(error);
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
