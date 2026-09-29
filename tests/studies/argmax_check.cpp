// SPDX-License-Identifier: Apache-2.0
//
// The bounded argmax against the full grid (D-172), on recovery panels: for every replicate of
// every scenario, every row of its iid bootstrap (B = 999) and of its jackknife, the argmax state
// (index, value, NaN count) from resample::detail::argmax_bounded must equal the full-grid
// reduction's bit for bit, over the panel's distinct rows as the engine runs them. Reports the rows
// compared, how many argmaxes lie on the grid's edge (the bound cases, where coarse searches fail
// first), and how many calls the guard sent to the full grid. Exits 1 on any difference.
//
//   study_argmax_check [--replicates R] [--scenarios ID,...]
// Defaults: the study subset of D-150, R = 1000 (the verification behind D-172). The slow CTest
// runs a reduced R on the subset, so the equality is re-checked on every CI run.
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "core/reducers/argmax.hpp"
#include "engine/surface_cache.hpp"
#include "resample/weights.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;

struct Counts {
    std::int64_t rows = 0, differing = 0, edge = 0, fallback = 0;
};

bool same_state(const reducers::ArgMax::State& a, const reducers::ArgMax::State& b) {
    return std::memcmp(&a.best, &b.best, sizeof a.best) == 0 && a.k == b.k && a.nan_count == b.nan_count;
}

}  // namespace

int main(int argc, char** argv) {
    std::uint32_t R = rc::kReplicates;
    std::vector<std::uint32_t> ids = {29, 37, 72, 4, 68, 49, 7, 43, 51};  // the study subset (D-150)
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--scenarios") == 0 && i + 1 < argc) {
            ids.clear();
            for (const char* c = argv[++i]; *c != '\0';) {
                char* end = nullptr;
                ids.push_back(static_cast<std::uint32_t>(std::strtoul(c, &end, 10)));
                c = *end == ',' ? end + 1 : end;
            }
        } else {
            std::fprintf(stderr, "usage: study_argmax_check [--replicates R] [--scenarios ID,...]\n");
            return 2;
        }
    }
    const Grid<2> g = rc::grid();
    const std::int64_t K = g.size();
    const auto primary = quadrature::parity_rule();
    const backends::CpuBackend all{}, serial{1};
    Counts total;
    std::printf("scenario  T    n     | rows compared  differing  on the grid edge  guard fallbacks\n");
    for (const std::uint32_t id : ids) {
        const rc::Scenario s = rc::scenario(id);
        const std::int64_t T = s.periods;
        engine::SurfaceRowCache cache;
        std::vector<Counts> per(R);
        all.parallel_for(R, [&](std::int64_t r) {
            Counts& c = per[static_cast<std::size_t>(r)];
            const auto d = rc::panel(s, static_cast<std::uint32_t>(r));
            std::vector<rc::Objective::Obs> obs(d.size());
            for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {s.obligors, d[t]};
            std::vector<double> L(static_cast<std::size_t>(T * K));
            engine::evaluate_surface_cached(serial, cache, rc::Objective{}, primary, obs.data(), T, g, L.data());
            const auto idx = resample::bootstrap_indices(rc::bootstrap_seed(s.id, static_cast<std::uint32_t>(r)),
                                                         resample::Scheme::IidBootstrap, rc::kBootstrapReplicates, T);
            const auto boot = resample::weights_from_indices(idx.data(), rc::kBootstrapReplicates, T, T);
            const auto jack = resample::jackknife_weights(T);
            for (const auto* W : {&boot, &jack}) {
                const std::int64_t B = static_cast<std::int64_t>(W->size()) / T;
                const auto cp = resample::compact_by_observation(obs.data(), T, L.data(), K, W->data(), B);
                std::vector<reducers::ArgMax::State> full(static_cast<std::size_t>(B)), bounded(full.size());
                engine::reduce_weighted(serial, reducers::ArgMax{}, cp.rows.data(), cp.distinct, K, cp.weights.data(), B,
                                        full.data());
                if (!resample::detail::argmax_bounded(serial, g, cp.rows.data(), cp.distinct, cp.weights.data(), B,
                                                      bounded.data())) {
                    ++c.fallback;
                    continue;  // the engine would use the full grid itself
                }
                for (std::int64_t b = 0; b < B; ++b) {
                    const auto& f = full[static_cast<std::size_t>(b)];
                    ++c.rows;
                    if (!same_state(f, bounded[static_cast<std::size_t>(b)])) ++c.differing;
                    std::int32_t i[2];
                    g.unflatten(f.k, i);
                    if (i[0] == 0 || i[0] == g.axis[0].n - 1 || i[1] == 0 || i[1] == g.axis[1].n - 1) ++c.edge;
                }
            }
        });
        Counts sc;
        for (const auto& c : per) {
            sc.rows += c.rows;
            sc.differing += c.differing;
            sc.edge += c.edge;
            sc.fallback += c.fallback;
        }
        std::printf("%-9u %-4lld %-5lld | %13lld  %9lld  %16lld  %15lld\n", id, static_cast<long long>(T),
                    static_cast<long long>(s.obligors), static_cast<long long>(sc.rows),
                    static_cast<long long>(sc.differing), static_cast<long long>(sc.edge),
                    static_cast<long long>(sc.fallback));
        total.rows += sc.rows;
        total.differing += sc.differing;
        total.edge += sc.edge;
        total.fallback += sc.fallback;
    }
    std::printf("total: %lld rows compared, %lld differing, %lld on the grid edge, %lld guard fallbacks\n",
                static_cast<long long>(total.rows), static_cast<long long>(total.differing),
                static_cast<long long>(total.edge), static_cast<long long>(total.fallback));
    return total.differing == 0 ? 0 : 1;
}
