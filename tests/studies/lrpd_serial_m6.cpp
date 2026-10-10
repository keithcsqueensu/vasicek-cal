// SPDX-License-Identifier: Apache-2.0
//
// S-35's M6 arm (studies/lrpd-serial-correlation/PREDICTION.md, the 2026-10-10 addendum, L16-L28; D-195,
// D-196, D-198): the AR(1)-aware binomial MLE (benchmark_models/ar1/) on S-35's panels: the recovery matrix under
// the AR(1)-factor DGP (dgp::simulate_panel_ar1) at each phi, on the recovery seed, scenario ids and replicates,
// the panels study_lrpd_serial fits. Per panel:
//   - the binomial-mixture MLE (calibrate on the recovery grid, with the shared row cache): the search's start at
//     phi = 0 and the reference of the first check before scoring;
//   - the M6 estimate (ar1::fit): PD, rho, phi, SEs, log-likelihood, the phi = 0 scan point, the doubled-rule and
//     placement checks, the scan's local maxima and the gap to the second, flags;
//   - its 3-D profile intervals for PD and rho, and for phi on the study subset (D-150) only;
//   - on the study subset, replicates 0-99, phi = 0.6 only: the parametric-bootstrap comparator, B = 199 panels
//     from the AR(1) DGP at the panel's own estimate, panel b of (s, r) on seed 0xB0071A61, scenario s,
//     replicate 256 r + b, each refitted for the point estimate from its own binomial MLE; percentile intervals
//     (type 7) for PD and rho.
// Nothing is scored here.
//
// Work units (owner, D-198): a fit job is one (scenario, replicate) at every phi; when it owes a comparator, it
// releases B refit tasks, which workers take before new fits, so no replicate is one long serial job. Each refit
// writes its own slot, and a replicate's percentiles and timings are assembled in b order once its last refit
// is in: the rows do not depend on the order or the thread count.
//
// --checkpoint DIR (P-11, D-179): each replicate is written to DIR when it is complete (its fits and, if owed,
// all B refits), so a partial set of refits is never taken for a finished one; a rerun with the same
// configuration skips the replicates on disk (tests/harness/checkpoint.hpp). --stop-after-units N exits after N
// replicates are written in this invocation, as a kill would (the resume test).
//
// --seed S runs the same design on another seed (the validation run, D-197: 0x7A11DA7E); the comparator is
// S-35's only and is not run then, unless --comparator (with --boot-panels B, for smoke tests). --boot-sample K,
// for runtime timing samples only (perf/README.md, D-187): K refits per comparator, their seconds scaled by B/K.
// Its comparator intervals would be wrong as results, so in that mode the tool writes no rows at all, only the
// timing sample (owner). --stop-after-units needs --checkpoint and exits before any output is written.
//
//   study_lrpd_serial_m6 [--phis 0,0.3,0.6] [--replicates R] [--scenarios ID,...] [--seed S] [--out FILE]
//                        [--threads N] [--timing-out FILE] [--boot-sample K] [--checkpoint DIR]
//                        [--stop-after-units N] [--comparator] [--boot-panels B]
//
// Doubles are written to 17 significant digits.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include "benchmark_models/ar1/search.hpp"
#include "generated/vcal_git_info.h"
#include "resample/bootstrap.hpp"
#include "tests/harness/checkpoint.hpp"
#include "tests/harness/dgp_provenance.hpp"
#include "tests/harness/timing_sample.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
namespace e = vcal::engine;
namespace ar1 = vcal::benchmark_models::ar1;
using Clock = std::chrono::steady_clock;

constexpr std::uint64_t kBootSeed = 0xB0071A61ull;          // the comparator's panels (PREDICTION.md)
constexpr std::uint32_t kBootPanels = 199;                  // B (--boot-panels: smoke tests only)
constexpr std::uint32_t kBootReplicates = 100;              // replicates 0-99
constexpr double kBootPhi = 0.6;                            // phi = 0.6 only
constexpr std::uint32_t kSubset[] = {29, 37, 72, 4, 68, 49, 7, 43, 51};  // the study subset (D-150)
const ar1::Box kBox{1e-4, 0.2, 1e-3, kDefaultRhoUpper};    // the recovery box (D-089, D-115)

enum Phase { kSimulate, kBinomial, kFit, kProfile, kProfilePhi, kBootstrap, kCacheFill, kPhases };
const char* const kPhaseNames[kPhases] = {"simulate", "binomial MLE", "M6 fit", "M6 profile PD rho",
                                          "M6 profile phi", "M6 bootstrap", "cache fill"};
using Seconds = std::array<double, kPhases>;

struct Timer {
    double* acc;
    Clock::time_point t0 = Clock::now();
    ~Timer() { *acc += std::chrono::duration<double>(Clock::now() - t0).count(); }
};

bool in_subset(std::uint32_t id) { return std::find(std::begin(kSubset), std::end(kSubset), id) != std::end(kSubset); }

struct Row {
    double binom_pd, binom_rho, binom_loglik;
    std::uint32_t binom_flags;
    ar1::Estimate3 est;
    std::uint32_t fitted;
    ar1::ProfileIntervals3 prof;
    double boot_lo[2], boot_hi[2];
    std::uint32_t boot_used, boot_excluded;  // comparator panels refitted with finite estimates / not
};

double nan() { return std::nan(""); }

std::vector<ar1::Obs> as_obs(const std::vector<std::int64_t>& d, std::int64_t n) {
    std::vector<ar1::Obs> obs(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {n, d[t]};
    return obs;
}

// The binomial-mixture MLE on the recovery grid, its surface from the shared row cache (fills timed apart).
e::Estimate2 binomial(const rc::Scenario& s, const std::vector<std::int64_t>& d, e::SurfaceRowCache& cache, double* sec,
                      Phase phase) {
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    std::vector<rc::Objective::Obs> obs(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {s.obligors, d[t]};
    std::vector<double> L;
    e::Estimate2 est{};
    const vcal::test::FillMeter fills;
    {
        Timer t{&sec[phase]};
        if (e::calibrate_cached(backends::CpuBackend{1}, cache, rc::Objective{}, primary, check, obs.data(), s.periods,
                                rc::grid(), L, est) != e::Status::Ok) {
            std::abort();
        }
    }
    sec[phase] -= fills.seconds();  // a fill is paid once per run, not per replicate (D-193)
    sec[kCacheFill] += fills.seconds();
    return est;
}

// One panel's binomial MLE, M6 estimate and profile intervals (the comparator comes from the refit tasks).
Row fit_panel(const rc::Scenario& s, const std::vector<std::int64_t>& d, e::SurfaceRowCache& cache, double* sec) {
    Row row{};
    row.boot_lo[0] = row.boot_lo[1] = row.boot_hi[0] = row.boot_hi[1] = nan();
    for (int a = 0; a < 3; ++a) row.prof.lo[a] = row.prof.hi[a] = nan();
    for (auto& f : row.prof.flags) f = e::kIntervalNotComputed;
    const auto T = s.periods;
    {
        const e::Estimate2 est = binomial(s, d, cache, sec, kBinomial);
        row.binom_pd = est.value[0];
        row.binom_rho = est.value[1];
        row.binom_loglik = est.loglik;
        row.binom_flags = est.flags;
    }
    const auto obs = as_obs(d, s.obligors);
    {
        Timer t{&sec[kFit]};
        row.fitted = ar1::fit(obs.data(), T, kBox, row.binom_pd, row.binom_rho, row.est) ? 1u : 0u;
    }
    if (!row.fitted) return row;
    {
        Timer t{&sec[kProfile]};
        const bool which[3] = {true, true, false};
        ar1::profile_intervals(obs.data(), T, kBox, row.est, which, row.prof);
    }
    if (in_subset(s.id)) {
        Timer t{&sec[kProfilePhi]};
        const bool which[3] = {false, false, true};
        ar1::ProfileIntervals3 p{};
        ar1::profile_intervals(obs.data(), T, kBox, row.est, which, p);
        row.prof.lo[2] = p.lo[2];
        row.prof.hi[2] = p.hi[2];
        row.prof.flags[2] = p.flags[2];
        row.prof.residual_max = std::fmax(row.prof.residual_max, p.residual_max);
    }
    return row;
}

// One comparator refit: panel b from the AR(1) DGP at the estimate, refitted exactly as the panel itself was (the
// scan starts at the bootstrap panel's own binomial MLE). Point estimate only (PREDICTION.md).
struct Refit {
    double pd, rho;
    std::uint32_t ok;
    Seconds sec;
};

Refit refit(const rc::Scenario& s, std::uint32_t r, std::uint32_t b, const ar1::Estimate3& est,
            e::SurfaceRowCache& cache) {
    Refit out{nan(), nan(), 0u, {}};
    const auto T = s.periods;
    std::vector<std::int64_t> n(static_cast<std::size_t>(T), s.obligors), db(n.size());
    {
        Timer t{&out.sec[kBootstrap]};
        const dgp::PanelSpec spec{kBootSeed, s.id, 256u * r + b, est.value[0], est.value[1], n.data(), T};
        if (dgp::simulate_panel_ar1(spec, est.value[2], db.data(), nullptr) != dgp::Status::Ok) std::abort();
    }
    const e::Estimate2 bin = binomial(s, db, cache, out.sec.data(), kBootstrap);
    const auto ob = as_obs(db, s.obligors);
    ar1::Estimate3 eb{};
    bool fitted;
    {
        Timer t{&out.sec[kBootstrap]};
        fitted = ar1::fit(ob.data(), T, kBox, bin.value[0], bin.value[1], eb);
    }
    if (fitted && std::isfinite(eb.value[0]) && std::isfinite(eb.value[1])) {
        out.pd = eb.value[0];
        out.rho = eb.value[1];
        out.ok = 1u;
    }
    return out;
}

// A replicate's checkpoint unit: a header and its P rows, byte for byte as held in memory.
struct UnitHeader {
    char magic[8];
    std::uint32_t scenario, replicate, rows;
    std::uint32_t row_size;
};

std::vector<std::uint8_t> pack(std::uint32_t id, std::uint32_t r, const Row* rows, std::uint32_t P) {
    UnitHeader h{{'M', '6', 'C', 'K', 'P', 'T', '0', '1'}, id, r, P, static_cast<std::uint32_t>(sizeof(Row))};
    std::vector<std::uint8_t> out(sizeof h + sizeof(Row) * P);
    std::memcpy(out.data(), &h, sizeof h);
    std::memcpy(out.data() + sizeof h, rows, sizeof(Row) * P);
    return out;
}

void unpack(const std::vector<std::uint8_t>& bytes, std::uint32_t id, std::uint32_t r, Row* rows, std::uint32_t P) {
    UnitHeader h{};
    if (bytes.size() != sizeof h + sizeof(Row) * P) throw std::runtime_error("checkpoint unit has the wrong size");
    std::memcpy(&h, bytes.data(), sizeof h);
    if (std::memcmp(h.magic, "M6CKPT01", 8) != 0 || h.scenario != id || h.replicate != r || h.rows != P ||
        h.row_size != sizeof(Row)) {
        throw std::runtime_error("checkpoint unit does not match its replicate");
    }
    std::memcpy(rows, bytes.data() + sizeof h, sizeof(Row) * P);
}

std::vector<double> parse_list(const char* p) {
    std::vector<double> out;
    while (*p) {
        char* end = nullptr;
        out.push_back(std::strtod(p, &end));
        if (end == p) std::abort();
        p = *end == ',' ? end + 1 : end;
    }
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    if (!vcal::test::announce_dgp_kernel("study_lrpd_serial_m6")) return 2;
    std::vector<double> phis = {0.0, 0.3, 0.6};
    std::uint32_t R = rc::kReplicates;
    std::vector<std::uint32_t> ids;
    std::uint64_t seed = rc::kSeed;
    std::string out = "lrpd_serial_m6.csv", timing_out, checkpoint_dir;
    unsigned threads = std::thread::hardware_concurrency();
    bool force_comparator = false;
    std::uint32_t boot_panels = kBootPanels, boot_sample = 0;
    std::size_t stop_after = 0;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--phis") == 0 && i + 1 < argc) {
            phis = parse_list(argv[++i]);
        } else if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--scenarios") == 0 && i + 1 < argc) {
            for (const double v : parse_list(argv[++i])) ids.push_back(static_cast<std::uint32_t>(v));
        } else if (std::strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            seed = std::strtoull(argv[++i], nullptr, 0);
        } else if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out = argv[++i];
        } else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--timing-out") == 0 && i + 1 < argc) {
            timing_out = argv[++i];
        } else if (std::strcmp(argv[i], "--boot-sample") == 0 && i + 1 < argc) {
            boot_sample = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--checkpoint") == 0 && i + 1 < argc) {
            checkpoint_dir = argv[++i];
        } else if (std::strcmp(argv[i], "--stop-after-units") == 0 && i + 1 < argc) {
            stop_after = static_cast<std::size_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--comparator") == 0) {
            force_comparator = true;
        } else if (std::strcmp(argv[i], "--boot-panels") == 0 && i + 1 < argc) {
            boot_panels = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else {
            std::fprintf(stderr, "usage: study_lrpd_serial_m6 [--phis 0,0.3,0.6] [--replicates R] [--scenarios ID,...] "
                                 "[--seed S] [--out FILE] [--threads N] [--timing-out FILE] [--boot-sample K] "
                                 "[--checkpoint DIR] [--stop-after-units N] [--comparator] [--boot-panels B]\n");
            return 2;
        }
    }
    if (ids.empty()) {
        for (std::uint32_t id = 0; id < rc::kScenarios; ++id) ids.push_back(id);
    }
    if (threads == 0) threads = 1;
    if (boot_sample > 0 && (timing_out.empty() || boot_sample > boot_panels)) {
        std::fprintf(stderr, "--boot-sample K is for timing samples: it needs --timing-out and K <= B\n");
        return 2;
    }
    if (stop_after > 0 && checkpoint_dir.empty()) {
        std::fprintf(stderr, "--stop-after-units N needs --checkpoint DIR (it is the resume test's kill)\n");
        return 2;
    }
    const bool boot = seed == rc::kSeed || force_comparator;  // the comparator belongs to S-35's run
    const std::uint32_t refits_per_job = boot_sample > 0 ? boot_sample : boot_panels;
    const auto P = static_cast<std::uint32_t>(phis.size());
    const std::size_t S = ids.size();
    std::uint32_t pb = P;  // the index of phi = 0.6 among phis, if any
    for (std::uint32_t p = 0; p < P; ++p) {
        if (phis[p] == kBootPhi) pb = p;
    }

    // The run's configuration, for the checkpoint (P-11): everything that determines the rows.
    std::string exe = vcal::test::file_sha256(argv[0]);
    if (exe.empty()) exe = vcal::test::file_sha256(std::string(argv[0]) + ".exe");
    std::ostringstream cfg;
    cfg << "tool=study_lrpd_serial_m6\nreplicates=" << R << "\nscenarios=";
    for (std::size_t k = 0; k < S; ++k) cfg << (k ? "," : "") << ids[k];
    cfg << "\nphis=";
    for (std::uint32_t p = 0; p < P; ++p) cfg << (p ? "," : "") << phis[p];
    cfg << "\nseed=" << seed << "\ncomparator=" << (boot ? 1 : 0) << "\nboot_panels=" << boot_panels
        << "\nboot_sample=" << boot_sample << "\ngit_commit=" << VCAL_GIT_COMMIT << "\ngit_dirty=" << VCAL_GIT_DIRTY
        << "\nexecutable_sha256=" << exe << "\ndgp_kernel=" << vcal::test::dgp_kernel() << "\n";
    if (!checkpoint_dir.empty() && exe.empty()) {
        std::fprintf(stderr, "--checkpoint: cannot read the executable %s to hash it\n", argv[0]);
        return 2;
    }
    std::unique_ptr<vcal::test::Checkpoint> ck;
    try {
        ck = std::make_unique<vcal::test::Checkpoint>(checkpoint_dir, cfg.str());
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "%s\n", ex.what());
        return 2;
    }

    std::vector<Row> rows(S * R * P);
    std::vector<Seconds> seconds(S * R);  // per job; a resumed job records none (timing samples do not resume)
    for (auto& v : seconds) v.fill(0.0);
    std::vector<std::uint32_t> fit_jobs;  // job = k * R + r, not yet checkpointed
    std::size_t resumed = 0;
    for (std::size_t k = 0; k < S; ++k) {
        for (std::uint32_t r = 0; r < R; ++r) {
            const auto job = static_cast<std::uint32_t>(k * R + r);
            const std::string unit = "s" + std::to_string(ids[k]) + "_r" + std::to_string(r);
            if (ck->has(unit)) {
                try {
                    unpack(ck->load(unit), ids[k], r, &rows[static_cast<std::size_t>(job) * P], P);
                } catch (const std::exception& ex) {
                    std::fprintf(stderr, "%s: %s\n", unit.c_str(), ex.what());
                    return 2;
                }
                ++resumed;
            } else {
                fit_jobs.push_back(job);
            }
        }
    }
    // Comparator slots, per job that owes one: B refits written by their tasks, assembled when the last is in.
    std::vector<std::vector<Refit>> slots(S * R);
    std::vector<std::uint32_t> refits_left(S * R, 0);
    std::mutex m;
    std::deque<std::pair<std::uint32_t, std::uint32_t>> refit_tasks;  // (job, b)
    std::size_t next_fit = 0, jobs_left = fit_jobs.size(), written = 0;
    e::SurfaceRowCache cache;
    const auto t0 = Clock::now();

    const auto finish_job = [&](std::uint32_t job) {  // with m unlocked, once per job
        const std::size_t k = job / R;
        const auto r = static_cast<std::uint32_t>(job % R);
        // Saved under the lock, so --stop-after-units N leaves exactly N units on disk: with the save outside it,
        // two threads could both save before either reached the exit (Copilot's finding on PR 37's S-10 tool).
        const std::lock_guard<std::mutex> lock(m);
        ck->save("s" + std::to_string(ids[k]) + "_r" + std::to_string(r),
                 pack(ids[k], r, &rows[static_cast<std::size_t>(job) * P], P));
        --jobs_left;
        if (stop_after > 0 && ck->enabled() && ++written == stop_after) {
            std::fprintf(stderr, "stopping after %zu replicates (--stop-after-units)\n", written);
            std::fflush(stderr);
            std::_Exit(3);
        }
    };
    // The comparator's intervals and timings, in b order, once every refit of the job is in.
    const auto assemble = [&](std::uint32_t job) {
        Row& row = rows[static_cast<std::size_t>(job) * P + pb];
        std::vector<double> pds, rhos;
        for (const Refit& f : slots[job]) {
            for (int ph = 0; ph < kPhases; ++ph) {
                const double scale = boot_sample > 0 ? static_cast<double>(boot_panels) / boot_sample : 1.0;
                seconds[job][static_cast<std::size_t>(ph)] += f.sec[static_cast<std::size_t>(ph)] *
                                                              (ph == kBootstrap ? scale : 1.0);
            }
            if (f.ok) {
                pds.push_back(f.pd);
                rhos.push_back(f.rho);
            } else {
                ++row.boot_excluded;
            }
        }
        row.boot_used = static_cast<std::uint32_t>(pds.size());
        if (!pds.empty()) {
            std::sort(pds.begin(), pds.end());
            std::sort(rhos.begin(), rhos.end());
            row.boot_lo[0] = resample::quantile_type7(pds, 0.025);
            row.boot_hi[0] = resample::quantile_type7(pds, 0.975);
            row.boot_lo[1] = resample::quantile_type7(rhos, 0.025);
            row.boot_hi[1] = resample::quantile_type7(rhos, 0.975);
        }
        slots[job].clear();
        slots[job].shrink_to_fit();
    };
    const auto worker = [&] {
        for (;;) {
            std::uint32_t job = 0, b = 0;
            bool is_refit = false;
            {
                std::unique_lock<std::mutex> lock(m);
                if (!refit_tasks.empty()) {
                    std::tie(job, b) = refit_tasks.front();
                    refit_tasks.pop_front();
                    is_refit = true;
                } else if (next_fit < fit_jobs.size()) {
                    job = fit_jobs[next_fit++];
                } else if (jobs_left == 0) {
                    return;
                } else {
                    lock.unlock();  // refits may still be released by fits in flight
                    std::this_thread::sleep_for(std::chrono::milliseconds(20));
                    continue;
                }
            }
            const std::size_t k = job / R;
            const auto r = static_cast<std::uint32_t>(job % R);
            const rc::Scenario s = rc::scenario(ids[k]);
            if (!is_refit) {
                Seconds sec{};
                for (std::uint32_t p = 0; p < P; ++p) {
                    std::vector<std::int64_t> n(static_cast<std::size_t>(s.periods), s.obligors), d(n.size());
                    {
                        Timer t{&sec[kSimulate]};
                        const dgp::PanelSpec spec{seed, s.id, r, s.pd, s.rho, n.data(), s.periods};
                        if (dgp::simulate_panel_ar1(spec, phis[p], d.data(), nullptr) != dgp::Status::Ok) std::abort();
                    }
                    rows[static_cast<std::size_t>(job) * P + p] = fit_panel(s, d, cache, sec.data());
                }
                const bool owes = boot && pb < P && in_subset(s.id) && r < kBootReplicates &&
                                  rows[static_cast<std::size_t>(job) * P + pb].fitted != 0u && refits_per_job > 0;
                std::unique_lock<std::mutex> lock(m);
                seconds[job] = sec;
                if (owes) {
                    slots[job].resize(refits_per_job);
                    refits_left[job] = refits_per_job;
                    for (std::uint32_t q = 0; q < refits_per_job; ++q) refit_tasks.emplace_back(job, q);
                } else {
                    lock.unlock();
                    finish_job(job);
                }
            } else {
                const Refit f = refit(s, r, b, rows[static_cast<std::size_t>(job) * P + pb].est, cache);
                std::unique_lock<std::mutex> lock(m);
                slots[job][b] = f;
                if (--refits_left[job] == 0) {
                    assemble(job);
                    lock.unlock();
                    finish_job(job);
                }
            }
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    const double wall = std::chrono::duration<double>(Clock::now() - t0).count();

    char buf[2048];
    std::snprintf(buf, sizeof buf, "# S-35 M6 arm (D-196, D-198): one row per (phi, scenario, replicate); seed 0x%llX; dgp_kernel %s\n",
                  static_cast<unsigned long long>(seed), vcal::test::dgp_kernel().c_str());
    std::string text = buf;
    text += "# profile flags: engine kInterval*; phi's interval on the study subset only; boot_*: the comparator "
            "(subset, replicates 0-99, phi = 0.6, S-35's seed only)\n"
            "phi,scenario,replicate,pd_true,rho_true,periods,obligors,binom_pd,binom_rho,binom_loglik,binom_flags,"
            "fitted,pd,rho,phi_hat,se_pd,se_rho,se_phi,loglik,phi0_pd,phi0_rho,phi0_loglik,doubled_check,"
            "placement_check,scan_local_maxima,scan_unconverged,second_maximum_gap,evaluations,scores,flags,"
            "prof_lo_pd,prof_hi_pd,prof_flags_pd,prof_lo_rho,prof_hi_rho,prof_flags_rho,prof_lo_phi,prof_hi_phi,"
            "prof_flags_phi,prof_residual,boot_lo_pd,boot_hi_pd,boot_lo_rho,boot_hi_rho,boot_used,boot_excluded\n";
    for (std::uint32_t p = 0; p < P; ++p) {
        for (std::size_t k = 0; k < S; ++k) {
            const rc::Scenario s = rc::scenario(ids[k]);
            for (std::uint32_t r = 0; r < R; ++r) {
                const Row& w = rows[(k * R + r) * P + p];
                const auto& x = w.est;
                std::snprintf(buf, sizeof buf,
                              "%.17g,%u,%u,%.17g,%.17g,%lld,%lld,%.17g,%.17g,%.17g,%u,%u,%.17g,%.17g,%.17g,%.17g,%.17g,"
                              "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%.17g,%lld,%lld,%u,%.17g,%.17g,%u,%.17g,"
                              "%.17g,%u,%.17g,%.17g,%u,%.17g,%.17g,%.17g,%.17g,%.17g,%u,%u\n",
                              phis[p], s.id, r, s.pd, s.rho, static_cast<long long>(s.periods),
                              static_cast<long long>(s.obligors), w.binom_pd, w.binom_rho, w.binom_loglik,
                              w.binom_flags, w.fitted, x.value[0], x.value[1], x.value[2], x.se[0], x.se[1],
                              x.se[2], x.loglik, x.phi0_value[0], x.phi0_value[1], x.phi0_loglik, x.doubled_check,
                              x.placement_check, x.scan_local_maxima, x.scan_unconverged, x.second_maximum_gap,
                              static_cast<long long>(x.evaluations), static_cast<long long>(x.scores), x.flags,
                              w.prof.lo[0], w.prof.hi[0], w.prof.flags[0], w.prof.lo[1], w.prof.hi[1], w.prof.flags[1],
                              w.prof.lo[2], w.prof.hi[2], w.prof.flags[2], w.prof.residual_max, w.boot_lo[0],
                              w.boot_hi[0], w.boot_lo[1], w.boot_hi[1], w.boot_used, w.boot_excluded);
                text += buf;
            }
        }
    }
    if (boot_sample > 0) {  // a timing sample: its comparator intervals are wrong as results, so no rows at all
        std::printf("--boot-sample %u: a timing sample; no rows written (%s untouched)\n", boot_sample, out.c_str());
    } else {
        std::ofstream f(out, std::ios::binary);
        f << text;
        f.close();
        if (!f) {
            std::fprintf(stderr, "cannot write %s\n", out.c_str());
            return 1;
        }
    }
    double total[kPhases] = {};
    for (const auto& v : seconds) {
        for (int ph = 0; ph < kPhases; ++ph) total[ph] += v[static_cast<std::size_t>(ph)];
    }
    double all = 0.0;
    for (int ph = 0; ph < kPhases; ++ph) all += total[ph];
    std::printf("%u replicates x %zu scenarios x %u phi in %.0f s (%u threads); CPU s %.0f:", R, S, P, wall, threads, all);
    for (int ph = 0; ph < kPhases; ++ph) std::printf(" %s %.0f (%.0f%%);", kPhaseNames[ph], total[ph], 100.0 * total[ph] / all);
    std::printf("\n");
    if (resumed > 0) std::printf("resumed: %zu of %zu replicates from %s\n", resumed, S * R, checkpoint_dir.c_str());
    if (!timing_out.empty()) {
        // Shares (perf/README.md): every phase runs on every replicate, except the comparator, on replicates 0-99
        // of 1,000 (0.1), its own phase (owner); with --boot-sample its seconds are K refits' scaled to B. The phi
        // interval runs on every replicate of the subset's scenarios only, which the per-cell means carry (other
        // scenarios record 0). The row cache's fills run once.
        std::vector<std::string> phases(kPhaseNames, kPhaseNames + kPhases), shares(kPhases, "1");
        shares[kBootstrap] = "0.1";
        shares[kCacheFill] = "once";
        vcal::test::TimingSample sample("study_lrpd_serial_m6", phases, shares);
        for (std::size_t k = 0; k < S; ++k) {
            const rc::Scenario s = rc::scenario(ids[k]);
            for (std::uint32_t r = 0; r < R; ++r) {
                const Seconds& v = seconds[k * R + r];
                sample.record(s.id, s.periods, s.obligors, r, std::vector<double>(v.begin(), v.end()));
            }
        }
        if (!sample.write(timing_out)) {
            std::fprintf(stderr, "cannot write %s\n", timing_out.c_str());
            return 1;
        }
    }
    return 0;
}
