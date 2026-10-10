// SPDX-License-Identifier: Apache-2.0
//
// A runtime timing sample (perf/README.md, "Projecting a full run", D-187): a study tool run on a few
// replicates of every scenario records, per (scenario, replicate) job, the thread-seconds of each of
// its phases, and perf/project_runtime.py projects the full run from them, weighting each (T, n) cell
// by its count in the matrix. Each phase carries its share of a full run's replicates that execute it
// ("1" for every replicate, "0.05" for replicates 0-49 of 1,000, or "once" for a one-time cost such
// as filling a surface-row cache, counted once rather than per replicate).
//
// One-time costs are recorded apart from per-replicate ones (D-193). A surface-row cache fill is paid once
// per distinct observation for the whole run, not per replicate, so a small sample that charges it to its
// few jobs projects too high (S-35, D-192). A tool measures each job's fills with FillMeter and records
// them in a phase with share "once", subtracted from the phase they happened in.
//
// File format (CSV): '#' comment lines, including one "# share,<phase>,<share>" line per phase and a
// "# tool,<name>" line; then the header scenario,periods,obligors,replicate,<phase>,...; one row per
// job, seconds as %.6f.
#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "engine/surface_cache.hpp"

namespace vcal::test {

// The calling thread's surface-row cache fills since construction: wrap the part of a job that may fill
// rows, then move seconds() from that phase into the job's "cache fill" phase.
class FillMeter {
public:
    FillMeter() : s0_(engine::SurfaceRowCache::thread_fills().seconds) {}
    double seconds() const { return engine::SurfaceRowCache::thread_fills().seconds - s0_; }

private:
    double s0_;
};

class TimingSample {
public:
    // phases and shares in the same order; shares as written to the file ("1", "0.05", "once").
    TimingSample(std::string tool, std::vector<std::string> phases, std::vector<std::string> shares)
        : tool_(std::move(tool)), phases_(std::move(phases)), shares_(std::move(shares)) {}

    std::size_t phases() const { return phases_.size(); }

    // Thread-safe; seconds has phases() values (one per phase, in order), or the call aborts.
    void record(std::uint32_t scenario, std::int64_t periods, std::int64_t obligors, std::uint32_t replicate,
                std::vector<double> seconds) {
        if (seconds.size() != phases_.size()) std::abort();
        Row r{scenario, periods, obligors, replicate, std::move(seconds)};
        const std::lock_guard<std::mutex> lock(mutex_);
        rows_.push_back(std::move(r));
    }

    // Writes the sample, rows in (scenario, replicate) order. False if the file cannot be written.
    bool write(const std::string& path) {
        const std::lock_guard<std::mutex> lock(mutex_);
        std::ofstream f(path, std::ios::binary);
        if (!f) return false;
        f << "# timing sample (D-187): thread-seconds per job and phase\n";
        f << "# tool," << tool_ << '\n';
        for (std::size_t p = 0; p < phases_.size(); ++p) f << "# share," << phases_[p] << ',' << shares_[p] << '\n';
        f << "scenario,periods,obligors,replicate";
        for (const auto& p : phases_) f << ',' << p;
        f << '\n' << std::fixed << std::setprecision(6);
        std::vector<const Row*> order;
        for (const auto& r : rows_) order.push_back(&r);
        std::sort(order.begin(), order.end(), [](const Row* a, const Row* b) { return less(*a, *b); });
        for (const Row* r : order) {
            f << r->scenario << ',' << r->periods << ',' << r->obligors << ',' << r->replicate;
            for (const double s : r->seconds) f << ',' << s;
            f << '\n';
        }
        f.close();
        return static_cast<bool>(f);
    }

private:
    struct Row {
        std::uint32_t scenario;
        std::int64_t periods, obligors;
        std::uint32_t replicate;
        std::vector<double> seconds;
    };
    static bool less(const Row& a, const Row& b) {
        return a.scenario != b.scenario ? a.scenario < b.scenario : a.replicate < b.replicate;
    }

    std::string tool_;
    std::vector<std::string> phases_, shares_;
    std::vector<Row> rows_;
    std::mutex mutex_;
};

}  // namespace vcal::test
