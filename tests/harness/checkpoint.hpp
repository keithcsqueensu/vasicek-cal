// SPDX-License-Identifier: Apache-2.0
//
// Checkpoint and resume for long study runs (P-11, D-179). A run writes each completed unit of work
// (a scenario) to its checkpoint directory as it finishes; a restarted run skips the units already
// written. What makes that safe:
//   - the directory records the run's configuration: the tool's name and arguments that affect the
//     results, the source commit and whether the tree was dirty, and the SHA-256 of the tool's own
//     executable. A run refuses to resume from a directory written under any other configuration, so
//     a checkpoint is only ever combined with output of the same code and settings (a rebuild of the
//     tool, even of the same source, is a different configuration: conservative by design);
//   - each unit is written to a temporary file and renamed into place, so a unit is either complete
//     on disk or absent, never partial, whenever the run is killed;
//   - a unit's bytes are exactly what the run would have kept in memory, and the final output is
//     built from the units after they are all present, so a resumed run's output is identical to an
//     uninterrupted one (tested by killing a run midway and resuming: study_parametric_bootstrap_resume).
#pragma once

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "tests/harness/sha256.hpp"

namespace vcal::test {

class Checkpoint {
public:
    // dir empty: checkpointing off. Otherwise the directory is created if needed; if it already holds
    // a configuration, it must equal `config` (else std::runtime_error).
    Checkpoint(const std::string& dir, const std::string& config) : dir_(dir) {
        if (dir_.empty()) return;
        namespace fs = std::filesystem;
        fs::create_directories(dir_);
        const fs::path cfg = fs::path(dir_) / "CONFIG";
        if (fs::exists(cfg)) {
            std::ifstream in(cfg, std::ios::binary);
            const std::string old((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            if (old != config) {
                throw std::runtime_error("checkpoint directory " + dir_ +
                                         " was written under a different configuration; refusing to resume\n--- "
                                         "recorded ---\n" + old + "--- this run ---\n" + config);
            }
            resumed_ = true;
        } else {
            write_atomic(cfg, std::vector<std::uint8_t>(config.begin(), config.end()));
        }
    }

    bool enabled() const { return !dir_.empty(); }
    bool resumed() const { return resumed_; }

    bool has(const std::string& unit) const {
        return enabled() && std::filesystem::exists(std::filesystem::path(dir_) / (unit + ".bin"));
    }

    std::vector<std::uint8_t> load(const std::string& unit) const {
        std::ifstream in(std::filesystem::path(dir_) / (unit + ".bin"), std::ios::binary);
        return std::vector<std::uint8_t>((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    }

    void save(const std::string& unit, const std::vector<std::uint8_t>& bytes) const {
        if (enabled()) write_atomic(std::filesystem::path(dir_) / (unit + ".bin"), bytes);
    }

private:
    static void write_atomic(const std::filesystem::path& path, const std::vector<std::uint8_t>& bytes) {
        const std::filesystem::path tmp = path.string() + ".tmp";
        {
            std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
            out.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
            out.flush();
            if (!out) throw std::runtime_error("cannot write checkpoint " + tmp.string());
        }
        std::filesystem::rename(tmp, path);  // replaces atomically on POSIX and on Windows (MoveFileEx)
    }

    std::string dir_;
    bool resumed_ = false;
};

// The SHA-256 of a file (the tool's executable, for the configuration), or "" if it cannot be read.
inline std::string file_sha256(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return "";
    const std::vector<std::uint8_t> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    return sha256_hex(bytes);
}

}  // namespace vcal::test
