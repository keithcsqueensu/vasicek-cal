// SPDX-License-Identifier: Apache-2.0
//
// The C ABI (include/vcal/vcal.h; M2c, D-138..D-144).
//
// Every entry point runs its body inside guard(), which clears the calling thread's error message,
// turns every C++ exception into a status and a message, and never lets one escape. Inputs are
// validated here, with messages naming the offending field, before anything reaches the engine;
// the engine's own checks remain as a second line, and a disagreement between the two is reported
// as VCAL_E_INTERNAL. Outputs are computed in full before any output struct is written.
//
// v0 serves one tuple: BinomialMixture<PrecisionF64> with the parity rule (D-118) and its doubled
// check, the ArgMax reducer and the CPU backend. The dispatch registry of ARCHITECTURE.md §3
// arrives with the second tuple.
#include "vcal/vcal.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "abi/build_info.hpp"
#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "core/special/lbinom.hpp"
#include "dgp/dgp.hpp"
#include "engine/calibrate.hpp"
#include "engine/profile.hpp"
#include "engine/refine.hpp"
#include "engine/surface.hpp"
#include "resample/bootstrap.hpp"
#include "resample/weights.hpp"

struct vcal_context {
    std::int32_t profile;
    std::int32_t n_threads;
};

namespace {

namespace e = vcal::engine;
namespace rs = vcal::resample;
using Objective = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
using Obs = Objective::Obs;

// ---- layout: every struct is its fields with no implicit padding (tests/abi checks offsets) ----

static_assert(sizeof(void*) == 8, "the vcal ABI is defined for 64-bit platforms");
static_assert(sizeof(vcal_context_options) == 16);
static_assert(sizeof(vcal_panel) == 32);
static_assert(sizeof(vcal_grid) == 56);
static_assert(sizeof(vcal_estimate) == 96);
static_assert(sizeof(vcal_profile_intervals) == 88);
static_assert(sizeof(vcal_resample_spec) == 88);
static_assert(sizeof(vcal_replicates) == 56);
static_assert(sizeof(vcal_percentile_intervals) == 80);
static_assert(sizeof(vcal_dgp_spec) == 56);

// ---- the ABI's constants are the engine's ----

static_assert(std::uint32_t{VCAL_FLAG_GRID_EDGE} == std::uint32_t{e::kFlagGridEdge});
static_assert(std::uint32_t{VCAL_FLAG_FLAT_SURFACE} == std::uint32_t{e::kFlagFlatSurface});
static_assert(std::uint32_t{VCAL_FLAG_QUADRATURE_UNCONVERGED} == std::uint32_t{e::kFlagQuadratureUnconverged});
static_assert(std::uint32_t{VCAL_FLAG_REFINEMENT_REJECTED} == std::uint32_t{e::kFlagRefinementRejected});
static_assert(std::uint32_t{VCAL_FLAG_NUMERIC} == std::uint32_t{e::kFlagNumeric});
static_assert(std::uint32_t{VCAL_FLAG_NEAR_BOUND} == std::uint32_t{e::kFlagNearBound});
static_assert(std::uint32_t{VCAL_INTERVAL_LOWER_TRUNCATED} == std::uint32_t{e::kIntervalLowerTruncated});
static_assert(std::uint32_t{VCAL_INTERVAL_UPPER_TRUNCATED} == std::uint32_t{e::kIntervalUpperTruncated});
static_assert(std::uint32_t{VCAL_INTERVAL_NOT_COMPUTED} == std::uint32_t{e::kIntervalNotComputed});
static_assert(VCAL_SCALE_LINEAR == static_cast<int>(vcal::AxisScale::Linear));
static_assert(VCAL_SCALE_LOG == static_cast<int>(vcal::AxisScale::Log));
static_assert(VCAL_SCALE_LOGIT == static_cast<int>(vcal::AxisScale::Logit));
static_assert(VCAL_SCALE_PROBIT == static_cast<int>(vcal::AxisScale::Probit));
static_assert(VCAL_RESAMPLE_IID_BOOTSTRAP == static_cast<int>(rs::Scheme::IidBootstrap));
static_assert(VCAL_RESAMPLE_BLOCK_BOOTSTRAP == static_cast<int>(rs::Scheme::BlockBootstrap));

// ---- errors ----

// The calling thread's message. A fixed buffer, so recording a failure cannot itself fail.
constexpr std::size_t kErrorCapacity = 1024;
thread_local char t_error[kErrorCapacity];

class AbiError : public std::runtime_error {
public:
    AbiError(vcal_status s, const std::string& message) : std::runtime_error(message), status(s) {}
    vcal_status status;
};

[[noreturn]] void invalid(const std::string& message) { throw AbiError(VCAL_E_INVALID_ARGUMENT, message); }

std::string num(double x) {
    char b[32];
    std::snprintf(b, sizeof b, "%.17g", x);
    return b;
}
std::string num(std::int64_t x) { return std::to_string(x); }
std::string num(std::uint64_t x) { return std::to_string(x); }
std::string num(std::int32_t x) { return std::to_string(x); }
std::string num(std::uint32_t x) { return std::to_string(x); }

vcal_status record(vcal_status status, const char* function, const char* message) noexcept {
    std::snprintf(t_error, kErrorCapacity, "%s: %s", function, message);
    return status;
}

template <class F>
vcal_status guard(const char* function, F&& body) noexcept {
    t_error[0] = '\0';
    try {
        body();
        return VCAL_OK;
    } catch (const AbiError& x) {
        return record(x.status, function, x.what());
    } catch (const std::bad_alloc&) {
        return record(VCAL_E_OUT_OF_MEMORY, function, "out of memory");
    } catch (const std::length_error&) {
        return record(VCAL_E_OUT_OF_MEMORY, function, "requested allocation is too large");
    } catch (const std::invalid_argument& x) {  // from resample/weights.hpp's own checks
        return record(VCAL_E_INVALID_ARGUMENT, function, x.what());
    } catch (const std::exception& x) {
        return record(VCAL_E_INTERNAL, function, x.what());
    } catch (...) {
        return record(VCAL_E_INTERNAL, function, "unknown exception");
    }
}

// ---- struct versioning (vcal.h, "Structs") ----

// The size of each struct in ABI 0.1, the smallest a caller may pass. Later minor versions append
// fields, so sizeof grows while these stay.
template <class T>
constexpr std::uint32_t kMinSize = static_cast<std::uint32_t>(sizeof(T));

template <class T>
void check_size(const T* p, const char* name) {
    if (p == nullptr) invalid(std::string(name) + " is NULL");
    if (p->struct_size < kMinSize<T>) {
        invalid(std::string(name) + ".struct_size = " + num(p->struct_size) + " is below " + num(kMinSize<T>) +
                ", the size in ABI 0.1: set it to sizeof the struct");
    }
}

// A copy of an input struct. Bytes a newer caller appended must be zero: a nonzero one is a
// setting this library does not know, and ignoring it silently would change the result.
template <class T>
T read_in(const T* p, const char* name) {
    check_size(p, name);
    if (p->struct_size > sizeof(T)) {
        const auto* bytes = reinterpret_cast<const unsigned char*>(p);
        for (std::size_t i = sizeof(T); i < p->struct_size; ++i) {
            if (bytes[i] != 0) {
                throw AbiError(VCAL_E_UNSUPPORTED, std::string(name) + " has a nonzero byte at offset " +
                                                       num(static_cast<std::uint64_t>(i)) +
                                                       ", beyond the fields this library knows (ABI 0.1)");
            }
        }
    }
    T v{};
    std::memcpy(&v, p, std::min<std::size_t>(p->struct_size, sizeof(T)));
    return v;
}

void check_reserved(std::uint32_t reserved, const char* name) {
    if (reserved != 0) {
        throw AbiError(VCAL_E_UNSUPPORTED, std::string(name) + ".reserved = " + num(reserved) + ": must be 0");
    }
}

// Writes the fields this library knows; the caller's struct_size is kept.
template <class T>
void write_out(T* p, T v) {
    v.struct_size = p->struct_size;
    std::memcpy(p, &v, std::min<std::size_t>(p->struct_size, sizeof(T)));
}

// Variable-size outputs (vcal.h, "Memory"). Returns false for a size query.
bool check_buffer(const void* buffer, std::int64_t capacity, std::int64_t needed, std::int64_t* required,
                  const char* name) {
    if (required != nullptr) *required = needed;
    if (buffer == nullptr) {
        if (required == nullptr) invalid(std::string(name) + " and required are both NULL");
        return false;
    }
    if (capacity < needed) {
        throw AbiError(VCAL_E_BUFFER_TOO_SMALL,
                       std::string(name) + ": capacity " + num(capacity) + " is below the " + num(needed) + " needed");
    }
    return true;
}

void copy_string(const char* s, char* buffer, std::int64_t capacity, std::int64_t* required) {
    const auto needed = static_cast<std::int64_t>(std::strlen(s)) + 1;
    if (check_buffer(buffer, capacity, needed, required, "buffer")) std::memcpy(buffer, s, static_cast<std::size_t>(needed));
}

// ---- inputs ----

const vcal_context& context_of(const vcal_context* c) {
    if (c == nullptr) invalid("context is NULL");
    return *c;
}

vcal::backends::CpuBackend backend(const vcal_context& c) { return vcal::backends::CpuBackend{c.n_threads}; }

std::vector<Obs> read_panel(const vcal_panel* p) {
    const vcal_panel v = read_in(p, "panel");
    check_reserved(v.reserved, "panel");
    if (v.n_periods < 1) invalid("panel.n_periods = " + num(v.n_periods) + ": need at least 1");
    if (v.n_obligors == nullptr) invalid("panel.n_obligors is NULL");
    if (v.n_defaults == nullptr) invalid("panel.n_defaults is NULL");
    std::vector<Obs> obs(static_cast<std::size_t>(v.n_periods));
    for (std::int64_t t = 0; t < v.n_periods; ++t) {
        const std::int64_t n = v.n_obligors[t];
        const std::int64_t d = v.n_defaults[t];
        if (n < 1 || n > vcal::special::kLbinomMaxN) {
            invalid("panel.n_obligors[" + num(t) + "] = " + num(n) + ": need 1 <= n <= 2^53");
        }
        if (d < 0 || d > n) {
            invalid("panel.n_defaults[" + num(t) + "] = " + num(d) + ": need 0 <= d <= n_obligors[" + num(t) +
                    "] = " + num(n));
        }
        obs[static_cast<std::size_t>(t)] = {n, d};
    }
    return obs;
}

vcal::Axis read_axis(double lo, double hi, std::int32_t points, std::int32_t scale, const char* name) {
    if (scale < VCAL_SCALE_LINEAR || scale > VCAL_SCALE_PROBIT) {
        invalid(std::string("grid.") + name + "_scale = " + num(scale) + " is not a VCAL_SCALE_* value");
    }
    // Both parameters are probabilities, whatever the axis scale.
    if (!(lo > 0.0 && lo < hi && hi < 1.0)) {
        invalid(std::string("grid.") + name + "_lo = " + num(lo) + ", " + name + "_hi = " + num(hi) +
                ": need 0 < lo < hi < 1");
    }
    const vcal::Axis a{lo, hi, points, static_cast<vcal::AxisScale>(scale)};
    if (const char* msg = vcal::axis_error(a)) invalid(std::string("grid, ") + name + " axis: " + msg);
    return a;
}

vcal::Grid<2> read_grid(const vcal_grid* g) {
    const vcal_grid v = read_in(g, "grid");
    check_reserved(v.reserved, "grid");
    return {{read_axis(v.pd_lo, v.pd_hi, v.pd_points, v.pd_scale, "pd"),
             read_axis(v.rho_lo, v.rho_hi, v.rho_points, v.rho_scale, "rho")}};
}

// T x K doubles must be addressable.
std::int64_t surface_cells(std::int64_t periods, const vcal::Grid<2>& g) {
    const std::int64_t K = g.size();
    if (K > std::numeric_limits<std::int64_t>::max() / 8 / periods) {
        invalid("the surface would have " + num(periods) + " x " + num(K) + " cells, too many to address");
    }
    return periods * K;
}

const vcal::quadrature::SplitRule<vcal::PrecisionF64>& primary_rule() {
    static const auto r = vcal::quadrature::parity_rule();
    return r;
}
const vcal::quadrature::SplitRule<vcal::PrecisionF64>& check_rule() {
    static const auto r = vcal::quadrature::parity_rule(true);
    return r;
}

std::vector<double> surface(const vcal_context& c, const std::vector<Obs>& obs, const vcal::Grid<2>& g) {
    const auto T = static_cast<std::int64_t>(obs.size());
    std::vector<double> L(static_cast<std::size_t>(surface_cells(T, g)));
    e::evaluate_surface(backend(c), Objective{}, primary_rule(), obs.data(), T, g, L.data());
    return L;
}

// ---- resampling specs ----

struct Resample {
    vcal_resample_spec spec;
    std::int64_t periods;
    std::int64_t replicates;
    std::int64_t block_length;  // resolved default
    double level;               // resolved default
};

const char* scheme_name(std::int32_t s) {
    switch (s) {
        case VCAL_RESAMPLE_IID_BOOTSTRAP: return "VCAL_RESAMPLE_IID_BOOTSTRAP";
        case VCAL_RESAMPLE_BLOCK_BOOTSTRAP: return "VCAL_RESAMPLE_BLOCK_BOOTSTRAP";
        case VCAL_RESAMPLE_JACKKNIFE: return "VCAL_RESAMPLE_JACKKNIFE";
        case VCAL_RESAMPLE_WALK_FORWARD: return "VCAL_RESAMPLE_WALK_FORWARD";
        case VCAL_RESAMPLE_INDICES: return "VCAL_RESAMPLE_INDICES";
        case VCAL_RESAMPLE_WEIGHTS: return "VCAL_RESAMPLE_WEIGHTS";
        default: return nullptr;
    }
}

// The smallest l >= 1 with l^3 >= T: ceil(T^(1/3)) in integers (D-133), so no cube root can round.
std::int64_t default_block_length(std::int64_t periods) {
    std::int64_t l = 1;
    while (l * l * l < periods) ++l;
    return l;
}

Resample read_resample(const vcal_resample_spec* p, std::int64_t periods) {
    Resample r{read_in(p, "spec"), periods, 0, 0, 0.95};
    const vcal_resample_spec& s = r.spec;
    const char* scheme = scheme_name(s.scheme);
    if (scheme == nullptr) invalid("spec.scheme = " + num(s.scheme) + " is not a VCAL_RESAMPLE_* value");

    // A field the scheme does not use must be zero: a stray value is a caller's mistake.
    const auto unused = [&](bool is_zero, const char* field) {
        if (!is_zero) invalid(std::string("spec.") + field + " is not used by " + scheme + " and must be 0");
    };
    const bool bootstrap = s.scheme == VCAL_RESAMPLE_IID_BOOTSTRAP || s.scheme == VCAL_RESAMPLE_BLOCK_BOOTSTRAP;
    const bool walk = s.scheme == VCAL_RESAMPLE_WALK_FORWARD;
    if (!bootstrap) unused(s.seed == 0, "seed");
    if (s.scheme != VCAL_RESAMPLE_BLOCK_BOOTSTRAP) unused(s.block_length == 0, "block_length");
    if (!walk) {
        unused(s.window == 0, "window");
        unused(s.first_end == 0, "first_end");
        unused(s.step == 0, "step");
    }
    if (s.scheme != VCAL_RESAMPLE_INDICES) {
        unused(s.draws == 0, "draws");
        unused(s.indices == nullptr, "indices");
    }
    if (s.scheme != VCAL_RESAMPLE_WEIGHTS) unused(s.weights == nullptr, "weights");
    if (s.scheme == VCAL_RESAMPLE_JACKKNIFE || walk) unused(s.replicates == 0, "replicates");

    if (s.level != 0.0) {
        if (!(s.level > 0.0 && s.level < 1.0)) invalid("spec.level = " + num(s.level) + ": need 0 < level < 1, or 0 for 0.95");
        r.level = s.level;
    }

    const std::int64_t max_cells = std::numeric_limits<std::int64_t>::max() / 8;
    switch (s.scheme) {
        case VCAL_RESAMPLE_IID_BOOTSTRAP:
        case VCAL_RESAMPLE_BLOCK_BOOTSTRAP:
            if (s.replicates < 1 || s.replicates > std::int64_t{0xFFFFFFFF}) {
                invalid("spec.replicates = " + num(s.replicates) + ": need 1 <= B < 2^32");
            }
            if (periods > (std::int64_t{1} << 31)) invalid("the bootstrap needs at most 2^31 periods");
            r.replicates = s.replicates;
            if (s.scheme == VCAL_RESAMPLE_BLOCK_BOOTSTRAP) {
                r.block_length = s.block_length == 0 ? default_block_length(periods) : s.block_length;
                if (r.block_length < 1 || r.block_length > periods) {
                    invalid("spec.block_length = " + num(s.block_length) + ": need 1 <= l <= n_periods = " +
                            num(periods) + ", or 0 for ceil(T^(1/3))");
                }
            }
            break;
        case VCAL_RESAMPLE_JACKKNIFE:
            if (periods < 2) invalid("the jackknife needs at least 2 periods");
            r.replicates = periods;
            break;
        case VCAL_RESAMPLE_WALK_FORWARD:
            if (s.step < 1) invalid("spec.step = " + num(s.step) + ": need step >= 1");
            if (s.first_end < 1 || s.first_end > periods) {
                invalid("spec.first_end = " + num(s.first_end) + ": need 1 <= first_end <= n_periods = " + num(periods));
            }
            if (s.window < 0 || s.window > s.first_end) {
                invalid("spec.window = " + num(s.window) + ": need 0 <= window <= first_end = " + num(s.first_end));
            }
            r.replicates = (periods - s.first_end) / s.step + 1;
            break;
        case VCAL_RESAMPLE_INDICES:
            if (s.replicates < 1) invalid("spec.replicates = " + num(s.replicates) + ": need at least 1");
            if (s.draws < 1) invalid("spec.draws = " + num(s.draws) + ": need at least 1");
            if (s.indices == nullptr) invalid("spec.indices is NULL");
            if (s.draws > max_cells / s.replicates) invalid("spec.replicates x spec.draws is too large to address");
            r.replicates = s.replicates;
            break;
        case VCAL_RESAMPLE_WEIGHTS:
            if (s.replicates < 1) invalid("spec.replicates = " + num(s.replicates) + ": need at least 1");
            if (s.weights == nullptr) invalid("spec.weights is NULL");
            r.replicates = s.replicates;
            break;
        default: break;
    }
    if (periods > max_cells / r.replicates) invalid("the weight matrix, replicates x n_periods, is too large to address");
    if (s.scheme == VCAL_RESAMPLE_WEIGHTS) {
        for (std::int64_t b = 0; b < r.replicates; ++b) {
            bool positive = false;
            for (std::int64_t t = 0; t < periods; ++t) {
                const double w = s.weights[b * periods + t];
                if (!(std::isfinite(w) && w >= 0.0)) {
                    invalid("spec.weights[" + num(b) + " * n_periods + " + num(t) + "] = " + num(w) +
                            ": weights must be finite and >= 0");
                }
                positive = positive || w > 0.0;
            }
            if (!positive) invalid("spec.weights row " + num(b) + " has no positive weight");
        }
    }
    return r;
}

// W, B x T row-major. For VCAL_RESAMPLE_WEIGHTS the caller's matrix is used in place (w_owned
// stays empty).
const double* weights_of(const Resample& r, std::vector<double>& w_owned) {
    const vcal_resample_spec& s = r.spec;
    const std::int64_t T = r.periods;
    switch (s.scheme) {
        case VCAL_RESAMPLE_IID_BOOTSTRAP:
        case VCAL_RESAMPLE_BLOCK_BOOTSTRAP: {
            const auto scheme = static_cast<rs::Scheme>(s.scheme);
            const auto idx = rs::bootstrap_indices(s.seed, scheme, static_cast<std::uint32_t>(r.replicates), T,
                                                   s.scheme == VCAL_RESAMPLE_BLOCK_BOOTSTRAP ? r.block_length : 1);
            w_owned = rs::weights_from_indices(idx.data(), r.replicates, T, T);
            break;
        }
        case VCAL_RESAMPLE_JACKKNIFE: w_owned = rs::jackknife_weights(T); break;
        case VCAL_RESAMPLE_WALK_FORWARD: w_owned = rs::walk_forward_weights(T, s.window, s.first_end, s.step); break;
        case VCAL_RESAMPLE_INDICES: w_owned = rs::weights_from_indices(s.indices, r.replicates, s.draws, T); break;
        case VCAL_RESAMPLE_WEIGHTS: return s.weights;
        default: throw std::logic_error("unreachable resampling scheme");
    }
    if (static_cast<std::int64_t>(w_owned.size()) != r.replicates * T) {
        throw std::logic_error("weight matrix has " + num(static_cast<std::int64_t>(w_owned.size())) +
                               " entries, expected " + num(r.replicates * T));
    }
    return w_owned.data();
}

}  // namespace

// =============================================================================================
extern "C" {

VCAL_API uint32_t VCAL_CALL vcal_abi_version(void) { return VCAL_ABI_VERSION; }

VCAL_API const char* VCAL_CALL vcal_status_string(vcal_status status) {
    switch (status) {
        case VCAL_OK: return "VCAL_OK";
        case VCAL_E_INVALID_ARGUMENT: return "VCAL_E_INVALID_ARGUMENT";
        case VCAL_E_BUFFER_TOO_SMALL: return "VCAL_E_BUFFER_TOO_SMALL";
        case VCAL_E_UNSUPPORTED: return "VCAL_E_UNSUPPORTED";
        case VCAL_E_NUMERIC: return "VCAL_E_NUMERIC";
        case VCAL_E_OUT_OF_MEMORY: return "VCAL_E_OUT_OF_MEMORY";
        case VCAL_E_FP_ENVIRONMENT: return "VCAL_E_FP_ENVIRONMENT";
        case VCAL_E_INTERNAL: return "VCAL_E_INTERNAL";
        default: return "unknown vcal_status";
    }
}

VCAL_API vcal_status VCAL_CALL vcal_last_error(char* buffer, int64_t capacity, int64_t* required) {
    // Not guarded: guard() would clear the message this call exists to read.
    const auto needed = static_cast<std::int64_t>(std::strlen(t_error)) + 1;
    if (required != nullptr) *required = needed;
    if (buffer == nullptr) return required != nullptr ? VCAL_OK : VCAL_E_INVALID_ARGUMENT;
    if (capacity < needed) return VCAL_E_BUFFER_TOO_SMALL;
    std::memcpy(buffer, t_error, static_cast<std::size_t>(needed));
    return VCAL_OK;
}

VCAL_API vcal_status VCAL_CALL vcal_build_info(char* buffer, int64_t capacity, int64_t* required) {
    return guard("vcal_build_info", [&] { copy_string(vcal::abi::build_info(), buffer, capacity, required); });
}

VCAL_API vcal_status VCAL_CALL vcal_context_create(const vcal_context_options* options, vcal_context** context) {
    return guard("vcal_context_create", [&] {
        if (context == nullptr) invalid("context is NULL");
        *context = nullptr;
        vcal_context_options o{};
        if (options != nullptr) {
            o = read_in(options, "options");
            check_reserved(o.reserved, "options");
        }
        if (o.profile != VCAL_PROFILE_PARITY) invalid("options.profile = " + num(o.profile) + " is not a VCAL_PROFILE_* value");
        if (o.n_threads < 0) invalid("options.n_threads = " + num(o.n_threads) + ": need >= 0 (0 = the OpenMP default)");
        *context = new vcal_context{o.profile, o.n_threads};
    });
}

VCAL_API vcal_status VCAL_CALL vcal_context_destroy(vcal_context* context) {
    return guard("vcal_context_destroy", [&] { delete context; });
}

VCAL_API vcal_status VCAL_CALL vcal_grid_default(vcal_grid* grid) {
    return guard("vcal_grid_default", [&] {
        check_size(grid, "grid");
        vcal_grid g{};
        g.pd_lo = 1e-4;
        g.pd_hi = 0.2;
        g.pd_points = 61;
        g.pd_scale = VCAL_SCALE_LOGIT;
        g.rho_lo = 1e-3;
        g.rho_hi = vcal::kDefaultRhoUpper;
        g.rho_points = 41;
        g.rho_scale = VCAL_SCALE_LOGIT;
        write_out(grid, g);
    });
}

VCAL_API vcal_status VCAL_CALL vcal_grid_values(const vcal_grid* grid, int32_t axis, double* values, int64_t capacity,
                                                int64_t* required) {
    return guard("vcal_grid_values", [&] {
        const vcal::Grid<2> g = read_grid(grid);
        if (axis != VCAL_AXIS_PD && axis != VCAL_AXIS_RHO) invalid("axis = " + num(axis) + " is not a VCAL_AXIS_* value");
        const vcal::Axis& a = g.axis[axis];
        if (!check_buffer(values, capacity, a.n, required, "values")) return;
        for (std::int32_t i = 0; i < a.n; ++i) values[i] = a.value_at(i);
    });
}

VCAL_API vcal_status VCAL_CALL vcal_calibrate(vcal_context* context, const vcal_panel* panel, const vcal_grid* grid,
                                              vcal_estimate* estimate, vcal_profile_intervals* profile) {
    return guard("vcal_calibrate", [&] {
        const vcal_context& c = context_of(context);
        const std::vector<Obs> obs = read_panel(panel);
        const vcal::Grid<2> g = read_grid(grid);
        check_size(estimate, "estimate");
        if (profile != nullptr) check_size(profile, "profile");
        const auto T = static_cast<std::int64_t>(obs.size());
        surface_cells(T, g);

        std::vector<double> L;
        e::Estimate2 est{};
        switch (e::calibrate(backend(c), Objective{}, primary_rule(), check_rule(), obs.data(), T, g, L, est)) {
            case e::Status::Ok: break;
            case e::Status::SurfaceUndefined:
                throw AbiError(VCAL_E_NUMERIC, "the log-likelihood is not finite at any grid point");
            default: throw std::logic_error("the engine rejected a grid or panel that the ABI accepted");
        }
        vcal_estimate out{};
        out.flags = est.flags;
        out.pd = est.value[0];
        out.rho = est.value[1];
        out.se_pd = est.se[0];
        out.se_rho = est.se[1];
        out.corr_pd_rho = est.corr;
        out.loglik = est.loglik;
        out.quad_check_max = est.quad_check_max;
        out.quad_check_total = est.quad_check_total;
        out.quad_check_flagged = est.quad_check_flagged;
        out.grid_index = est.grid_index;
        out.nan_count = est.nan_count;

        vcal_profile_intervals prof{};
        if (profile != nullptr) {
            const e::ProfileIntervals2 p = e::profile_intervals(Objective{}, primary_rule(), obs.data(), T, g, L, est);
            prof.pd_flags = p.flags[0];
            prof.rho_flags = p.flags[1];
            prof.pd_lo = p.lo[0];
            prof.pd_hi = p.hi[0];
            prof.rho_lo = p.lo[1];
            prof.rho_hi = p.hi[1];
            prof.loglik_max = p.loglik_max;
            prof.pd_at_max = p.max_at[0];
            prof.rho_at_max = p.max_at[1];
            prof.residual_max = p.residual_max;
            prof.evaluations = p.evaluations;
        }
        write_out(estimate, out);
        if (profile != nullptr) write_out(profile, prof);
    });
}

VCAL_API vcal_status VCAL_CALL vcal_surface(vcal_context* context, const vcal_panel* panel, const vcal_grid* grid,
                                            double* surface_out, int64_t capacity, int64_t* required) {
    return guard("vcal_surface", [&] {
        const vcal_context& c = context_of(context);
        const std::vector<Obs> obs = read_panel(panel);
        const vcal::Grid<2> g = read_grid(grid);
        const auto T = static_cast<std::int64_t>(obs.size());
        if (!check_buffer(surface_out, capacity, surface_cells(T, g), required, "surface")) return;
        e::evaluate_surface(backend(c), Objective{}, primary_rule(), obs.data(), T, g, surface_out);
    });
}

VCAL_API vcal_status VCAL_CALL vcal_resample(vcal_context* context, const vcal_panel* panel, const vcal_grid* grid,
                                             const vcal_resample_spec* spec, vcal_replicates* replicates,
                                             vcal_percentile_intervals* intervals, int64_t* required) {
    return guard("vcal_resample", [&] {
        const vcal_context& c = context_of(context);
        const std::vector<Obs> obs = read_panel(panel);
        const vcal::Grid<2> g = read_grid(grid);
        const auto T = static_cast<std::int64_t>(obs.size());
        const Resample r = read_resample(spec, T);
        const std::int64_t B = r.replicates;
        vcal_replicates arrays{};
        if (replicates != nullptr) {
            arrays = read_in(replicates, "replicates");
            check_reserved(arrays.reserved, "replicates");
        }
        if (intervals != nullptr) check_size(intervals, "intervals");
        if (required != nullptr) *required = B;
        if (replicates == nullptr && intervals == nullptr) {
            if (required == nullptr) invalid("replicates, intervals and required are all NULL");
            return;
        }
        if (replicates != nullptr && arrays.capacity < B) {
            throw AbiError(VCAL_E_BUFFER_TOO_SMALL, "replicates.capacity = " + num(arrays.capacity) + " is below the " +
                                                        num(B) + " replicates");
        }

        std::vector<double> w_owned;
        const double* W = weights_of(r, w_owned);
        const std::vector<double> L = surface(c, obs, g);
        std::vector<rs::Replicate2> reps(static_cast<std::size_t>(B));
        rs::replicate_estimates(backend(c), g, L.data(), T, W, B, reps.data());

        vcal_percentile_intervals iv{};
        if (intervals != nullptr) {
            const rs::PercentileInterval p0 = rs::percentile_interval(reps.data(), B, 0, r.level);
            const rs::PercentileInterval p1 = rs::percentile_interval(reps.data(), B, 1, r.level);
            iv.level = r.level;
            iv.pd_lo = p0.lo;
            iv.pd_hi = p0.hi;
            iv.rho_lo = p1.lo;
            iv.rho_hi = p1.hi;
            iv.replicates = B;
            iv.pd_excluded = p0.excluded;
            iv.rho_excluded = p1.excluded;
            for (const auto& x : reps) iv.grid_edge += (x.flags & e::kFlagGridEdge) ? 1 : 0;
        }
        if (replicates != nullptr) {
            for (std::int64_t b = 0; b < B; ++b) {
                const rs::Replicate2& x = reps[static_cast<std::size_t>(b)];
                if (arrays.pd != nullptr) arrays.pd[b] = x.value[0];
                if (arrays.rho != nullptr) arrays.rho[b] = x.value[1];
                if (arrays.grid_loglik != nullptr) arrays.grid_loglik[b] = x.surface_max;
                if (arrays.grid_index != nullptr) arrays.grid_index[b] = x.grid_index;
                if (arrays.flags != nullptr) arrays.flags[b] = x.flags;
            }
        }
        if (intervals != nullptr) write_out(intervals, iv);
    });
}

VCAL_API vcal_status VCAL_CALL vcal_resample_weights(vcal_context* context, const vcal_resample_spec* spec,
                                                     int64_t n_periods, double* weights, int64_t capacity,
                                                     int64_t* required) {
    return guard("vcal_resample_weights", [&] {
        context_of(context);
        if (n_periods < 1) invalid("n_periods = " + num(n_periods) + ": need at least 1");
        const Resample r = read_resample(spec, n_periods);
        if (!check_buffer(weights, capacity, r.replicates * n_periods, required, "weights")) return;
        std::vector<double> w_owned;
        const double* W = weights_of(r, w_owned);
        std::copy(W, W + r.replicates * n_periods, weights);
    });
}

VCAL_API vcal_status VCAL_CALL vcal_dgp_simulate(vcal_context* context, const vcal_dgp_spec* spec, int64_t* n_defaults,
                                                 double* factors, int64_t capacity, int64_t* required) {
    return guard("vcal_dgp_simulate", [&] {
        context_of(context);
        const vcal_dgp_spec s = read_in(spec, "spec");
        check_reserved(s.reserved, "spec");
        if (!(s.pd > 0.0 && s.pd < 1.0)) invalid("spec.pd = " + num(s.pd) + ": need 0 < pd < 1");
        if (!(s.rho > 0.0 && s.rho < 1.0)) invalid("spec.rho = " + num(s.rho) + ": need 0 < rho < 1");
        if (s.n_periods < 1 || s.n_periods > (std::int64_t{1} << 32)) {
            invalid("spec.n_periods = " + num(s.n_periods) + ": need 1 <= n_periods <= 2^32");
        }
        if (s.n_obligors == nullptr) invalid("spec.n_obligors is NULL");
        for (std::int64_t t = 0; t < s.n_periods; ++t) {
            if (s.n_obligors[t] < 0 || s.n_obligors[t] > vcal::dgp::kMaxObligors) {
                invalid("spec.n_obligors[" + num(t) + "] = " + num(s.n_obligors[t]) + ": need 0 <= n <= 2^33 - 2^17");
            }
        }
        if (n_defaults == nullptr && factors != nullptr) invalid("n_defaults is NULL but factors is not");
        if (!check_buffer(n_defaults, capacity, s.n_periods, required, "n_defaults")) return;
        const vcal::dgp::PanelSpec ps{s.seed, s.scenario, s.replicate, s.pd, s.rho, s.n_obligors, s.n_periods};
        switch (vcal::dgp::simulate_panel(ps, n_defaults, factors)) {
            case vcal::dgp::Status::Ok: break;
            case vcal::dgp::Status::FpEnvironment:
                throw AbiError(VCAL_E_FP_ENVIRONMENT,
                               "flush-to-zero or denormals-are-zero is on in this thread; the DGP needs IEEE "
                               "subnormals to be bitwise reproducible (D-109)");
            case vcal::dgp::Status::FactorStreamExhausted:
                throw AbiError(VCAL_E_NUMERIC, "the factor draw exhausted its Philox block range (D-109)");
            default: throw std::logic_error("the DGP rejected a spec that the ABI accepted");
        }
    });
}

}  // extern "C"
