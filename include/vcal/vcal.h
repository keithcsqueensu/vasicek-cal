/* SPDX-License-Identifier: Apache-2.0
 *
 * vasicek-cal C ABI, version 0.1 (M2c; D-138..D-144). The only public header.
 *
 * Plain C99: fixed-width integers, doubles and pointers, nothing else. tests/abi compiles a
 * consumer as C99 and as C11 so that anything C++-only fails the build.
 *
 * Conventions
 * - Names: every function and type starts with vcal_, every macro and enumerator with VCAL_.
 *   Functions use the C calling convention (__cdecl on Windows).
 * - Versioning: vcal_abi_version() returns (major << 16) | minor. Adding a field or a function
 *   bumps the minor version; any change that could break an existing caller bumps the major.
 *   Callers should check that the library's major equals VCAL_ABI_VERSION_MAJOR and its minor is
 *   at least VCAL_ABI_VERSION_MINOR.
 * - Structs: every public struct starts with `uint32_t struct_size`, which the CALLER sets to
 *   sizeof(the struct) as the caller compiled it, for inputs and outputs alike. Later versions
 *   only append fields, so struct sizes only grow:
 *     - a struct_size below this version's size is VCAL_E_INVALID_ARGUMENT;
 *     - an input struct larger than the library knows is accepted if every byte past the known
 *       fields is zero (fields a newer caller left at their defaults), otherwise
 *       VCAL_E_UNSUPPORTED: the library never silently ignores a setting;
 *     - an output struct larger than the library knows has only its known fields written.
 *   Fields named `reserved` must be zero on input and are written as zero on output. Structs
 *   have no implicit padding (tests/abi checks the layout).
 * - Memory: the caller allocates every buffer. The library never returns memory for the
 *   caller to free, and keeps no caller pointer after a call returns. Variable-size outputs
 *   take (buffer, capacity, required): with a NULL buffer the call only stores the required
 *   length in *required and returns VCAL_OK; with a buffer smaller than that it stores the
 *   length and returns VCAL_E_BUFFER_TOO_SMALL without writing the buffer. Lengths count
 *   elements (bytes for strings, including the terminating NUL).
 * - Errors: every function returns a vcal_status, except vcal_abi_version and
 *   vcal_status_string, which cannot fail. A failing call stores a message for the calling
 *   thread, which vcal_last_error copies out. Every other function returning a vcal_status clears
 *   it on entry; vcal_abi_version, vcal_status_string and vcal_last_error leave it alone. No
 *   C++ exception crosses the boundary.
 * - Platforms: 64-bit only (the layouts are checked there).
 * - Threads: a context must not be used by two threads at the same time; distinct contexts
 *   are independent. Results are bitwise identical for any n_threads.
 * - Scope: one model in v0, the one-factor binomial-mixture MLE under the `parity` profile
 *   (docs/methodology/binomial_mixture_mle.md), on the CPU backend. Later objectives, profiles
 *   and backends arrive as appended fields or new functions.
 */
#ifndef VCAL_VCAL_H
#define VCAL_VCAL_H

#include <stdint.h>

#if defined(VCAL_STATIC)
#define VCAL_API
#elif defined(_WIN32)
#if defined(VCAL_BUILDING)
#define VCAL_API __declspec(dllexport)
#else
#define VCAL_API __declspec(dllimport)
#endif
#elif defined(__GNUC__)
#define VCAL_API __attribute__((visibility("default")))
#else
#define VCAL_API
#endif

#if defined(_WIN32)
#define VCAL_CALL __cdecl
#else
#define VCAL_CALL
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define VCAL_ABI_VERSION_MAJOR 0
#define VCAL_ABI_VERSION_MINOR 1
#define VCAL_ABI_VERSION ((VCAL_ABI_VERSION_MAJOR << 16) | VCAL_ABI_VERSION_MINOR)

/* ---- status codes ------------------------------------------------------------------------- */

typedef int32_t vcal_status;
enum {
    VCAL_OK = 0,
    VCAL_E_INVALID_ARGUMENT = 1,  /* a NULL pointer, a value out of range, an inconsistent panel */
    VCAL_E_BUFFER_TOO_SMALL = 2,  /* *required holds the length needed */
    VCAL_E_UNSUPPORTED = 3,       /* a nonzero field this library does not know (newer caller) */
    VCAL_E_NUMERIC = 4,           /* no finite log-likelihood anywhere on the grid, or similar */
    VCAL_E_OUT_OF_MEMORY = 5,
    VCAL_E_FP_ENVIRONMENT = 6,    /* flush-to-zero or denormals-are-zero is on (the DGP refuses) */
    VCAL_E_INTERNAL = 99          /* a bug: please report it with vcal_last_error's message */
};

/* ---- enumerations (int32_t fields) -------------------------------------------------------- */

enum { VCAL_PROFILE_PARITY = 0 }; /* textbook, scipy-replicable (D-035) */

enum { VCAL_SCALE_LINEAR = 0, VCAL_SCALE_LOG = 1, VCAL_SCALE_LOGIT = 2, VCAL_SCALE_PROBIT = 3 };

enum { VCAL_AXIS_PD = 0, VCAL_AXIS_RHO = 1 };

enum {
    VCAL_RESAMPLE_IID_BOOTSTRAP = 1,   /* periods drawn with replacement (D-132) */
    VCAL_RESAMPLE_BLOCK_BOOTSTRAP = 2, /* moving blocks, non-circular (D-133) */
    VCAL_RESAMPLE_JACKKNIFE = 3,       /* delete-one; T replicates */
    VCAL_RESAMPLE_WALK_FORWARD = 4,    /* rolling or expanding windows */
    VCAL_RESAMPLE_INDICES = 5,         /* a caller-supplied B x m period-index matrix (D-134) */
    VCAL_RESAMPLE_WEIGHTS = 6          /* a caller-supplied B x T weight matrix */
};

/* Estimate flags (uint32_t). */
enum {
    VCAL_FLAG_GRID_EDGE = 1u << 0,              /* argmax on the edge of the grid: no SEs */
    VCAL_FLAG_FLAT_SURFACE = 1u << 1,           /* curvature not negative definite: no SEs */
    VCAL_FLAG_QUADRATURE_UNCONVERGED = 1u << 2, /* the doubled-rule check failed somewhere (D-120) */
    VCAL_FLAG_REFINEMENT_REJECTED = 1u << 3,    /* quadratic step left the stencil: grid point kept */
    VCAL_FLAG_NUMERIC = 1u << 4,                /* NaN somewhere on the surface */
    VCAL_FLAG_NEAR_BOUND = 1u << 5              /* within 2 SEs of a grid bound: use profile intervals */
};

/* Interval flags (uint32_t), per parameter. */
enum {
    VCAL_INTERVAL_LOWER_TRUNCATED = 1u << 0, /* lower end is the grid bound (D-130) */
    VCAL_INTERVAL_UPPER_TRUNCATED = 1u << 1,
    VCAL_INTERVAL_NOT_COMPUTED = 1u << 2     /* no usable maximum (flat or numeric estimate) */
};

/* ---- library information and errors ------------------------------------------------------- */

VCAL_API uint32_t VCAL_CALL vcal_abi_version(void);

/* Build provenance as "key=value" lines: library and ABI version, git commit and whether the
 * tree was dirty, compiler, C++ standard, build type, OpenMP, CUDA, profiles. */
VCAL_API vcal_status VCAL_CALL vcal_build_info(char* buffer, int64_t capacity, int64_t* required);

/* The calling thread's message from its last failing call ("" after a successful one). Does
 * not change the stored message. */
VCAL_API vcal_status VCAL_CALL vcal_last_error(char* buffer, int64_t capacity, int64_t* required);

/* A static name for a status code, e.g. "VCAL_E_INVALID_ARGUMENT"; never NULL, never freed. */
VCAL_API const char* VCAL_CALL vcal_status_string(vcal_status status);

/* ---- context ------------------------------------------------------------------------------ */

typedef struct vcal_context vcal_context; /* opaque (D-045) */

typedef struct vcal_context_options {
    uint32_t struct_size;
    int32_t profile;   /* VCAL_PROFILE_PARITY */
    int32_t n_threads; /* CPU threads; 0 = the OpenMP default. Results do not depend on it. */
    uint32_t reserved;
} vcal_context_options;

/* options may be NULL for the defaults (parity, n_threads = 0). */
VCAL_API vcal_status VCAL_CALL vcal_context_create(const vcal_context_options* options, vcal_context** context);
/* NULL is accepted and ignored. */
VCAL_API vcal_status VCAL_CALL vcal_context_destroy(vcal_context* context);

/* ---- data and grid ------------------------------------------------------------------------ */

/* A default-count panel: period t has n_obligors[t] > 0 obligors and 0 <= n_defaults[t] <= n_obligors[t]. */
typedef struct vcal_panel {
    uint32_t struct_size;
    uint32_t reserved;
    int64_t n_periods;
    const int64_t* n_obligors; /* n_periods values, n <= 2^53 */
    const int64_t* n_defaults; /* n_periods values */
} vcal_panel;

/* The parameter grid: two axes, each uniform in its scaled coordinate (logit for the default),
 * with both ends included and at least 3 points. Grid point k = i_pd * rho_points + i_rho. */
typedef struct vcal_grid {
    uint32_t struct_size;
    uint32_t reserved;
    double pd_lo, pd_hi;
    double rho_lo, rho_hi;
    int32_t pd_points, pd_scale;   /* VCAL_SCALE_* */
    int32_t rho_points, rho_scale;
} vcal_grid;

/* Fills the grid used by the recovery study (D-115): PD [1e-4, 0.2] x 61 and rho [1e-3, 0.5]
 * x 41, both logit. grid->struct_size must be set. */
VCAL_API vcal_status VCAL_CALL vcal_grid_default(vcal_grid* grid);

/* The natural-scale values of one axis (VCAL_AXIS_*), exactly as the engine uses them. */
VCAL_API vcal_status VCAL_CALL vcal_grid_values(const vcal_grid* grid, int32_t axis, double* values, int64_t capacity,
                                                int64_t* required);

/* ---- calibration -------------------------------------------------------------------------- */

typedef struct vcal_estimate {
    uint32_t struct_size;
    uint32_t flags; /* VCAL_FLAG_* */
    double pd, rho;
    double se_pd, se_rho, corr_pd_rho; /* observed information at the estimate; NaN if unavailable */
    double loglik;                     /* log-likelihood at the estimate, including log C(n, d) */
    double quad_check_max;             /* max_t |l_t(rule) - l_t(doubled rule)| at the estimate */
    double quad_check_total;           /* sum_t of the same */
    int64_t quad_check_flagged;        /* periods above the D-120 threshold */
    int64_t grid_index;                /* grid argmax before refinement */
    int64_t nan_count;                 /* NaN surface values seen */
} vcal_estimate;

/* 95% profile-likelihood intervals (threshold chi2_1(0.95)/2 = 1.9207294103470630, truncated
 * at the grid bounds; D-128..D-130). The recommended intervals for this model (D-137). */
typedef struct vcal_profile_intervals {
    uint32_t struct_size;
    uint32_t pd_flags;  /* VCAL_INTERVAL_* */
    uint32_t rho_flags;
    uint32_t reserved;
    double pd_lo, pd_hi;
    double rho_lo, rho_hi;
    double loglik_max;          /* the maximum polished off the grid */
    double pd_at_max, rho_at_max;
    double residual_max;        /* max |profile(endpoint) - threshold| over solved endpoints */
    int64_t evaluations;        /* panel log-likelihood evaluations used */
} vcal_profile_intervals;

/* Fits the panel on the grid (ARCHITECTURE.md §5.1). profile may be NULL; if not, the profile
 * intervals are computed too (roughly as costly again as the fit). */
VCAL_API vcal_status VCAL_CALL vcal_calibrate(vcal_context* context, const vcal_panel* panel, const vcal_grid* grid,
                                              vcal_estimate* estimate, vcal_profile_intervals* profile);

/* The per-period log-likelihood surface L, n_periods x (pd_points * rho_points), row-major:
 * L[t * K + k] = l_t at grid point k. For cell-by-cell comparison with an independent script. */
VCAL_API vcal_status VCAL_CALL vcal_surface(vcal_context* context, const vcal_panel* panel, const vcal_grid* grid,
                                            double* surface, int64_t capacity, int64_t* required);

/* ---- resampling --------------------------------------------------------------------------- */

/* Fields a scheme does not use must be zero / NULL. */
typedef struct vcal_resample_spec {
    uint32_t struct_size;
    int32_t scheme;         /* VCAL_RESAMPLE_* */
    uint64_t seed;          /* bootstrap schemes: the stream is (seed, scheme, replicate, draw) */
    int64_t replicates;     /* bootstrap, INDICES, WEIGHTS: B >= 1 (bootstrap: B < 2^32) */
    int64_t block_length;   /* BLOCK_BOOTSTRAP: 1 <= l <= T; 0 = ceil(T^(1/3)) */
    int64_t window;         /* WALK_FORWARD: window length; 0 = expanding from period 0 */
    int64_t first_end;      /* WALK_FORWARD: first window ends before this period, 1 <= first_end <= T */
    int64_t step;           /* WALK_FORWARD: >= 1; windows end at first_end + b * step <= T */
    int64_t draws;          /* INDICES: indices per replicate */
    double level;           /* percentile interval level; 0 = 0.95 */
    const int64_t* indices; /* INDICES: replicates x draws, each in [0, T) */
    const double* weights;  /* WEIGHTS: replicates x T, finite and >= 0 */
} vcal_resample_spec;

/* Per-replicate outputs as caller-owned arrays of `capacity` elements; any may be NULL. */
typedef struct vcal_replicates {
    uint32_t struct_size;
    uint32_t reserved;
    int64_t capacity;
    double* pd;
    double* rho;
    double* grid_loglik;  /* the replicate's weighted log-likelihood at its grid argmax */
    int64_t* grid_index;
    uint32_t* flags;      /* VCAL_FLAG_* (edge, flat, rejected, numeric); no SEs per replicate */
} vcal_replicates;

/* Percentile intervals (Hyndman-Fan type 7; D-135). Provided for comparison: not recommended
 * for rho, where they undercover (D-137). */
typedef struct vcal_percentile_intervals {
    uint32_t struct_size;
    uint32_t reserved;
    double level;
    double pd_lo, pd_hi;
    double rho_lo, rho_hi;
    int64_t replicates;
    int64_t pd_excluded, rho_excluded; /* non-finite replicate values, left out */
    int64_t grid_edge;                 /* replicates on the grid edge, kept at their grid value */
} vcal_percentile_intervals;

/* Replicate estimates from W x L (D-134): the quadrature runs once, each replicate is a
 * reweighting. replicates and intervals may each be NULL; with both NULL the call is a size
 * query and stores the replicate count B in *required. A non-NULL replicates needs
 * capacity >= B. */
VCAL_API vcal_status VCAL_CALL vcal_resample(vcal_context* context, const vcal_panel* panel, const vcal_grid* grid,
                                             const vcal_resample_spec* spec, vcal_replicates* replicates,
                                             vcal_percentile_intervals* intervals, int64_t* required);

/* The B x n_periods weight matrix W the scheme implies (counts for the bootstrap schemes), so
 * that any replicate can be rebuilt elsewhere. */
VCAL_API vcal_status VCAL_CALL vcal_resample_weights(vcal_context* context, const vcal_resample_spec* spec,
                                                     int64_t n_periods, double* weights, int64_t capacity,
                                                     int64_t* required);

/* ---- synthetic data ----------------------------------------------------------------------- */

/* One panel of the one-factor DGP (docs/methodology/dgp.md): a pure function of
 * (seed, scenario, replicate, pd, rho, n_obligors), bitwise identical on every platform. */
typedef struct vcal_dgp_spec {
    uint32_t struct_size;
    uint32_t reserved;
    uint64_t seed;
    uint32_t scenario;
    uint32_t replicate;
    double pd;                 /* 0 < pd < 1 */
    double rho;                /* 0 < rho < 1 */
    int64_t n_periods;         /* 1 .. 2^32 */
    const int64_t* n_obligors; /* n_periods values, 0 <= n <= 2^33 - 2^17 */
} vcal_dgp_spec;

/* n_defaults (and factors, the Z_t, if not NULL) receive n_periods values. */
VCAL_API vcal_status VCAL_CALL vcal_dgp_simulate(vcal_context* context, const vcal_dgp_spec* spec, int64_t* n_defaults,
                                                 double* factors, int64_t capacity, int64_t* required);

#ifdef __cplusplus
}
#endif

#endif /* VCAL_VCAL_H */
