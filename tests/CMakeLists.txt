# SPDX-License-Identifier: Apache-2.0

# In-house runner (D-029). Test sources include it as "tests/harness/vcal_test.hpp".
add_library(vcal_test_harness STATIC harness/runner.cpp harness/golden.cpp harness/sha256.cpp)
target_include_directories(vcal_test_harness PUBLIC "${PROJECT_SOURCE_DIR}")
# Every test can include core/, so every test links its one compiled part, the erfcx kernels (P-7, D-212).
target_link_libraries(vcal_test_harness PUBLIC vcal::core)
target_compile_definitions(vcal_test_harness PRIVATE VCAL_GOLDEN_DIR="${CMAKE_CURRENT_SOURCE_DIR}/golden")
vcal_target_defaults(vcal_test_harness)

# vcal_add_test_executable(<name> SOURCES ... [LIBS ...])
#   One executable per suite, registered as one CTest entry of the same name.
function(vcal_add_test_executable name)
    cmake_parse_arguments(ARG "" "" "SOURCES;LIBS" ${ARGN})
    add_executable(${name} ${ARG_SOURCES})
    target_link_libraries(${name} PRIVATE vcal_test_harness ${ARG_LIBS})
    vcal_target_defaults(${name})
    add_test(NAME ${name} COMMAND ${name})
endfunction()

# --- harness self-tests -------------------------------------------------------
vcal_add_test_executable(harness_selftest SOURCES harness/selftest.cpp)

# Deliberately failing checks: the runner must exit non-zero and count exactly
# one failed check per test (VCAL_REQUIRE must abort, not continue).
add_executable(harness_expected_failures harness/expected_failures.cpp)
target_link_libraries(harness_expected_failures PRIVATE vcal_test_harness)
vcal_target_defaults(harness_expected_failures)
add_test(NAME harness_expected_failures_exit_code COMMAND harness_expected_failures)
set_tests_properties(harness_expected_failures_exit_code PROPERTIES WILL_FAIL TRUE)
add_test(NAME harness_expected_failures_count COMMAND harness_expected_failures)
set_tests_properties(harness_expected_failures_count PROPERTIES
    PASS_REGULAR_EXPRESSION "7 of 7 tests failed \\(7 failed checks\\)")
add_test(NAME harness_unknown_test_name COMMAND harness_selftest no_such_test)
set_tests_properties(harness_unknown_test_name PROPERTIES
    PASS_REGULAR_EXPRESSION "unknown test: no_such_test")

# --- layering (D-050) ---------------------------------------------------------
add_test(NAME layering_check
    COMMAND ${CMAKE_COMMAND} -DROOT=${PROJECT_SOURCE_DIR} -P ${PROJECT_SOURCE_DIR}/cmake/LayeringCheck.cmake)
# The fixture tree holds one violation per rule; the check must find all of them.
add_test(NAME layering_check_fires
    COMMAND ${CMAKE_COMMAND} -DROOT=${CMAKE_CURRENT_SOURCE_DIR}/layering/fixture
            -P ${PROJECT_SOURCE_DIR}/cmake/LayeringCheck.cmake)
set_tests_properties(layering_check_fires PROPERTIES
    PASS_REGULAR_EXPRESSION "found 7 layering violation\\(s\\)")

# --- C++17 zone (D-058, D-062, D-064) -----------------------------------------
# Every header in core/, engine/ and backends/cuda/ must compile as C++17. Host
# compiler for now; the D-059 nvcc job will compile the same generated TU.
# (.cuh files need nvcc and join the check with that job.)
include(cxx17/check_headers.cmake)
file(GLOB_RECURSE _vcal_device_visible CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/core/*.hpp"
    "${PROJECT_SOURCE_DIR}/engine/*.hpp"
    "${PROJECT_SOURCE_DIR}/backends/cuda/*.hpp")
list(SORT _vcal_device_visible)
list(LENGTH _vcal_device_visible _vcal_n)
message(STATUS "vcal: C++17 header check covers ${_vcal_n} header(s)")
vcal_add_cxx17_header_check(cxx17_header_check
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/cxx17_header_check.cpp"
    HEADERS ${_vcal_device_visible})

# --- nvcc (K-1, D-059, D-100): device-visible code must compile as device code ------------------
if(VCAL_ENABLE_CUDA)
    file(GLOB_RECURSE _vcal_cuda_headers CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/backends/cuda/*.cuh")
    list(SORT _vcal_cuda_headers)
    vcal_add_cxx17_header_check(cxx17_header_check_nvcc LANGUAGE CUDA
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/cxx17_header_check_nvcc.cu"
        HEADERS ${_vcal_device_visible} ${_vcal_cuda_headers})
    add_library(cuda_device_compile_check OBJECT cuda/device_compile_check.cu)
    target_link_libraries(cuda_device_compile_check PRIVATE vcal::engine)
    vcal_cuda_target_defaults(cuda_device_compile_check)
    # ...and a host-only callee reached from device code must break the build.
    add_library(cuda_host_call_fixture OBJECT EXCLUDE_FROM_ALL cuda/fixture/host_call_from_device.cu)
    target_include_directories(cuda_host_call_fixture PRIVATE "${PROJECT_SOURCE_DIR}")
    vcal_cuda_target_defaults(cuda_host_call_fixture)
endif()

# Compile-fail tests build an EXCLUDE_FROM_ALL target and assert the diagnostic.
# They share the build tree, so they must not run concurrently.
function(vcal_add_compile_fail_test name target regex)
    add_test(NAME ${name} COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${target})
    set_tests_properties(${name} PROPERTIES PASS_REGULAR_EXPRESSION "${regex}" RESOURCE_LOCK build_tree)
endfunction()

vcal_add_cxx17_header_check(cxx17_header_check_fixture EXCLUDE_FROM_ALL
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/cxx17_header_check_fixture.cpp"
    HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/cxx17/fixture/uses_concepts.hpp")
vcal_add_compile_fail_test(cxx17_header_check_fires cxx17_header_check_fixture
    "uses_concepts\\.hpp[^\n]*error")

if(VCAL_ENABLE_CUDA)
    vcal_add_compile_fail_test(cuda_device_check_fires cuda_host_call_fixture
        "host_call_from_device\\.cu\\([0-9]+\\): error: calling a __host__ function[^\n]*is not allowed")
endif()

# --- policy contracts: static_assert messages must reach the user (D-063) ---------
add_library(contract_fail_precision OBJECT EXCLUDE_FROM_ALL contract_fail/precision_accum_float.cpp)
target_include_directories(contract_fail_precision PRIVATE "${PROJECT_SOURCE_DIR}")
vcal_target_defaults(contract_fail_precision CXX17)
vcal_add_compile_fail_test(contract_fail_precision_message contract_fail_precision
    "Precision::accum_t must be double")

# --- unit ---------------------------------------------------------------------
vcal_add_test_executable(unit_precision SOURCES unit/precision_test.cpp LIBS vcal::core)
vcal_add_test_executable(unit_special SOURCES unit/special_test.cpp LIBS vcal::core)
vcal_add_test_executable(unit_quadrature SOURCES unit/quadrature_test.cpp LIBS vcal::core)
vcal_add_test_executable(unit_grid SOURCES unit/grid_test.cpp LIBS vcal::core)
vcal_add_test_executable(unit_reducer SOURCES unit/reducer_test.cpp LIBS vcal::core)
vcal_add_test_executable(unit_engine SOURCES unit/engine_test.cpp LIBS vcal::engine vcal::backend_cpu)
vcal_add_test_executable(unit_profile SOURCES unit/profile_test.cpp LIBS vcal::engine vcal::backend_cpu)
# S-48: the HAC (Newey-West) sandwich SE, a native option (D-209).
vcal_add_test_executable(unit_hac SOURCES unit/hac_test.cpp LIBS vcal::engine vcal::backend_cpu)
vcal_add_test_executable(unit_vasicek_rate SOURCES unit/vasicek_rate_test.cpp LIBS vcal::engine vcal::backend_cpu vcal::dgp)
vcal_add_test_executable(unit_moments SOURCES unit/moments_test.cpp LIBS vcal::engine vcal::dgp)
vcal_add_test_executable(unit_posterior SOURCES unit/posterior_test.cpp LIBS vcal::engine vcal::backend_cpu vcal::dgp)
# M3: the grid-Bayesian estimator's reference results, for its regression test and scipy replication.
add_executable(posterior_reference recovery/posterior_reference.cpp)
target_link_libraries(posterior_reference PRIVATE vcal_test_harness vcal::engine vcal::backend_cpu vcal::dgp)
vcal_target_defaults(posterior_reference)
# M3: the method of moments' reference results, for its regression test and scipy replication.
add_executable(mom_reference recovery/mom_reference.cpp)
target_link_libraries(mom_reference PRIVATE vcal_test_harness vcal::engine vcal::dgp)
vcal_target_defaults(mom_reference)
# M3: the Vasicek-rate MLE's recovery on its own model; writes or replays (slow) its summary.
add_executable(rate_recovery_harness recovery/rate_recovery.cpp)
target_link_libraries(rate_recovery_harness PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(rate_recovery_harness)
add_test(NAME recovery_rate_replay COMMAND rate_recovery_harness)
set_tests_properties(recovery_rate_replay PROPERTIES LABELS slow)
vcal_add_test_executable(unit_surface_cache SOURCES unit/surface_cache_test.cpp LIBS vcal::engine vcal::backend_cpu vcal::resample vcal::dgp)
vcal_add_test_executable(unit_resample SOURCES unit/resample_test.cpp LIBS vcal::resample vcal::backend_cpu)
# S-4b, the Bartlett-corrected profile interval (D-182, D-184). One executable, two CTest entries: the
# reproduction of the S-4 run (995 bootstrap fits) is labelled slow.
add_executable(unit_bartlett unit/bartlett_test.cpp)
target_link_libraries(unit_bartlett PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(unit_bartlett)
add_test(NAME unit_bartlett COMMAND unit_bartlett bartlett_factor_is_the_mean_profile_statistic
         bartlett_is_deterministic_and_independent_of_the_thread_count
         bartlett_factor_on_the_bound_gives_the_profile_interval bartlett_rejects_an_invalid_spec)
add_test(NAME bartlett_s4_run COMMAND unit_bartlett bartlett_reproduces_the_s4_run)
set_tests_properties(bartlett_s4_run PROPERTIES LABELS slow)
vcal_add_test_executable(unit_ref SOURCES unit/ref_test.cpp LIBS vcal::ref)
# The AR(1)-factor likelihood (M6, D-195): phi = 0 bit for bit, ref's dense filter on adversarial panels, the
# placement check, the score. The smoothness sweep (ref at 41 points per line) is slow.
add_executable(unit_ar1_filter unit/ar1_filter_test.cpp)
target_link_libraries(unit_ar1_filter PRIVATE vcal_test_harness vcal::benchmark_models vcal::ref)
vcal_target_defaults(unit_ar1_filter)
add_test(NAME unit_ar1_filter
         COMMAND unit_ar1_filter ar1_phi_zero_is_the_per_period_likelihood_bit_for_bit
                 ar1_matches_the_dense_reference_on_adversarial_panels
                 ar1_doubled_nodes_cannot_see_placement_error_and_the_placement_check_does
                 ar1_score_matches_richardson_differences ar1_likelihood_goldens_reproduce
                 ar1_filter_is_the_same_under_every_kernel)
add_test(NAME unit_ar1_filter_smooth COMMAND unit_ar1_filter ar1_likelihood_is_smooth_in_phi_and_rho)
set_tests_properties(unit_ar1_filter_smooth PROPERTIES LABELS slow)
# P-18 (D-199): the filter suite under each transition kernel forced; a kernel this machine lacks skips the entry
# visibly (exit 77).
foreach(_k scalar sse2 avx2)
    add_test(NAME unit_ar1_filter_kernel_${_k}
             COMMAND unit_ar1_filter ar1_forced_kernel_is_available ar1_phi_zero_is_the_per_period_likelihood_bit_for_bit
                     ar1_matches_the_dense_reference_on_adversarial_panels ar1_score_matches_richardson_differences
                     ar1_likelihood_goldens_reproduce)
    set_tests_properties(unit_ar1_filter_kernel_${_k} PROPERTIES ENVIRONMENT "VCAL_AR1_KERNEL=${_k}" SKIP_RETURN_CODE 77)
endforeach()
# P-7 (D-212): the filter's committed likelihood goldens under each erfcx kernel forced.
foreach(_k scalar sse2 avx2)
    add_test(NAME unit_ar1_filter_erfcx_${_k}
             COMMAND unit_ar1_filter ar1_phi_zero_is_the_per_period_likelihood_bit_for_bit ar1_likelihood_goldens_reproduce)
    set_tests_properties(unit_ar1_filter_erfcx_${_k} PROPERTIES ENVIRONMENT "VCAL_ERFCX_KERNEL=${_k}")
endforeach()
# The AR(1) filter's transition kernels (P-18, D-199): every kernel against a plain reference, bit for bit; each
# kernel also forced by VCAL_AR1_KERNEL, a kernel this machine lacks reported as skipped (exit 77), never as passed.
add_executable(unit_ar1_transition unit/ar1_transition_test.cpp)
target_link_libraries(unit_ar1_transition PRIVATE vcal_test_harness vcal::benchmark_models)
vcal_target_defaults(unit_ar1_transition)
add_test(NAME unit_ar1_transition COMMAND unit_ar1_transition)
foreach(_k scalar sse2 avx2 none-such)
    add_test(NAME unit_ar1_kernel_${_k}
             COMMAND unit_ar1_transition ar1_kernel_selection_follows_vcal_ar1_kernel ar1_kernels_equal_the_reference_bit_for_bit)
    set_tests_properties(unit_ar1_kernel_${_k} PROPERTIES ENVIRONMENT "VCAL_AR1_KERNEL=${_k}" SKIP_RETURN_CODE 77)
endforeach()

# The erfcx lane kernels (P-7, D-212): every kernel against the scalar erfcx bit for bit, the integrand's node
# batches, and log_contrib with its derivatives under every kernel against no kernel; each kernel also forced by
# VCAL_ERFCX_KERNEL, a kernel this machine lacks reported as skipped (exit 77), never as passed. The fit suites with
# committed goldens run again under each forced kernel.
add_executable(unit_erfcx_batch unit/erfcx_batch_test.cpp)
target_link_libraries(unit_erfcx_batch PRIVATE vcal_test_harness vcal::core)
vcal_target_defaults(unit_erfcx_batch)
add_test(NAME unit_erfcx_batch COMMAND unit_erfcx_batch)
foreach(_k scalar sse2 avx2 none-such)
    add_test(NAME unit_erfcx_kernel_${_k} COMMAND unit_erfcx_batch erfcx_kernel_selection_follows_vcal_erfcx_kernel)
    set_tests_properties(unit_erfcx_kernel_${_k} PROPERTIES ENVIRONMENT "VCAL_ERFCX_KERNEL=${_k}" SKIP_RETURN_CODE 77)
endforeach()

# The AR(1) MLE's search and profile intervals (D-195, D-196): recovery on DGP panels with known phi, ref's
# independent surface at the estimate, the binomial MLE at phi = 0, the box, determinism, the profile ends.
add_executable(unit_ar1_search unit/ar1_search_test.cpp)
target_link_libraries(unit_ar1_search PRIVATE vcal_test_harness vcal::benchmark_models vcal::ref vcal::backend_cpu vcal::dgp)
vcal_target_defaults(unit_ar1_search)
add_test(NAME unit_ar1_search
         COMMAND unit_ar1_search ar1_fit_is_a_maximum_of_the_dense_reference ar1_search_at_phi_zero_is_the_binomial_mle
                 ar1_fit_holds_a_bound_and_flags_it ar1_fit_is_deterministic
                 ar1_profile_ends_are_where_the_profile_crosses_the_threshold)
# Slow (about 80 s): 39 cold-started 2-D maxima on each of four panels, and three T = 200 fits.
add_test(NAME unit_ar1_search_recovery
         COMMAND unit_ar1_search ar1_fit_is_not_beaten_on_a_fine_phi_grid
                 ar1_fit_recovers_parameters_on_long_informative_panels)
set_tests_properties(unit_ar1_search_recovery PROPERTIES LABELS slow)
vcal_add_test_executable(unit_dgp SOURCES unit/dgp_test.cpp LIBS vcal::dgp)
# The Bernoulli kernel forced at start-up through VCAL_DGP_KERNEL (P-17, D-194), and once with an invalid value:
# the selection must follow the variable, and an unavailable one must refuse. Only the tests that choose their
# kernel explicitly run here (the others would refuse on a machine lacking the forced kernel); the plain
# unit_dgp entry already runs the reference panels and the count comparison under every available kernel.
# A kernel the machine lacks (AVX2 on an older CPU, SSE2 off x86-64) makes its entry exit 77: reported as
# skipped, never as passed.
foreach(_k scalar sse2 avx2 none-such)
    add_test(NAME unit_dgp_kernel_${_k}
             COMMAND unit_dgp bernoulli_kernel_selection_follows_vcal_dgp_kernel
                     reference_panels_are_bitwise_identical_to_the_mirror
                     ar1_reference_panels_are_bitwise_identical_to_the_mirror)
    set_tests_properties(unit_dgp_kernel_${_k} PROPERTIES ENVIRONMENT "VCAL_DGP_KERNEL=${_k}" SKIP_RETURN_CODE 77)
endforeach()

# --- M1.7 cross-reference (D-107, D-115) -----------------------------------------------------------
# One executable, three CTest entries. Estimates against ref take seconds per fit, so that entry is
# labelled "slow" and the Debug test preset skips it.
add_executable(unit_xref xref/xref_test.cpp)
target_link_libraries(unit_xref PRIVATE vcal_test_harness vcal::engine vcal::backend_cpu vcal::ref vcal::dgp)
vcal_target_defaults(unit_xref)
add_test(NAME xref_period_loglik COMMAND unit_xref xref_period_loglik_core_vs_ref)
add_test(NAME xref_core_regression COMMAND unit_xref xref_core_estimates_match_regression_values)
# P-7 (D-212): the committed regression estimates under each erfcx kernel forced (skipped where unavailable).
foreach(_k scalar sse2 avx2)
    add_test(NAME xref_core_regression_erfcx_${_k} COMMAND unit_xref xref_core_estimates_match_regression_values)
    set_tests_properties(xref_core_regression_erfcx_${_k} PROPERTIES ENVIRONMENT "VCAL_ERFCX_KERNEL=${_k}")
endforeach()
add_test(NAME xref_surface_smoothness COMMAND unit_xref xref_surface_is_smooth_for_zero_default_panels)
add_test(NAME validation_goldens COMMAND unit_xref validation_goldens_are_current)
add_test(NAME xref_estimates COMMAND unit_xref xref_estimates_core_vs_ref xref_estimates_agree_on_the_bound)
add_test(NAME xref_profile COMMAND unit_xref xref_profile_intervals_core_vs_ref)
set_tests_properties(xref_estimates xref_profile PROPERTIES LABELS slow)

# --- M1.8 recovery harness (D-121) -------------------------------------------------------------------
# recovery_harness is a tool, not a test: `recovery_harness --write` refits R = 1000 replicates of 81
# scenarios (about half an hour on 24 threads) and rewrites tests/golden/recovery/ and
# docs/methodology/recovery_results.md. unit_recovery replays a subset and enforces the exit criteria.
add_executable(recovery_harness recovery/recovery_harness.cpp)
target_link_libraries(recovery_harness PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
target_compile_definitions(recovery_harness PRIVATE VCAL_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
vcal_target_defaults(recovery_harness)
add_executable(unit_recovery recovery/recovery_test.cpp)
target_link_libraries(unit_recovery PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
target_compile_definitions(unit_recovery PRIVATE VCAL_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
vcal_target_defaults(unit_recovery)
add_test(NAME recovery_statistics COMMAND unit_recovery recovery_summary_statistics_on_a_constructed_sample)
add_test(NAME recovery_exit_criteria COMMAND unit_recovery recovery_goldens_meet_the_exit_criteria)
add_test(NAME recovery_replay COMMAND unit_recovery recovery_replay)
set_tests_properties(recovery_replay PROPERTIES LABELS slow)

# --- studies (D-148) ---------------------------------------------------------------------------------
# S-1 (D-300): the calibration is invariant under z -> -z within rounding. A test-only mirrored
# objective against parity on 34 fixed panels; a few minutes, so labelled slow.
add_executable(study_z_sign_invariance studies/z_sign_invariance.cpp)
target_link_libraries(study_z_sign_invariance PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(study_z_sign_invariance)
add_test(NAME study_z_sign_invariance COMMAND study_z_sign_invariance)
set_tests_properties(study_z_sign_invariance PROPERTIES LABELS slow)

# The shared jackknife run (S-3, S-5, S-21; D-155): a study tool, with a smoke test on two
# replicates of one scenario (the pinning run is the tool itself, at R = 1,000).
add_executable(study_jackknife_run studies/jackknife_run.cpp)
target_link_libraries(study_jackknife_run PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(study_jackknife_run)
add_test(NAME study_jackknife_run_smoke COMMAND study_jackknife_run --replicates 2 --scenarios 29)
set_tests_properties(study_jackknife_run_smoke PROPERTIES LABELS slow)
# Spike (not merged without review): caching surface rows by (n, d) across a scenario's replicates.
add_executable(surface_cache_spike studies/surface_cache_spike.cpp)
target_link_libraries(surface_cache_spike PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(surface_cache_spike)
# The bounded argmax against the full grid (D-172): every bootstrap and jackknife row of recovery
# panels, bit for bit. The full verification is R = 1000 on the study subset; this slow test
# re-checks a reduced R on every CI run.
add_executable(study_argmax_check studies/argmax_check.cpp)
target_link_libraries(study_argmax_check PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(study_argmax_check)
add_test(NAME study_argmax_check COMMAND study_argmax_check --replicates 40)
# P-9's measurement (perf/p9_repeat_rate.md): repeated panels in the recovery and S-10 runs. Built so it
# keeps compiling; run by hand.
add_executable(p9_repeat_rate ${PROJECT_SOURCE_DIR}/perf/p9_repeat_rate.cpp)
target_link_libraries(p9_repeat_rate PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(p9_repeat_rate)
# P-17's benchmark (perf/p17_simd_philox.md): panel simulation with each Bernoulli kernel, every panel
# checked against simulate_panel. Run by hand; a short run is a test, so the check runs in CI.
add_executable(p17_simd_philox ${PROJECT_SOURCE_DIR}/perf/p17_simd_philox.cpp)
target_link_libraries(p17_simd_philox PRIVATE vcal::dgp)
vcal_target_defaults(p17_simd_philox)
add_test(NAME p17_simd_philox_identity COMMAND p17_simd_philox --seconds 0.04)
# P-7's benchmark (perf/p7_erfcx.md, D-212): log_contrib and its derivatives under each erfcx kernel, every kernel's
# values checked against the no-kernel path bit for bit first. A short run is a CTest entry.
add_executable(p7_erfcx ${PROJECT_SOURCE_DIR}/perf/p7_erfcx.cpp)
target_link_libraries(p7_erfcx PRIVATE vcal::core)
vcal_target_defaults(p7_erfcx)
add_test(NAME p7_erfcx_identity COMMAND p7_erfcx --seconds 0.04)
set_tests_properties(study_argmax_check PROPERTIES LABELS slow)
# Its committed summary against the reviewed labels, kept apart (D-155): fast, always run.
vcal_add_test_executable(unit_study_jackknife SOURCES studies/jackknife_reviewed_test.cpp LIBS vcal::engine)

# S-34 (D-160): a study tool, with a smoke test on two replicates of one scenario.
add_executable(study_severe_period studies/severe_period.cpp)
target_link_libraries(study_severe_period PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(study_severe_period)
add_test(NAME study_severe_period_smoke COMMAND study_severe_period --replicates 2 --scenarios 29)
set_tests_properties(study_severe_period_smoke PROPERTIES LABELS slow)

# The estimator-comparison pass (S-8, S-27, S-28; studies/mle-vs-mom/PREDICTION.md): every panel fitted by
# the binomial MLE, the method of moments and the Vasicek-rate MLE under each zero-default treatment.
add_executable(study_estimator_pass studies/estimator_pass.cpp)
target_link_libraries(study_estimator_pass PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(study_estimator_pass)
add_test(NAME study_estimator_pass_smoke COMMAND study_estimator_pass --replicates 2 --scenarios 29,4,90)
set_tests_properties(study_estimator_pass_smoke PROPERTIES LABELS slow)

# S-9 (studies/bayes-coverage/PREDICTION.md): the grid posterior on every recovery panel, four arms.
add_executable(study_bayes_coverage studies/bayes_coverage.cpp)
target_link_libraries(study_bayes_coverage PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(study_bayes_coverage)
add_test(NAME study_bayes_coverage_smoke COMMAND study_bayes_coverage --replicates 2 --scenarios 0,29,43)
set_tests_properties(study_bayes_coverage_smoke PROPERTIES LABELS slow)

# S-15 (studies/bayes-sbc/PREDICTION.md): simulation-based calibration of the grid posterior.
add_executable(study_bayes_sbc studies/bayes_sbc.cpp)
target_link_libraries(study_bayes_sbc PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(study_bayes_sbc)
add_test(NAME study_bayes_sbc_smoke COMMAND study_bayes_sbc --draws 3)
set_tests_properties(study_bayes_sbc_smoke PROPERTIES LABELS slow)

# S-10 and S-4 (studies/parametric-bootstrap, studies/bartlett-profile): the parametric bootstrap, the
# profile LR statistic at the truth and S-4b's factor from the same bootstrap panels.
# S-35 (D-190): the recovery matrix under the AR(1)-factor DGP; scored by
# studies/lrpd-serial-correlation/compare.py. The smoke test runs two replicates of one scenario.
# S-35's M6 arm (PREDICTION.md, the 2026-10-10 addendum; D-196): the AR(1)-aware MLE on S-35's panels, and with
# --seed the validation run (D-197). A tool; the smoke test runs two replicates of one subset scenario on its own
# seed, with a 3-panel comparator.
add_executable(study_lrpd_serial_m6 studies/lrpd_serial_m6.cpp)
target_link_libraries(study_lrpd_serial_m6 PRIVATE vcal_test_harness vcal::benchmark_models vcal::resample vcal::backend_cpu
                                                   vcal::dgp)
vcal_target_defaults(study_lrpd_serial_m6)
target_include_directories(study_lrpd_serial_m6 PRIVATE "${PROJECT_BINARY_DIR}")  # the checkpoint's source commit
add_dependencies(study_lrpd_serial_m6 vcal_git_info)
add_test(NAME study_lrpd_serial_m6_smoke
         COMMAND study_lrpd_serial_m6 --scenarios 72 --replicates 2 --phis 0,0.6 --seed 0x5AFE --comparator --boot-panels 3
                 --out m6_smoke.csv)
add_test(NAME study_lrpd_serial_m6_resume
         COMMAND "${CMAKE_COMMAND}" "-DTOOL=$<TARGET_FILE:study_lrpd_serial_m6>"
                 "-DWORK=${CMAKE_CURRENT_BINARY_DIR}/m6_resume_test" -P "${CMAKE_CURRENT_SOURCE_DIR}/studies/m6_resume_test.cmake")
set_tests_properties(study_lrpd_serial_m6_resume PROPERTIES LABELS slow)

# The AR(1) likelihood goldens (tests/golden/ar1/likelihood.csv; D-195): engine and ref at fixed parameters on
# hand-made panels, no fits. Run by hand to regenerate; unit_ar1_filter checks the engine against them.
add_executable(ar1_goldens studies/ar1_goldens.cpp)
target_link_libraries(ar1_goldens PRIVATE vcal::benchmark_models vcal::ref)
vcal_target_defaults(ar1_goldens)
# S-35's M6 arm, model-based reference numbers (D-196): the expected Fisher information of the AR(1) likelihood at
# the true parameters (scores only, no fit). Run by hand; a two-panel run is a smoke test.
add_executable(study_ar1_information studies/ar1_information.cpp)
target_link_libraries(study_ar1_information PRIVATE vcal_test_harness vcal::benchmark_models vcal::backend_cpu vcal::dgp)
vcal_target_defaults(study_ar1_information)
add_test(NAME study_ar1_information_smoke COMMAND study_ar1_information --panels 2 --out ar1_information_smoke.csv)
add_executable(study_lrpd_serial studies/lrpd_serial.cpp)
target_link_libraries(study_lrpd_serial PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(study_lrpd_serial)
add_test(NAME study_lrpd_serial_smoke COMMAND study_lrpd_serial --replicates 2 --scenarios 29,2 --phis 0,0.6
         --out lrpd_smoke.csv --rate-out lrpd_smoke_rate.csv)
# S-48 (D-209): S-35's panels refitted, with the HAC sandwich SE at three lags. Run by hand; a two-panel
# run is a smoke test.
add_executable(study_hac_se studies/hac_se.cpp)
target_link_libraries(study_hac_se PRIVATE vcal_test_harness vcal::backend_cpu vcal::dgp)
vcal_target_defaults(study_hac_se)
add_test(NAME study_hac_se_smoke COMMAND study_hac_se --replicates 2 --scenarios 29,2 --phis 0,0.6
         --out hac_se_smoke.csv)
add_executable(study_parametric_bootstrap studies/parametric_bootstrap.cpp)
target_link_libraries(study_parametric_bootstrap PRIVATE vcal_test_harness vcal::resample vcal::backend_cpu vcal::dgp)
vcal_target_defaults(study_parametric_bootstrap)
add_test(NAME study_parametric_bootstrap_smoke COMMAND study_parametric_bootstrap --replicates 2 --scenarios 29)
set_tests_properties(study_parametric_bootstrap_smoke PROPERTIES LABELS slow)
# P-11 (D-179): the checkpoint records the source commit (the generated build-info header).
target_include_directories(study_parametric_bootstrap PRIVATE "${PROJECT_BINARY_DIR}")
add_dependencies(study_parametric_bootstrap vcal_git_info)
add_test(NAME study_parametric_bootstrap_resume
         COMMAND "${CMAKE_COMMAND}" "-DTOOL=$<TARGET_FILE:study_parametric_bootstrap>"
                 "-DWORK=${CMAKE_CURRENT_BINARY_DIR}/resume_test" -P "${CMAKE_CURRENT_SOURCE_DIR}/studies/resume_test.cmake")
set_tests_properties(study_parametric_bootstrap_resume PROPERTIES LABELS slow)

# --- M2c C ABI (D-138..D-144) -----------------------------------------------------------------------
# The consumer is plain C: abi_c_test is C11, and abi_c99_header_check compiles vcal.h as C99 with
# pedantic errors (C11 on MSVC, which has no C99 mode), so anything C++-only in the header breaks
# the build. abi_engine_equivalence checks the numbers against direct engine calls, bit for bit.
function(vcal_c_target_defaults target standard)
    set_target_properties(${target} PROPERTIES C_STANDARD ${standard} C_STANDARD_REQUIRED ON C_EXTENSIONS OFF)
    if(CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        set(_opts /W4 /utf-8)
        if(VCAL_WARNINGS_AS_ERRORS)
            list(APPEND _opts /WX)
        endif()
    else()
        set(_opts -Wall -Wextra -Wpedantic -Wconversion -Wshadow -pedantic-errors)
        if(VCAL_WARNINGS_AS_ERRORS)
            list(APPEND _opts -Werror)
        endif()
    endif()
    target_compile_options(${target} PRIVATE "$<$<COMPILE_LANGUAGE:C>:${_opts}>")
endfunction()

add_executable(abi_c_test abi/abi_c_test.c)
target_link_libraries(abi_c_test PRIVATE vcal::vcal)
vcal_c_target_defaults(abi_c_test 11)
# C cannot include tests/tolerances.hpp, so the one tolerance abi_c_test uses is read from it and
# passed in: the register stays the single source of the value.
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/tolerances.hpp" _tol_line REGEX "TOL_RESAMPLE_VS_REFIT_REL = ")
string(REGEX REPLACE ".*TOL_RESAMPLE_VS_REFIT_REL = ([0-9.eE+-]+).*" "\\1" _tol_value "${_tol_line}")
if(NOT _tol_value MATCHES "^[0-9.eE+-]+$")
    message(FATAL_ERROR "abi_c_test: could not read TOL_RESAMPLE_VS_REFIT_REL from tolerances.hpp")
endif()
target_compile_definitions(abi_c_test PRIVATE "TOL_RESAMPLE_VS_REFIT_REL=${_tol_value}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tolerances.hpp")
add_test(NAME abi_c COMMAND abi_c_test "${CMAKE_CURRENT_SOURCE_DIR}/golden/dgp/reference_panels.csv")

add_library(abi_c99_header_check OBJECT abi/c99_header_check.c)
target_include_directories(abi_c99_header_check PRIVATE "${PROJECT_SOURCE_DIR}/include")
if(CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    vcal_c_target_defaults(abi_c99_header_check 11)
else()
    vcal_c_target_defaults(abi_c99_header_check 99)
endif()

vcal_add_test_executable(abi_engine_equivalence SOURCES abi/abi_equivalence_test.cpp
    LIBS vcal::vcal vcal::resample vcal::backend_cpu vcal::dgp)

# The README's C example is code, not prose: its first ```c block is compiled as C11 and run, and
# its output must equal the ```text block that follows it. Editing README.md re-runs configure.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/README.md")
file(READ "${PROJECT_SOURCE_DIR}/README.md" _vcal_readme)
# The blocks hold semicolons, so they are kept in named variables, never in a CMake list.
foreach(_kind c text)
    string(FIND "${_vcal_readme}" "```${_kind}\n" _begin)
    if(_begin EQUAL -1)
        message(FATAL_ERROR "README.md: no ```${_kind} block for the readme_c_example test")
    endif()
    string(LENGTH "```${_kind}\n" _len)
    math(EXPR _begin "${_begin} + ${_len}")
    string(SUBSTRING "${_vcal_readme}" ${_begin} -1 _vcal_readme)
    string(FIND "${_vcal_readme}" "```" _end)
    string(SUBSTRING "${_vcal_readme}" 0 ${_end} _vcal_readme_${_kind})
    string(SUBSTRING "${_vcal_readme}" ${_end} -1 _vcal_readme)
endforeach()
file(CONFIGURE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/readme_example.c" CONTENT "${_vcal_readme_c}" @ONLY)
file(CONFIGURE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/readme_example.txt" CONTENT "${_vcal_readme_text}" @ONLY)
add_executable(readme_c_example "${CMAKE_CURRENT_BINARY_DIR}/readme_example.c")
target_link_libraries(readme_c_example PRIVATE vcal::vcal)
vcal_c_target_defaults(readme_c_example 11)
add_test(NAME readme_c_example COMMAND "${CMAKE_COMMAND}" "-DEXE=$<TARGET_FILE:readme_c_example>"
    "-DEXPECTED=${CMAKE_CURRENT_BINARY_DIR}/readme_example.txt" -P "${PROJECT_SOURCE_DIR}/cmake/RunAndCompare.cmake")

# Windows finds DLLs next to the executable: copy vcal.dll there.
if(WIN32)
    foreach(_t abi_c_test abi_engine_equivalence readme_c_example)
        add_custom_command(TARGET ${_t} POST_BUILD
            COMMAND "${CMAKE_COMMAND}" -E copy_if_different "$<TARGET_FILE:vcal>" "$<TARGET_FILE_DIR:${_t}>"
            VERBATIM)
    endforeach()
endif()

# The export table is exactly the functions vcal.h declares (cmake/AbiSymbolCheck.cmake), and the
# check fires: against a fixture header with one function the library lacks, it must name that
# function as missing and the real ones as extra.
set(_vcal_symbol_tool "")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_NM)
    set(_vcal_symbol_tool "${CMAKE_NM}")
    set(_vcal_symbol_kind nm)
elseif(MSVC)
    get_filename_component(_vcal_msvc_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
    find_program(VCAL_DUMPBIN dumpbin HINTS "${_vcal_msvc_bin}")
    if(VCAL_DUMPBIN)
        set(_vcal_symbol_tool "${VCAL_DUMPBIN}")
        set(_vcal_symbol_kind dumpbin)
    endif()
endif()
if(_vcal_symbol_tool)
    foreach(_case "abi_exports;${PROJECT_SOURCE_DIR}/include/vcal/vcal.h"
                  "abi_exports_fires;${CMAKE_CURRENT_SOURCE_DIR}/abi/fixture/exports_fixture.h")
        list(GET _case 0 _name)
        list(GET _case 1 _header)
        add_test(NAME ${_name} COMMAND "${CMAKE_COMMAND}" "-DLIBRARY=$<TARGET_FILE:vcal>" "-DHEADER=${_header}"
            "-DTOOL=${_vcal_symbol_tool}" "-DKIND=${_vcal_symbol_kind}"
            -P "${PROJECT_SOURCE_DIR}/cmake/AbiSymbolCheck.cmake")
    endforeach()
    set_tests_properties(abi_exports_fires PROPERTIES
        PASS_REGULAR_EXPRESSION "missing: vcal_not_exported[
]+ *extra: vcal_build_info;vcal_calibrate;")
endif()
if(NOT TEST abi_exports)
    message(WARNING "vcal: no nm/dumpbin found; the abi_exports test is not registered")
endif()

# --- CI prose filter (D-145, D-147) --------------------------------------------------------------------
# Prose-only changes build nothing, so every Markdown file the code reads must be in the plan job's
# code_md list. ci_path_filter_sync derives that list from the code; _fires runs it against a workflow
# that omits README.md and lists a file nothing reads, and must report both; _fires_paths against one
# with a workflow-level path filter, which would leave the required ci-ok check unreported.
add_test(NAME ci_path_filter_sync
    COMMAND "${CMAKE_COMMAND}" "-DROOT=${PROJECT_SOURCE_DIR}" "-DWORKFLOW=${PROJECT_SOURCE_DIR}/.github/workflows/ci.yml"
            -P "${PROJECT_SOURCE_DIR}/cmake/CiPathFilterCheck.cmake")
add_test(NAME ci_path_filter_sync_fires
    COMMAND "${CMAKE_COMMAND}" "-DROOT=${PROJECT_SOURCE_DIR}" "-DWORKFLOW=${CMAKE_CURRENT_SOURCE_DIR}/ci_paths/fixture_ci.yml"
            -P "${PROJECT_SOURCE_DIR}/cmake/CiPathFilterCheck.cmake")
set_tests_properties(ci_path_filter_sync_fires PROPERTIES
    PASS_REGULAR_EXPRESSION "read but not re-included: README\\.md[\r\n]+ *re-included but not read by any test: STATE\\.md")
add_test(NAME ci_path_filter_sync_fires_paths
    COMMAND "${CMAKE_COMMAND}" "-DROOT=${PROJECT_SOURCE_DIR}" "-DWORKFLOW=${CMAKE_CURRENT_SOURCE_DIR}/ci_paths/fixture_paths_ci.yml"
            -P "${PROJECT_SOURCE_DIR}/cmake/CiPathFilterCheck.cmake")
set_tests_properties(ci_path_filter_sync_fires_paths PROPERTIES
    PASS_REGULAR_EXPRESSION "workflow-level path filter")

# --- tolerance register (ARCHITECTURE.md §7) -------------------------------------
add_test(NAME tolerance_register_sync
    COMMAND ${CMAKE_COMMAND} -DHEADER=${CMAKE_CURRENT_SOURCE_DIR}/tolerances.hpp
            -DDOC=${PROJECT_SOURCE_DIR}/docs/methodology/tolerances.md
            -P ${PROJECT_SOURCE_DIR}/cmake/ToleranceRegisterCheck.cmake)
