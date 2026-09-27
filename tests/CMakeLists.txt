# SPDX-License-Identifier: Apache-2.0

# In-house runner (D-029). Test sources include it as "tests/harness/vcal_test.hpp".
add_library(vcal_test_harness STATIC harness/runner.cpp harness/golden.cpp harness/sha256.cpp)
target_include_directories(vcal_test_harness PUBLIC "${PROJECT_SOURCE_DIR}")
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
vcal_add_test_executable(unit_resample SOURCES unit/resample_test.cpp LIBS vcal::resample vcal::backend_cpu)
vcal_add_test_executable(unit_ref SOURCES unit/ref_test.cpp LIBS vcal::ref)
vcal_add_test_executable(unit_dgp SOURCES unit/dgp_test.cpp LIBS vcal::dgp)

# --- M1.7 cross-reference (D-107, D-115) -----------------------------------------------------------
# One executable, three CTest entries. Estimates against ref take seconds per fit, so that entry is
# labelled "slow" and the Debug test preset skips it.
add_executable(unit_xref xref/xref_test.cpp)
target_link_libraries(unit_xref PRIVATE vcal_test_harness vcal::engine vcal::backend_cpu vcal::ref vcal::dgp)
vcal_target_defaults(unit_xref)
add_test(NAME xref_period_loglik COMMAND unit_xref xref_period_loglik_core_vs_ref)
add_test(NAME xref_core_regression COMMAND unit_xref xref_core_estimates_match_regression_values)
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
