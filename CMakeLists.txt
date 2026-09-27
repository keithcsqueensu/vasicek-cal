# SPDX-License-Identifier: Apache-2.0
cmake_minimum_required(VERSION 3.25)  # D-047

project(vasicek_cal
    VERSION 0.0.1
    DESCRIPTION "Calibration engine for Vasicek/ASRF PD models"
    LANGUAGES C CXX)  # C: the ABI's pure-C consumer tests (D-144)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# D-025: C++20 by default, no GNU extensions. Device-visible code (core/, engine/,
# backends/cuda/) is C++17 (D-058, D-062); tests/cxx17 enforces it.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(VCAL_WARNINGS_AS_ERRORS "Treat compiler warnings as errors (ON in presets and CI, D-049)" OFF)
option(VCAL_BUILD_TESTS "Build the test suite" ${PROJECT_IS_TOP_LEVEL})

include(VcalTargetDefaults)
include(VcalCuda)  # defines VCAL_ENABLE_CUDA (D-030); architectures per D-026/D-040

# core/: header-only, host+device, C++17 (D-058). Consumers include it as "core/...".
add_library(vcal_core INTERFACE)
add_library(vcal::core ALIAS vcal_core)
target_include_directories(vcal_core INTERFACE "${PROJECT_SOURCE_DIR}")
target_compile_features(vcal_core INTERFACE cxx_std_17)

# engine/: backend-agnostic surface evaluation, reduction, refinement; C++17 zone (D-062).
add_library(vcal_engine INTERFACE)
add_library(vcal::engine ALIAS vcal_engine)
target_link_libraries(vcal_engine INTERFACE vcal_core)

# ref/: independent FP64 reference, host-only, standard library only (D-031); the layering
# check keeps it from including anything else in the repository (D-050).
add_library(vcal_ref STATIC ref/ref.cpp)
add_library(vcal::ref ALIAS vcal_ref)
target_include_directories(vcal_ref PUBLIC "${PROJECT_SOURCE_DIR}")
vcal_target_defaults(vcal_ref)

# dgp/: synthetic data, CPU only (D-034), bitwise identical across platforms (D-052..D-057):
# its own maths, strict FP semantics, no contraction; shares no code with core (D-050).
add_library(vcal_dgp STATIC dgp/det_math.cpp dgp/dgp.cpp)
add_library(vcal::dgp ALIAS vcal_dgp)
target_include_directories(vcal_dgp PUBLIC "${PROJECT_SOURCE_DIR}")
vcal_target_defaults(vcal_dgp FP_STRICT)
set_target_properties(vcal_dgp PROPERTIES POSITION_INDEPENDENT_CODE ON)  # linked into the shared vcal

# resample/: weight matrices W and replicate estimates from W x L (M2b, D-132..D-135). Header-only;
# Philox comes from dgp/philox.hpp (header-only, integer arithmetic).
add_library(vcal_resample INTERFACE)
add_library(vcal::resample ALIAS vcal_resample)
target_link_libraries(vcal_resample INTERFACE vcal_engine)

# backends/cpu/: OpenMP `parallel for` only (D-028). Without OpenMP it runs serially and gives
# identical results (determinism does not depend on the thread count, §6).
option(VCAL_ENABLE_OPENMP "Parallelise the CPU backend with OpenMP" ON)
add_library(vcal_backend_cpu INTERFACE)
add_library(vcal::backend_cpu ALIAS vcal_backend_cpu)
target_include_directories(vcal_backend_cpu INTERFACE "${PROJECT_SOURCE_DIR}")
if(VCAL_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(vcal_backend_cpu INTERFACE OpenMP::OpenMP_CXX)
        message(STATUS "vcal: CPU backend uses OpenMP ${OpenMP_CXX_VERSION}")
    else()
        message(WARNING "vcal: OpenMP not found; the CPU backend runs serially")
    endif()
endif()

# abi/: the C ABI v0 (M2c, D-138..D-144). A shared library whose only exports are the vcal_*
# functions declared in include/vcal/vcal.h, the one public header. Symbols are hidden by default;
# on Linux a version script also hides the standard library's template instantiations, and the
# abi_exports test compares the export table with the header.
add_library(vcal SHARED abi/vcal.cpp abi/build_info.cpp)
add_library(vcal::vcal ALIAS vcal)
target_include_directories(vcal
    PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>" "$<INSTALL_INTERFACE:include>"
    PRIVATE "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/generated")
target_link_libraries(vcal PRIVATE vcal_resample vcal_backend_cpu vcal_dgp)
target_compile_definitions(vcal PRIVATE VCAL_BUILDING VCAL_VERSION="${PROJECT_VERSION}" VCAL_BUILD_TYPE="$<CONFIG>")
set_target_properties(vcal PROPERTIES
    C_VISIBILITY_PRESET hidden CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON
    VERSION ${PROJECT_VERSION} SOVERSION 0)  # SOVERSION follows the ABI major version
vcal_target_defaults(vcal)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_options(vcal PRIVATE "LINKER:--version-script=${PROJECT_SOURCE_DIR}/abi/vcal.map" "LINKER:--no-undefined")
    set_property(TARGET vcal APPEND PROPERTY LINK_DEPENDS "${PROJECT_SOURCE_DIR}/abi/vcal.map")
endif()
# The git commit for vcal_build_info, refreshed at every build (cmake/GitInfo.cmake).
add_custom_target(vcal_git_info
    COMMAND "${CMAKE_COMMAND}" "-DSOURCE_DIR=${PROJECT_SOURCE_DIR}"
            "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/generated/vcal_git_info.h" -P "${PROJECT_SOURCE_DIR}/cmake/GitInfo.cmake"
    BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/generated/vcal_git_info.h"
    COMMENT "vcal: git commit for build info"
    VERBATIM)
add_dependencies(vcal vcal_git_info)

include(GNUInstallDirs)
install(TARGETS vcal
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES include/vcal/vcal.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/vcal)

if(VCAL_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
