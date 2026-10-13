# SPDX-License-Identifier: Apache-2.0
#
# Link-time and profile-guided optimisation for local study builds (P-13, D-403; perf/p13_lto_pgo.md).
# Off by default: CI's presets, the shared library's release builds and every golden use neither. Both leave
# floating-point semantics alone (no contraction, no reassociation), so every result is bit for bit the plain
# build's; perf/p13_pgo.py checks that on the study tools it trains on.
#
#   -DVCAL_LTO=ON                     link-time optimisation (CMake's IPO: GCC/Clang -flto, MSVC /GL + /LTCG)
#   -DVCAL_PGO=GENERATE|USE           instrumented build, then the optimised build from its profiles
#   -DVCAL_PGO_DIR=<dir>              where the profiles live (default <build>/pgo)
#
# perf/p13_pgo.py runs the whole cycle: GENERATE, build, the training runs, (Clang: merge), USE, rebuild.
# CPU builds only: with VCAL_ENABLE_CUDA these options are refused.

option(VCAL_LTO "Link-time optimisation for local study builds (P-13, D-403)" OFF)
set(VCAL_PGO "OFF" CACHE STRING "Profile-guided optimisation for local study builds: OFF, GENERATE or USE (P-13, D-403)")
set_property(CACHE VCAL_PGO PROPERTY STRINGS OFF GENERATE USE)
set(VCAL_PGO_DIR "${PROJECT_BINARY_DIR}/pgo" CACHE PATH "Profile directory for VCAL_PGO (P-13, D-403)")

if(NOT VCAL_PGO MATCHES "^(OFF|GENERATE|USE)$")
    message(FATAL_ERROR "VCAL_PGO must be OFF, GENERATE or USE (got '${VCAL_PGO}')")
endif()
if((VCAL_LTO OR NOT VCAL_PGO STREQUAL "OFF") AND VCAL_ENABLE_CUDA)
    message(FATAL_ERROR "VCAL_LTO and VCAL_PGO are for CPU study builds only (P-13, D-403); configure with "
                        "VCAL_ENABLE_CUDA=OFF")
endif()

if(VCAL_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _vcal_ipo_ok OUTPUT _vcal_ipo_msg LANGUAGES CXX)
    if(NOT _vcal_ipo_ok)
        message(FATAL_ERROR "VCAL_LTO: this toolchain does not support link-time optimisation: ${_vcal_ipo_msg}")
    endif()
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    message(STATUS "vcal: link-time optimisation on (P-13)")
endif()

if(NOT VCAL_PGO STREQUAL "OFF")
    file(TO_CMAKE_PATH "${VCAL_PGO_DIR}" _vcal_pgo_dir)
    file(MAKE_DIRECTORY "${_vcal_pgo_dir}")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if(VCAL_PGO STREQUAL "GENERATE")
            # Atomic counter updates: the study tools are multithreaded, and racing counters corrupt the profile.
            add_compile_options(-fprofile-generate=${_vcal_pgo_dir} -fprofile-update=atomic)
            add_link_options(-fprofile-generate=${_vcal_pgo_dir})
        else()
            # Targets the training did not run have no profile: that is expected, not an error.
            add_compile_options(-fprofile-use=${_vcal_pgo_dir} -fprofile-correction -Wno-missing-profile)
            add_link_options(-fprofile-use=${_vcal_pgo_dir})
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        if(VCAL_PGO STREQUAL "GENERATE")
            add_compile_options(-fprofile-generate=${_vcal_pgo_dir})
            add_link_options(-fprofile-generate=${_vcal_pgo_dir})
        else()
            # p13_pgo.py merges the raw profiles into default.profdata first (llvm-profdata merge).
            add_compile_options(-fprofile-use=${_vcal_pgo_dir}/default.profdata -Wno-profile-instr-unprofiled
                                -Wno-profile-instr-out-of-date -Wno-backend-plugin)
            add_link_options(-fprofile-use=${_vcal_pgo_dir}/default.profdata)
        endif()
    elseif(MSVC)
        # MSVC's PGO works on link-time code generation: /GL to compile, /LTCG with a per-binary .pgd to link.
        # The instrumented binaries write <name>!N.pgc next to the .pgd; the USE link merges them.
        add_compile_options(/GL)
        if(VCAL_PGO STREQUAL "GENERATE")
            add_link_options(/LTCG "/GENPROFILE:PGD=${_vcal_pgo_dir}/$<TARGET_FILE_BASE_NAME:$<TARGET_PROPERTY:NAME>>.pgd")
        else()
            add_link_options(/LTCG "/USEPROFILE:PGD=${_vcal_pgo_dir}/$<TARGET_FILE_BASE_NAME:$<TARGET_PROPERTY:NAME>>.pgd")
        endif()
    else()
        message(FATAL_ERROR "VCAL_PGO: unsupported compiler '${CMAKE_CXX_COMPILER_ID}' (GCC, Clang or MSVC)")
    endif()
    message(STATUS "vcal: profile-guided optimisation ${VCAL_PGO}, profiles in ${_vcal_pgo_dir} (P-13)")
endif()
