# SPDX-License-Identifier: Apache-2.0
#
# CUDA toolchain selection. No CUDA sources exist before M4; this module only
# settles the option and the architecture list so the decisions live in one place.
#
# D-030: VCAL_ENABLE_CUDA defaults ON only when a CUDA compiler is found, so a
#        CPU-only build needs no toolkit.
# D-026: sm_120 SASS only when nvcc >= 12.8; no other version-driven floor for it.
# D-040: SASS sm_89 (+ sm_120) plus PTX for compute_80, so A100/H100 can JIT.
# D-105: explicit architecture numbers only. D-106: Ninja generator only.
# D-058: CUDA TUs are C++17.
# D-065: floor nvcc 11.8 (first release with sm_89). nvcc 11.8 accepts only
#        GCC <= 11, Clang <= 14.0, MSVC 19.1x-19.3x (VS 2017 - VS 2022 17.9); nvcc
#        itself rejects other host compilers. Confirmed as the supported floor (D-094).

include(CheckLanguage)
check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
    set(_vcal_cuda_default ON)
else()
    set(_vcal_cuda_default OFF)
endif()
option(VCAL_ENABLE_CUDA "Build the CUDA backend (default: ON if a CUDA compiler is found)" ${_vcal_cuda_default})
set(VCAL_CUDA_ARCHITECTURES "" CACHE STRING
    "Override the CUDA architecture list (empty = automatic per D-026/D-040)")

if(NOT VCAL_ENABLE_CUDA)
    message(STATUS "vcal: CUDA backend disabled")
    return()
endif()

# D-106: CUDA builds use Ninja. The CUDA toolkit's Visual Studio (MSBuild) integration is tied
# to specific Visual Studio releases and may be absent for the one installed, so the Visual
# Studio generator is rejected rather than left to fail obscurely.
if(CMAKE_GENERATOR MATCHES "^Visual Studio")
    message(FATAL_ERROR "vcal: CUDA builds require the Ninja generator (D-106); "
        "configure from a developer shell with -G Ninja or use the cuda-release preset")
endif()

# D-105: architectures are explicit numbers (optionally -real/-virtual), never all, all-major
# or native. The build must not depend on which GPU, if any, the build machine has.
function(_vcal_check_cuda_archs what archs)
    foreach(_a IN LISTS archs)
        if(NOT _a MATCHES "^[0-9]+(-real|-virtual)?$")
            message(FATAL_ERROR "vcal: ${what} entry '${_a}' is not an explicit architecture number "
                "(e.g. 89-real); all, all-major and native are not allowed (D-105)")
        endif()
    endforeach()
endfunction()
if(DEFINED CACHE{CMAKE_CUDA_ARCHITECTURES})
    _vcal_check_cuda_archs(CMAKE_CUDA_ARCHITECTURES "$CACHE{CMAKE_CUDA_ARCHITECTURES}")
endif()
if(VCAL_CUDA_ARCHITECTURES)
    _vcal_check_cuda_archs(VCAL_CUDA_ARCHITECTURES "${VCAL_CUDA_ARCHITECTURES}")
endif()

if(NOT CMAKE_CUDA_COMPILER)
    message(FATAL_ERROR "vcal: VCAL_ENABLE_CUDA=ON but no CUDA compiler was found")
endif()

enable_language(CUDA)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_EXTENSIONS OFF)

if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 11.8)
    message(FATAL_ERROR
        "vcal: nvcc ${CMAKE_CUDA_COMPILER_VERSION} cannot target sm_89 (needs >= 11.8; D-065). "
        "Configure with -DVCAL_ENABLE_CUDA=OFF for a CPU-only build.")
endif()

if(VCAL_CUDA_ARCHITECTURES)
    set(_vcal_archs ${VCAL_CUDA_ARCHITECTURES})
else()
    set(_vcal_archs 89-real)
    if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 12.8)
        list(APPEND _vcal_archs 120-real)
    else()
        message(STATUS "vcal: nvcc ${CMAKE_CUDA_COMPILER_VERSION} < 12.8, sm_120 SASS not built (D-026)")
    endif()
    list(APPEND _vcal_archs 80-virtual)
endif()

set(CMAKE_CUDA_ARCHITECTURES "${_vcal_archs}")
message(STATUS "vcal: CUDA ${CMAKE_CUDA_COMPILER_VERSION}, architectures: ${_vcal_archs}")
