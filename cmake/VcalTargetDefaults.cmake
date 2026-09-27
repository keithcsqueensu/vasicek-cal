# SPDX-License-Identifier: Apache-2.0
#
# vcal_target_defaults(<target> [CXX17] [FP_STRICT])
#   Applies the project-wide compile settings to a compiled target:
#   - C++20 (D-025), or C++17 with CXX17 for device-visible code (D-058)
#   - FP_STRICT: strict FP semantics for bitwise-reproducible code (dgp/, D-056): MSVC
#     /fp:strict instead of /fp:precise; GCC/Clang an explicit -fno-fast-math
#   - warning set, optionally as errors (D-049)
#   - FP contraction off (D-048): results must not depend on whether the
#     compiler chose to fuse a*b+c into an FMA. Fast-math is never enabled.
#   Flags are restricted to CXX so they never reach nvcc.

function(vcal_target_defaults target)
    cmake_parse_arguments(ARG "CXX17;FP_STRICT" "" "" ${ARGN})
    if(ARG_CXX17)
        set_target_properties(${target} PROPERTIES
            CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
    else()
        target_compile_features(${target} PRIVATE cxx_std_20)
    endif()

    if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        # MSVC does not contract under /fp:precise (contraction needs /fp:contract).
        if(ARG_FP_STRICT)
            set(_fp /fp:strict)
        else()
            set(_fp /fp:precise)
        endif()
        set(_opts /W4 /permissive- /utf-8 /Zc:__cplusplus /Zc:preprocessor ${_fp})
        if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            # clang-cl contracts within expressions by default.
            list(APPEND _opts /clang:-ffp-contract=off)
        endif()
        if(VCAL_WARNINGS_AS_ERRORS)
            list(APPEND _opts /WX)
        endif()
    else()
        set(_opts -Wall -Wextra -Wpedantic -Wconversion -Wshadow -ffp-contract=off)
        if(ARG_FP_STRICT)
            list(APPEND _opts -fno-fast-math)
        endif()
        if(VCAL_WARNINGS_AS_ERRORS)
            list(APPEND _opts -Werror)
        endif()
    endif()

    target_compile_options(${target} PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${_opts}>")
endfunction()

# vcal_cuda_target_defaults(<target>)
#   CUDA translation units: C++17 (D-058, D-094). A host-only function called from device
#   code, and any other nvcc warning, is an error (K-1, D-100). Device-side FP contraction
#   (--fmad) is decided with the CUDA backend in M4.
function(vcal_cuda_target_defaults target)
    set_target_properties(${target} PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON CUDA_EXTENSIONS OFF)
    target_compile_options(${target} PRIVATE
        "$<$<COMPILE_LANGUAGE:CUDA>:--Werror=cross-execution-space-call,all-warnings>")
endfunction()
