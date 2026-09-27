# SPDX-License-Identifier: Apache-2.0
#
# Enforces the layering table in ARCHITECTURE.md §2 (D-050).
# Usage: cmake -DROOT=<source tree> -P LayeringCheck.cmake
#
# Rules, applied to every #include in the named directory:
#   everywhere : no include path containing ".."
#   core/      : must not include engine/ backends/ resample/ dgp/ abi/ ref/
#                monitoring/ macro/ benchmark_models/ tests/
#   engine/    : must not include backends/ resample/ dgp/ abi/ ref/
#   dgp/       : must not include core/ engine/ backends/ abi/ ref/ (own math, D-031/Q14)
#   resample/  : must not include abi/ ref/
#   abi/       : must not include ref/ tests/ (the independent reference stays out of the library)
#   ref/       : may include only standard headers (no '/') and ref/... (D-031)

if(NOT DEFINED ROOT)
    message(FATAL_ERROR "LayeringCheck: pass -DROOT=<source tree>")
endif()

set(_forbid_core   "engine/;backends/;resample/;dgp/;abi/;ref/;monitoring/;macro/;benchmark_models/;tests/")
set(_forbid_engine "backends/;resample/;dgp/;abi/;ref/")
set(_forbid_dgp    "core/;engine/;backends/;abi/;ref/")
# resample/ (M2b) builds weight matrices and replicate estimates on top of engine/, with Philox
# from dgp/; it must not reach into the ABI or the independent reference.
set(_forbid_resample "abi/;ref/")
# abi/ (M2c) wraps everything above for C callers, but never the independent reference or tests.
set(_forbid_abi "ref/;tests/")

set(_violations 0)

function(_vcal_scan dir mode)
    if(NOT IS_DIRECTORY "${ROOT}/${dir}")
        return()
    endif()
    file(GLOB_RECURSE _files
        "${ROOT}/${dir}/*.hpp" "${ROOT}/${dir}/*.h"   "${ROOT}/${dir}/*.cpp"
        "${ROOT}/${dir}/*.cc"  "${ROOT}/${dir}/*.cu"  "${ROOT}/${dir}/*.cuh"
        "${ROOT}/${dir}/*.inl")
    set(_count ${_violations})
    foreach(_file IN LISTS _files)
        file(STRINGS "${_file}" _lines REGEX "^[ \t]*#[ \t]*include[ \t]*[<\"]")
        foreach(_line IN LISTS _lines)
            string(REGEX REPLACE "^[ \t]*#[ \t]*include[ \t]*[<\"]([^>\"]+)[>\"].*$" "\\1" _inc "${_line}")
            set(_bad "")
            if(_inc MATCHES "\\.\\.")
                set(_bad "relative '..' include")
            elseif(mode STREQUAL "ref")
                if(_inc MATCHES "/" AND NOT _inc MATCHES "^ref/")
                    set(_bad "ref/ may include only standard headers and ref/...")
                endif()
            else()
                foreach(_prefix IN LISTS _forbid_${mode})
                    string(FIND "${_inc}" "${_prefix}" _pos)
                    if(_pos EQUAL 0)
                        set(_bad "${dir}/ must not include ${_prefix}")
                    endif()
                endforeach()
            endif()
            if(_bad)
                file(RELATIVE_PATH _rel "${ROOT}" "${_file}")
                message(STATUS "layering violation: ${_rel}: #include \"${_inc}\" -- ${_bad}")
                math(EXPR _count "${_count} + 1")
            endif()
        endforeach()
    endforeach()
    set(_violations ${_count} PARENT_SCOPE)
endfunction()

_vcal_scan(core   core)
_vcal_scan(engine engine)
_vcal_scan(dgp    dgp)
_vcal_scan(resample resample)
_vcal_scan(abi    abi)
_vcal_scan(ref    ref)

if(_violations GREATER 0)
    message(FATAL_ERROR "found ${_violations} layering violation(s)")
endif()
message(STATUS "layering check passed")
