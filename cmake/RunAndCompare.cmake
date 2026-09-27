# SPDX-License-Identifier: Apache-2.0
#
# Runs EXE and requires its standard output to equal the contents of EXPECTED (line endings
# normalised, trailing whitespace ignored). Used for documentation examples whose stated output
# must stay true.
# Usage: cmake -DEXE=<program> -DEXPECTED=<file> -P RunAndCompare.cmake
cmake_minimum_required(VERSION 3.25)  # script mode: set policies explicitly

execute_process(COMMAND "${EXE}" OUTPUT_VARIABLE _out ERROR_VARIABLE _err RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "${EXE} exited with ${_rc}\n${_out}${_err}")
endif()
file(READ "${EXPECTED}" _expected)
foreach(_v _out _expected)
    string(REPLACE "\r\n" "\n" ${_v} "${${_v}}")
    string(STRIP "${${_v}}" ${_v})
endforeach()
if(NOT _out STREQUAL _expected)
    message(FATAL_ERROR "output differs from the documented output\n--- documented\n${_expected}\n--- actual\n${_out}")
endif()
message(STATUS "output matches the documented output:\n${_out}")
