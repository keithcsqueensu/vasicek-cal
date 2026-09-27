# SPDX-License-Identifier: Apache-2.0
#
# Writes OUTPUT with the git commit of SOURCE_DIR and whether tracked files differ from it, for
# vcal_build_info (D-142). Runs at every build, so the commit is never stale, and rewrites OUTPUT
# only when it changes, so an unchanged tree recompiles nothing.
# Usage: cmake -DSOURCE_DIR=<tree> -DOUTPUT=<header> -P GitInfo.cmake

cmake_minimum_required(VERSION 3.25)  # script mode: set policies explicitly

set(_commit "unknown")
set(_dirty "unknown")
find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${SOURCE_DIR}" rev-parse HEAD
        OUTPUT_VARIABLE _out RESULT_VARIABLE _rc OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(_rc EQUAL 0 AND _out MATCHES "^[0-9a-f]+$")
        set(_commit "${_out}")
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${SOURCE_DIR}" status --porcelain --untracked-files=no
            OUTPUT_VARIABLE _status RESULT_VARIABLE _rc ERROR_QUIET)
        if(_rc EQUAL 0)
            if(_status STREQUAL "")
                set(_dirty "false")
            else()
                set(_dirty "true")
            endif()
        endif()
    endif()
endif()

set(_content "// GENERATED at build time by cmake/GitInfo.cmake.
#define VCAL_GIT_COMMIT \"${_commit}\"
#define VCAL_GIT_DIRTY \"${_dirty}\"
")
set(_old "")
if(EXISTS "${OUTPUT}")
    file(READ "${OUTPUT}" _old)
endif()
if(NOT _old STREQUAL _content)
    file(WRITE "${OUTPUT}" "${_content}")
endif()
