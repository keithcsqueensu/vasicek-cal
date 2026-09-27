# SPDX-License-Identifier: Apache-2.0
#
# CI starts no run for prose-only commits (D-145): the workflow's push and pull_request `paths`
# filters exclude Markdown, then re-include the Markdown files that tests or build scripts read.
# This check derives that re-include list from the code and fails if the workflow's list differs,
# in either direction: a file read but not listed (edits to it would skip CI), or a file listed but
# no longer read.
#
# "Read" means referenced outside comments in what CI executes: tests/ (C, C++, CUDA, CMake, TOML
# and JSON), cmake/, tools/, validation/, the top-level CMakeLists.txt and CMakePresets.json.
# Comments (C/C++, CMake, Python and TOML) and Python docstrings are stripped first,
# so mentions for readers do not count. A reference counts only if it names an existing file under
# ROOT, which leaves out the links recovery_harness writes into its generated Markdown.
#
# Usage: cmake -DROOT=<source tree> -DWORKFLOW=<ci.yml> -P CiPathFilterCheck.cmake
cmake_minimum_required(VERSION 3.25)  # script mode: set policies explicitly

foreach(_v ROOT WORKFLOW)
    if(NOT DEFINED ${_v})
        message(FATAL_ERROR "CiPathFilterCheck: pass -D${_v}=...")
    endif()
endforeach()

# Removes every span from `open` to the next `close` (block comments, docstrings).
function(_strip_spans var open close)
    set(_s "${${var}}")
    set(_out "")
    string(LENGTH "${open}" _open_len)
    string(LENGTH "${close}" _close_len)
    while(TRUE)
        string(FIND "${_s}" "${open}" _a)
        if(_a EQUAL -1)
            break()
        endif()
        string(SUBSTRING "${_s}" 0 ${_a} _head)
        string(APPEND _out "${_head}")
        math(EXPR _from "${_a} + ${_open_len}")
        string(SUBSTRING "${_s}" ${_from} -1 _s)
        string(FIND "${_s}" "${close}" _b)
        if(_b EQUAL -1)
            set(_s "")
            break()
        endif()
        math(EXPR _from "${_b} + ${_close_len}")
        string(SUBSTRING "${_s}" ${_from} -1 _s)
    endwhile()
    set(${var} "${_out}${_s}" PARENT_SCOPE)
endfunction()

# --- what the code reads ----------------------------------------------------------------------
file(GLOB_RECURSE _sources
    "${ROOT}/tests/*.cpp" "${ROOT}/tests/*.c" "${ROOT}/tests/*.hpp" "${ROOT}/tests/*.h"
    "${ROOT}/tests/*.cu" "${ROOT}/tests/*.cuh" "${ROOT}/tests/*.toml" "${ROOT}/tests/*.json"
    "${ROOT}/tests/*.cmake" "${ROOT}/tests/CMakeLists.txt" "${ROOT}/cmake/*.cmake"
    "${ROOT}/tools/*.py" "${ROOT}/validation/*.py")
list(APPEND _sources "${ROOT}/CMakeLists.txt" "${ROOT}/CMakePresets.json")

set(_read "")
foreach(_file IN LISTS _sources)
    file(READ "${_file}" _text)
    if(_file MATCHES "\\.(cpp|c|hpp|h|cu|cuh)$")
        _strip_spans(_text "/*" "*/")
        string(REGEX REPLACE "//[^\n]*" "" _text "${_text}")
    elseif(_file MATCHES "\\.json$")
        # JSON has no comments: every reference counts
    elseif(_file MATCHES "\\.(py|toml)$")
        _strip_spans(_text "\"\"\"" "\"\"\"")
        _strip_spans(_text "'''" "'''")
        string(REGEX REPLACE "#[^\n]*" "" _text "${_text}")
    else()  # CMake
        string(REGEX REPLACE "#[^\n]*" "" _text "${_text}")
    endif()
    string(REGEX MATCHALL "[A-Za-z0-9_][A-Za-z0-9_./-]*\\.md" _refs "${_text}")
    foreach(_ref IN LISTS _refs)
        if(EXISTS "${ROOT}/${_ref}" AND NOT IS_DIRECTORY "${ROOT}/${_ref}")
            list(APPEND _read "${_ref}")
        endif()
    endforeach()
endforeach()
list(REMOVE_DUPLICATES _read)
list(SORT _read)

# --- what the workflow re-includes ------------------------------------------------------------
file(STRINGS "${WORKFLOW}" _lines)
set(_section "")
set(_push "")
set(_pr "")
foreach(_line IN LISTS _lines)
    if(_line MATCHES "^  (push|pull_request):")
        set(_section "${CMAKE_MATCH_1}")
    elseif(_line MATCHES "^ ? ?[^ ]")
        set(_section "")
    elseif(_section AND _line MATCHES "^ +- '([^'!*][^']*)'")
        if(_section STREQUAL "push")
            list(APPEND _push "${CMAKE_MATCH_1}")
        else()
            list(APPEND _pr "${CMAKE_MATCH_1}")
        endif()
    endif()
endforeach()
list(SORT _push)
list(SORT _pr)
if(NOT _push STREQUAL _pr)
    message(FATAL_ERROR "the push and pull_request path filters re-include different files\n"
                        "  push: ${_push}\n  pull_request: ${_pr}")
endif()

set(_missing ${_read})
if(_push)
    list(REMOVE_ITEM _missing ${_push})
endif()
set(_stale ${_push})
if(_read)
    list(REMOVE_ITEM _stale ${_read})
endif()
if(_missing OR _stale)
    message(FATAL_ERROR "the CI path filter does not match the Markdown files the code reads\n"
                        "  read but not re-included: ${_missing}\n"
                        "  re-included but not read by any test: ${_stale}")
endif()
list(LENGTH _read _n)
message(STATUS "CI path filter re-includes exactly the ${_n} Markdown files the code reads: ${_read}")
