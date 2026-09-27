# SPDX-License-Identifier: Apache-2.0
#
# Prose-only changes build nothing in CI (D-145, D-147): the workflow's plan job treats *.md,
# LICENSE and NOTICE as prose, except the Markdown files that tests or build scripts read, which
# it lists in a bash array `code_md=( ... )`. This check derives that list from the code and fails
# if plan's list differs, in either direction: a file read but not listed (edits to it would skip
# CI), or a file listed but no longer read. It also fails if the workflow has a workflow-level
# `paths` or `paths-ignore` filter: a run that never starts reports no check, and main's required
# check (ci-ok) would wait forever (D-147).
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

# --- what the workflow's plan job lists as code ------------------------------------------------
file(STRINGS "${WORKFLOW}" _lines)
set(_in_on FALSE)
set(_in_list FALSE)
set(_found FALSE)
set(_filters "")
set(_listed "")
foreach(_line IN LISTS _lines)
    if(_line MATCHES "^on:")
        set(_in_on TRUE)
    elseif(_line MATCHES "^[^ #]")
        set(_in_on FALSE)
    endif()
    if(_in_on AND _line MATCHES "^ +paths(-ignore)?:")
        list(APPEND _filters "${_line}")
    endif()
    if(_line MATCHES "code_md=\\(")
        set(_in_list TRUE)
        set(_found TRUE)
    elseif(_in_list AND _line MATCHES "^ *\\)")
        set(_in_list FALSE)
    elseif(_in_list AND _line MATCHES "^ *'([^']+)'")
        list(APPEND _listed "${CMAKE_MATCH_1}")
    endif()
endforeach()
if(_filters)
    message(FATAL_ERROR "the workflow has a workflow-level path filter (${_filters}); "
                        "skipped runs would never report the required ci-ok check (D-147)")
endif()
if(NOT _found)
    message(FATAL_ERROR "no code_md=( ... ) list found in ${WORKFLOW}")
endif()
list(SORT _listed)

set(_missing ${_read})
if(_listed)
    list(REMOVE_ITEM _missing ${_listed})
endif()
set(_stale ${_listed})
if(_read)
    list(REMOVE_ITEM _stale ${_read})
endif()
if(_missing OR _stale)
    message(FATAL_ERROR "the plan job's code_md list does not match the Markdown files the code reads\n"
                        "  read but not re-included: ${_missing}\n"
                        "  re-included but not read by any test: ${_stale}")
endif()
list(LENGTH _read _n)
message(STATUS "the plan job counts exactly the ${_n} Markdown files the code reads as code: ${_read}")
