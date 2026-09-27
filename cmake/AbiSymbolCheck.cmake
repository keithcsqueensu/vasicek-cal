# SPDX-License-Identifier: Apache-2.0
#
# The shared library must export exactly the functions include/vcal/vcal.h declares: nothing
# missing, and nothing else, such as a leaked C++ symbol (D-142).
# Usage: cmake -DLIBRARY=<file> -DHEADER=<vcal.h> -DTOOL=<nm or dumpbin> -DKIND=nm|dumpbin
#              -P AbiSymbolCheck.cmake
# Declarations are the header lines of the form "VCAL_API <type> VCAL_CALL vcal_<name>(".

cmake_minimum_required(VERSION 3.25)  # script mode: set policies (IN_LIST needs CMP0057)

foreach(_v LIBRARY HEADER TOOL KIND)
    if(NOT DEFINED ${_v})
        message(FATAL_ERROR "AbiSymbolCheck: pass -D${_v}=...")
    endif()
endforeach()

file(STRINGS "${HEADER}" _decls REGEX "^VCAL_API .* VCAL_CALL vcal_[a-z0-9_]+[(]")
set(_declared "")
foreach(_line IN LISTS _decls)
    string(REGEX REPLACE "^.* VCAL_CALL (vcal_[a-z0-9_]+)[(].*$" "\\1" _name "${_line}")
    list(APPEND _declared "${_name}")
endforeach()

if(KIND STREQUAL "nm")
    execute_process(COMMAND "${TOOL}" -D --defined-only "${LIBRARY}"
        OUTPUT_VARIABLE _out RESULT_VARIABLE _rc ERROR_VARIABLE _err)
    set(_pattern "^[0-9a-fA-F]+ [A-Za-z] ([^ ]+)$")
elseif(KIND STREQUAL "dumpbin")
    execute_process(COMMAND "${TOOL}" /NOLOGO /EXPORTS "${LIBRARY}"
        OUTPUT_VARIABLE _out RESULT_VARIABLE _rc ERROR_VARIABLE _err)
    set(_pattern "^ +[0-9]+ +[0-9A-Fa-f]+ +[0-9A-Fa-f]+ +([^ ]+)")
else()
    message(FATAL_ERROR "AbiSymbolCheck: KIND must be nm or dumpbin")
endif()
if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "AbiSymbolCheck: ${TOOL} failed (${_rc}): ${_err}")
endif()

string(REPLACE "\r" "" _out "${_out}")
string(REPLACE "\n" ";" _lines "${_out}")
set(_exported "")
# Symbols the GNU linker defines itself (section bounds, init/fini) are not the library's.
set(_linker_defined _edata _end __bss_start _init _fini)
foreach(_line IN LISTS _lines)
    if(_line MATCHES "${_pattern}")
        if(NOT CMAKE_MATCH_1 IN_LIST _linker_defined)
            list(APPEND _exported "${CMAKE_MATCH_1}")
        endif()
    endif()
endforeach()

list(SORT _declared)
list(SORT _exported)
set(_missing ${_declared})
if(_exported)
    list(REMOVE_ITEM _missing ${_exported})
endif()
set(_extra ${_exported})
if(_declared)
    list(REMOVE_ITEM _extra ${_declared})
endif()
list(LENGTH _declared _n)
if(_missing OR _extra OR _n EQUAL 0)
    message(FATAL_ERROR "exports differ from vcal.h (${_n} declared)\n  missing: ${_missing}\n  extra: ${_extra}")
endif()
message(STATUS "exports match vcal.h: ${_n} functions")
