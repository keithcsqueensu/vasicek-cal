# SPDX-License-Identifier: Apache-2.0
#
# Fails if the TOL_ ids in the test-side register and the documented register differ.
# Usage: cmake -DHEADER=<tests/tolerances.hpp> -DDOC=<docs/methodology/tolerances.md> -P ...

foreach(_var HEADER DOC)
    if(NOT DEFINED ${_var})
        message(FATAL_ERROR "ToleranceRegisterCheck: pass -D${_var}=<path>")
    endif()
endforeach()

function(_vcal_tol_ids file out)
    file(READ "${file}" _text)
    string(REGEX MATCHALL "TOL_[A-Z0-9_]+" _ids "${_text}")
    list(REMOVE_DUPLICATES _ids)
    list(SORT _ids)
    set(${out} "${_ids}" PARENT_SCOPE)
endfunction()

_vcal_tol_ids("${HEADER}" _code)
_vcal_tol_ids("${DOC}" _doc)

set(_only_code ${_code})
list(REMOVE_ITEM _only_code ${_doc})
set(_only_doc ${_doc})
list(REMOVE_ITEM _only_doc ${_code})

if(_only_code OR _only_doc)
    message(FATAL_ERROR "tolerance register out of sync\n"
        "  only in ${HEADER}: ${_only_code}\n"
        "  only in ${DOC}: ${_only_doc}")
endif()
list(LENGTH _code _n)
message(STATUS "tolerance register in sync (${_n} ids)")
