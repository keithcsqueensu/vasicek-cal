# SPDX-License-Identifier: Apache-2.0
#
# D-198 (as P-11's resume_test.cmake for S-10): study_lrpd_serial_m6 stopped midway and resumed from its
# checkpoint writes output identical to an uninterrupted run's, and the output does not depend on the thread
# count or on how refit tasks interleave. Three subset scenarios x 2 replicates at phi = 0.6, on the smoke seed
# with a 3-panel comparator: (1) uninterrupted; (2) the same on one thread; (3) with a checkpoint directory,
# stopped abruptly after two replicates are written; (4) the same command again, which must resume them and
# write a byte-identical CSV; (5) a different configuration pointed at the same directory, refused.
#
#   cmake -DTOOL=<study_lrpd_serial_m6> -DWORK=<scratch dir> -P m6_resume_test.cmake

file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
set(ARGS --replicates 2 --scenarios 72,51,4 --phis 0.6 --seed 0x5AFE --comparator --boot-panels 3)

execute_process(COMMAND "${TOOL}" ${ARGS} --out "${WORK}/uninterrupted.csv" RESULT_VARIABLE rc OUTPUT_QUIET ERROR_QUIET)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "uninterrupted run failed (${rc})")
endif()

execute_process(COMMAND "${TOOL}" ${ARGS} --threads 1 --out "${WORK}/one_thread.csv" RESULT_VARIABLE rc OUTPUT_QUIET
                ERROR_QUIET)
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${WORK}/uninterrupted.csv" "${WORK}/one_thread.csv"
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "the one-thread run's output differs from the multi-thread run's")
endif()

execute_process(COMMAND "${TOOL}" ${ARGS} --checkpoint "${WORK}/ck" --stop-after-units 2 --out "${WORK}/stopped.csv"
                RESULT_VARIABLE rc OUTPUT_QUIET ERROR_VARIABLE err)
if(NOT rc EQUAL 3)
    message(FATAL_ERROR "the stopped run should exit with 3 after two replicates, got ${rc}: ${err}")
endif()
if(EXISTS "${WORK}/stopped.csv")
    message(FATAL_ERROR "the stopped run wrote output; it should have stopped before")
endif()

execute_process(COMMAND "${TOOL}" ${ARGS} --checkpoint "${WORK}/ck" --out "${WORK}/resumed.csv"
                RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "resumed run failed (${rc}): ${err}")
endif()
if(NOT out MATCHES "resumed: 2 of 6 replicates")
    message(FATAL_ERROR "the resumed run did not take two replicates from the checkpoint: ${out}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${WORK}/uninterrupted.csv" "${WORK}/resumed.csv"
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "the resumed run's output differs from the uninterrupted run's")
endif()

execute_process(COMMAND "${TOOL}" --replicates 3 --scenarios 72,51,4 --phis 0.6 --seed 0x5AFE --comparator
                        --boot-panels 3 --checkpoint "${WORK}/ck" --out "${WORK}/other.csv"
                RESULT_VARIABLE rc OUTPUT_QUIET ERROR_VARIABLE err)
if(NOT rc EQUAL 2 OR NOT err MATCHES "different configuration")
    message(FATAL_ERROR "a different configuration was not refused (${rc}): ${err}")
endif()
message(STATUS "M6 resume: identical output after a stop and resume, and on one thread; a different configuration refused")
