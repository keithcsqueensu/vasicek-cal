# SPDX-License-Identifier: Apache-2.0
#
# P-11 (D-179): a study run stopped midway and resumed from its checkpoint writes output identical to
# an uninterrupted run's. study_parametric_bootstrap on 3 scenarios x 2 replicates: (1) uninterrupted;
# (2) with a checkpoint directory, stopped abruptly after the first scenario (--stop-after-scenarios 1,
# which exits as a kill would, after that scenario's unit is on disk); (3) the same command again,
# which must resume (skipping the first scenario) and write a byte-identical CSV; (4) a different
# configuration pointed at the same directory, which must be refused.
#
#   cmake -DTOOL=<study_parametric_bootstrap> -DWORK=<scratch dir> -P resume_test.cmake

file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
set(ARGS --replicates 2 --scenarios 72,51,4)

execute_process(COMMAND "${TOOL}" ${ARGS} --out "${WORK}/uninterrupted.csv" RESULT_VARIABLE rc OUTPUT_QUIET ERROR_QUIET)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "uninterrupted run failed (${rc})")
endif()

execute_process(COMMAND "${TOOL}" ${ARGS} --checkpoint "${WORK}/ck" --stop-after-scenarios 1 --out "${WORK}/stopped.csv"
                RESULT_VARIABLE rc OUTPUT_QUIET ERROR_VARIABLE err)
if(NOT rc EQUAL 3)
    message(FATAL_ERROR "the stopped run should exit with 3 after one scenario, got ${rc}: ${err}")
endif()
if(EXISTS "${WORK}/stopped.csv")
    message(FATAL_ERROR "the stopped run wrote output; it should have stopped before")
endif()

execute_process(COMMAND "${TOOL}" ${ARGS} --checkpoint "${WORK}/ck" --out "${WORK}/resumed.csv"
                RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "resumed run failed (${rc}): ${err}")
endif()
if(NOT out MATCHES "resumed: 1 of 3 scenarios")
    message(FATAL_ERROR "the resumed run did not take the first scenario from the checkpoint: ${out}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${WORK}/uninterrupted.csv" "${WORK}/resumed.csv"
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "the resumed run's output differs from the uninterrupted run's")
endif()

execute_process(COMMAND "${TOOL}" --replicates 3 --scenarios 72,51,4 --checkpoint "${WORK}/ck" --out "${WORK}/other.csv"
                RESULT_VARIABLE rc OUTPUT_QUIET ERROR_VARIABLE err)
if(NOT rc EQUAL 2 OR NOT err MATCHES "different configuration")
    message(FATAL_ERROR "a different configuration was not refused (${rc}): ${err}")
endif()
message(STATUS "resume: identical output after a stop and resume; a different configuration refused")
