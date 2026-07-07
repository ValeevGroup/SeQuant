# Runs cost_analysis on ${SRC}/${NAME}.json and diffs the report against
# ${SRC}/${NAME}.md.expected (volatile lines stripped). If ${DUMP} is set, also
# checks that the --dump_tree file for that result exists.
execute_process(COMMAND ${EXE} --driver ${SRC}/${NAME}.json
                WORKING_DIRECTORY ${SRC} RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "cost_analysis exited ${rc}")
endif()

function(strip_volatile in out)
  file(STRINGS ${in} lines)
  set(kept "")
  foreach(line IN LISTS lines)
    if(NOT line MATCHES "^_Generated:" AND NOT line MATCHES "^_SeQuant:")
      string(APPEND kept "${line}\n")
    endif()
  endforeach()
  # Normalize trailing blank lines: the committed .md.expected fixture gets
  # its final blank line trimmed by the repo's end-of-file-fixer pre-commit
  # hook, while the tool's live (gitignored) report always ends in one; only
  # the trailing-newline count should differ, so collapse both sides equally.
  string(REGEX REPLACE "\n+$" "\n" kept "${kept}")
  file(WRITE ${out} "${kept}")
endfunction()

strip_volatile(${SRC}/${NAME}.md ${SRC}/${NAME}.md.filtered)
strip_volatile(${SRC}/${NAME}.md.expected ${SRC}/${NAME}.md.expected.filtered)
execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files
                ${SRC}/${NAME}.md.filtered ${SRC}/${NAME}.md.expected.filtered
                RESULT_VARIABLE diff)
if(NOT diff EQUAL 0)
  message(FATAL_ERROR "report differs from expected")
endif()

if(DUMP)
  if(NOT EXISTS ${SRC}/${DUMP}.tree.txt)
    message(FATAL_ERROR "expected dump file ${DUMP}.tree.txt not found")
  endif()
endif()
