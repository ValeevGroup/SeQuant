# Runs cost_analysis on ${SRC}/${NAME}.json and diffs the report against
# ${SRC}/${NAME}.md.expected (volatile lines stripped). If ${DUMP} is set, also
# checks that the --dump_tree file for that result exists.
execute_process(COMMAND "${EXE}" --driver "${SRC}/${NAME}.json"
                WORKING_DIRECTORY "${SRC}" RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "cost_analysis exited ${rc}")
endif()

function(strip_volatile in out)
  file(STRINGS "${in}" lines)
  set(kept "")
  foreach(line IN LISTS lines)
    if(NOT line MATCHES "^_SeQuant:")
      string(APPEND kept "${line}\n")
    endif()
  endforeach()
  # Normalize trailing blank lines: the committed .md.expected fixture gets
  # its final blank line trimmed by the repo's end-of-file-fixer pre-commit
  # hook, while the tool's live (gitignored) report always ends in one; only
  # the trailing-newline count should differ, so collapse both sides equally.
  string(REGEX REPLACE "\n+$" "\n" kept "${kept}")
  file(WRITE "${out}" "${kept}")
endfunction()

# file(STRINGS) treats a missing input as empty, so a report that was never
# written would slip through strip_volatile as an empty file and misreport as a
# content diff. Fail clearly instead.
if(NOT EXISTS "${SRC}/${NAME}.md")
  message(FATAL_ERROR "report ${NAME}.md was not produced")
endif()
strip_volatile("${SRC}/${NAME}.md" "${SRC}/${NAME}.md.filtered")
strip_volatile("${SRC}/${NAME}.md.expected" "${SRC}/${NAME}.md.expected.filtered")
# compare_files: 0 == identical, 1 == differ, anything else == error; report
# the last distinctly so a broken comparison isn't misread as a content diff.
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
                "${SRC}/${NAME}.md.filtered" "${SRC}/${NAME}.md.expected.filtered"
                RESULT_VARIABLE diff)
if(diff EQUAL 1)
  message(FATAL_ERROR "report differs from expected")
elseif(NOT diff EQUAL 0)
  message(FATAL_ERROR "compare_files failed (rc=${diff})")
endif()

if(DUMP)
  if(NOT EXISTS "${SRC}/${DUMP}.tree.txt")
    message(FATAL_ERROR "expected dump file ${DUMP}.tree.txt not found")
  endif()
endif()
