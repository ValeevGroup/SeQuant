include(FindOrFetchCLI11)
include(FindOrFetchJSON)

add_executable(cost_analysis cost_analysis.cpp report.cpp)
set_target_properties(cost_analysis PROPERTIES
    CXX_SCAN_FOR_MODULES OFF
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
target_link_libraries(cost_analysis PRIVATE
    SeQuant::mbpt nlohmann_json::nlohmann_json CLI11::CLI11)
target_set_compiler_flags(cost_analysis)

# The tool writes its report (and any dump-tree file) next to the driver, so
# the tests run against a build-tree copy of examples/ rather than the source
# directory -- keeps the source tree untouched and lets read-only / out-of-source
# builds pass. The copy is refreshed each ctest run by this setup fixture.
set(_examples_src "${CMAKE_CURRENT_LIST_DIR}/examples")
set(_examples_bin "${CMAKE_CURRENT_BINARY_DIR}/examples")
add_test(
    NAME "sequant/cost_analysis/prepare"
    COMMAND "${CMAKE_COMMAND}" -E copy_directory "${_examples_src}" "${_examples_bin}")
set_tests_properties("sequant/cost_analysis/prepare"
    PROPERTIES FIXTURES_SETUP "COST_ANALYSIS_EXAMPLES")

foreach(NAME IN ITEMS "ccsd_r2" "df_r1")
    # --omit-revision drops the volatile git-revision line so a plain
    # compare_files suffices; drive the build-tree copy prepared above.
    add_test(
        NAME "sequant/cost_analysis/${NAME}/generate"
        COMMAND cost_analysis --driver "${_examples_bin}/${NAME}.json" --omit-revision)
    set_tests_properties("sequant/cost_analysis/${NAME}/generate"
        PROPERTIES FIXTURES_REQUIRED "COST_ANALYSIS_EXAMPLES"
                   FIXTURES_SETUP "COST_ANALYSIS_${NAME}")

    add_test(
        NAME "sequant/cost_analysis/${NAME}/verify"
        COMMAND "${CMAKE_COMMAND}" -E compare_files "${NAME}.md" "${NAME}.md.expected"
        WORKING_DIRECTORY "${_examples_bin}")
    set_tests_properties("sequant/cost_analysis/${NAME}/verify"
        PROPERTIES FIXTURES_REQUIRED "COST_ANALYSIS_${NAME}")
endforeach()

# ccsd_r2 exercises --dump_tree; diff the tree dump against its fixture.
add_test(
    NAME "sequant/cost_analysis/ccsd_r2/dump_tree"
    COMMAND "${CMAKE_COMMAND}" -E compare_files "R2.tree.txt" "R2.tree.txt.expected"
    WORKING_DIRECTORY "${_examples_bin}")
set_tests_properties("sequant/cost_analysis/ccsd_r2/dump_tree"
    PROPERTIES FIXTURES_REQUIRED "COST_ANALYSIS_ccsd_r2")
