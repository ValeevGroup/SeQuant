include(FindOrFetchCLI11)
include(FindOrFetchJSON)

add_executable(cost_analysis cost_analysis.cpp report.cpp)
set_target_properties(cost_analysis PROPERTIES
    CXX_SCAN_FOR_MODULES OFF
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
target_link_libraries(cost_analysis PRIVATE
    SeQuant::mbpt nlohmann_json::nlohmann_json CLI11::CLI11)
target_set_warning_flags(cost_analysis)

foreach(NAME IN ITEMS "ccsd_r2" "df_r1")
    add_test(
        NAME "sequant/cost_analysis/${NAME}/generate"
        COMMAND cost_analysis --driver "${CMAKE_CURRENT_LIST_DIR}/examples/${NAME}.json")
    set_tests_properties("sequant/cost_analysis/${NAME}/generate"
        PROPERTIES FIXTURES_SETUP "COST_ANALYSIS_${NAME}")

    # Strip the volatile git-revision line before diffing; -B ignores the extra
    # trailing blank line the live report carries but the fixture doesn't.
    add_test(
        NAME "sequant/cost_analysis/${NAME}/verify"
        COMMAND bash -c "diff -B <(sed '/^_SeQuant:/d' '${NAME}.md') <(sed '/^_SeQuant:/d' '${NAME}.md.expected')"
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/examples")
    set_tests_properties("sequant/cost_analysis/${NAME}/verify"
        PROPERTIES FIXTURES_REQUIRED "COST_ANALYSIS_${NAME}")
endforeach()

# ccsd_r2 exercises --dump_tree; check the tree file was produced.
add_test(
    NAME "sequant/cost_analysis/ccsd_r2/dump_tree"
    COMMAND bash -c "test -f R2.tree.txt"
    WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/examples")
set_tests_properties("sequant/cost_analysis/ccsd_r2/dump_tree"
    PROPERTIES FIXTURES_REQUIRED "COST_ANALYSIS_ccsd_r2")
