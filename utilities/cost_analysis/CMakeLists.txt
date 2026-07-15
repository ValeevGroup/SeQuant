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
    # --omit-revision drops the volatile git-revision line so a plain
    # compare_files suffices; the tool writes the report next to the driver.
    add_test(
        NAME "sequant/cost_analysis/${NAME}/generate"
        COMMAND cost_analysis --driver "${CMAKE_CURRENT_LIST_DIR}/examples/${NAME}.json" --omit-revision)
    set_tests_properties("sequant/cost_analysis/${NAME}/generate"
        PROPERTIES FIXTURES_SETUP "COST_ANALYSIS_${NAME}")

    add_test(
        NAME "sequant/cost_analysis/${NAME}/verify"
        COMMAND "${CMAKE_COMMAND}" -E compare_files "${NAME}.md" "${NAME}.md.expected"
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/examples")
    set_tests_properties("sequant/cost_analysis/${NAME}/verify"
        PROPERTIES FIXTURES_REQUIRED "COST_ANALYSIS_${NAME}")
endforeach()

# ccsd_r2 exercises --dump_tree; diff the tree dump against its fixture.
add_test(
    NAME "sequant/cost_analysis/ccsd_r2/dump_tree"
    COMMAND "${CMAKE_COMMAND}" -E compare_files "R2.tree.txt" "R2.tree.txt.expected"
    WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/examples")
set_tests_properties("sequant/cost_analysis/ccsd_r2/dump_tree"
    PROPERTIES FIXTURES_REQUIRED "COST_ANALYSIS_ccsd_r2")
