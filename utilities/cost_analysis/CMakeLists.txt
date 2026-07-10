include(FindOrFetchCLI11)
include(FindOrFetchJSON)

add_executable(cost_analysis cost_analysis.cpp report.cpp)
set_target_properties(cost_analysis PROPERTIES
    CXX_SCAN_FOR_MODULES OFF
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
target_link_libraries(cost_analysis PRIVATE
    SeQuant::mbpt nlohmann_json::nlohmann_json CLI11::CLI11)
target_set_warning_flags(cost_analysis)

add_test(NAME cost_analysis_ccsd_r2
    COMMAND ${CMAKE_COMMAND}
    -DEXE=$<TARGET_FILE:cost_analysis>
    -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/examples
    -DNAME=ccsd_r2 -DDUMP=R2
    -P ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.cmake)

add_test(NAME cost_analysis_df_r1
    COMMAND ${CMAKE_COMMAND}
    -DEXE=$<TARGET_FILE:cost_analysis>
    -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/examples
    -DNAME=df_r1
    -P ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.cmake)
