include(FindOrFetchCLI11)
include(FindOrFetchJSON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

add_executable(cost_analysis cost_analysis.cpp report.cpp)
set_target_properties(cost_analysis PROPERTIES CXX_SCAN_FOR_MODULES OFF)
target_link_libraries(cost_analysis PRIVATE
    SeQuant::mbpt nlohmann_json::nlohmann_json CLI11::CLI11)
target_set_warning_flags(cost_analysis)
