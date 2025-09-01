include(FindOrFetchCLI11)
include(FindOrFetchSpdlog)
include(FindOrFetchJSON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

add_executable(external_interface
    external_interface.cpp
    processing.cpp
    utils.cpp
)

target_link_libraries(external_interface
    PRIVATE
        SeQuant::SeQuant
        nlohmann_json::nlohmann_json
        CLI11::CLI11
        Boost::headers
        spdlog::spdlog
)

target_set_warning_flags(external_interface)
