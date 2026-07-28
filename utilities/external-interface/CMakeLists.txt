include(FindOrFetchCLI11)
include(FindOrFetchSpdlog)
include(FindOrFetchJSON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

add_executable(external_interface
    external_interface.cpp
    processing.cpp
    utils.cpp
)
set_target_properties(external_interface PROPERTIES CXX_SCAN_FOR_MODULES OFF)

target_link_libraries(external_interface
    PRIVATE
        SeQuant::mbpt
        SeQuant::export
        nlohmann_json::nlohmann_json
        CLI11::CLI11
        Boost::headers
        spdlog::spdlog
)

target_set_warning_flags(external_interface)

foreach(DRIVER IN ITEMS "ccsd" "nevpt2")
    add_test(
        NAME "sequant/external_interface/${DRIVER}/generate"
        COMMAND external_interface --driver "${CMAKE_CURRENT_LIST_DIR}/examples/${DRIVER}.json"
    )
    set_tests_properties("sequant/external_interface/${DRIVER}/generate" PROPERTIES FIXTURES_SETUP "SEQUANT_TEST_EXTERNAL_INTERFACE_${DRIVER}")

    add_test(
        NAME "sequant/external_interface/${DRIVER}/verify"
        COMMAND "${CMAKE_COMMAND}" -E compare_files "${DRIVER}.itfaa" "${DRIVER}.itfaa.expected"
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/examples/"
    )
    set_tests_properties("sequant/external_interface/${DRIVER}/verify" PROPERTIES FIXTURES_REQUIRED "SEQUANT_TEST_EXTERNAL_INTERFACE_${DRIVER}")
endforeach()
