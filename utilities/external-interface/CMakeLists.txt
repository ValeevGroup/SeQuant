include(FindOrFetchCLI11)
include(FindOrFetchSpdlog)
include(FindOrFetchJSON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

add_executable(external_interface
    external_interface.cpp
    processing.cpp
    utils.cpp

    canonicalize_step.cpp
    cse_step.cpp
    density_fitting_step.cpp
    execution_context.cpp
    executor.cpp
    export_step.cpp
    index_batching_step.cpp
    optimization_step.cpp
    output_step.cpp
    processing_step.cpp
    processing_step_factory.cpp
    projection_step.cpp
    read_input_step.cpp
    simplify_step.cpp
    spintracing_step.cpp
    to_export_tree_step.cpp
    validate_step.cpp
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


# CCSD
add_test(
	NAME "sequant/external_interface/ccsd/generate"
	COMMAND external_interface --driver "${CMAKE_CURRENT_LIST_DIR}/examples/ccsd.json"
)
set_tests_properties("sequant/external_interface/ccsd/generate" PROPERTIES FIXTURES_SETUP "SEQUANT_TEST_EXTERNAL_INTERFACE_CCSD")

add_test(
	NAME "sequant/external_interface/ccsd/verify"
	COMMAND "${CMAKE_COMMAND}" -E compare_files "ccsd.itfaa" "ccsd.itfaa.expected"
	WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/examples/"
)
set_tests_properties("sequant/external_interface/ccsd/verify" PROPERTIES FIXTURES_REQUIRED "SEQUANT_TEST_EXTERNAL_INTERFACE_CCSD")


# NEVPT2 v1
# Without CSE
add_test(
	NAME "sequant/external_interface/nevpt2_v1/generate"
	COMMAND external_interface --driver "${CMAKE_CURRENT_LIST_DIR}/examples/nevpt2_v1.json"
)
set_tests_properties("sequant/external_interface/nevpt2_v1/generate" PROPERTIES FIXTURES_SETUP "SEQUANT_TEST_EXTERNAL_INTERFACE_NEVPT2v1")

add_test(
	NAME "sequant/external_interface/nevpt2_v1/verify"
	COMMAND "${CMAKE_COMMAND}" -E compare_files "nevpt2_v1.itfaa" "nevpt2.itfaa.expected"
	WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/examples/"
)
set_tests_properties("sequant/external_interface/nevpt2_v1/verify" PROPERTIES FIXTURES_REQUIRED "SEQUANT_TEST_EXTERNAL_INTERFACE_NEVPT2v1")

# With CSE
add_test(
	NAME "sequant/external_interface/nevpt2_v1_cse/generate"
	COMMAND external_interface --driver "${CMAKE_CURRENT_LIST_DIR}/examples/nevpt2_v1_cse.json"
)
set_tests_properties("sequant/external_interface/nevpt2_v1_cse/generate" PROPERTIES FIXTURES_SETUP "SEQUANT_TEST_EXTERNAL_INTERFACE_NEVPT2v1CSE")

add_test(
	NAME "sequant/external_interface/nevpt2_v1_cse/verify"
	COMMAND "${CMAKE_COMMAND}" -E compare_files "nevpt2_v1_cse_merged.itfaa" "nevpt2_cse_merged.itfaa.expected"
	WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/examples/"
)
set_tests_properties("sequant/external_interface/nevpt2_v1_cse/verify" PROPERTIES FIXTURES_REQUIRED "SEQUANT_TEST_EXTERNAL_INTERFACE_NEVPT2v1CSE")


# NEVPT2 v2
add_test(
	NAME "sequant/external_interface/nevpt2_v2/generate"
	COMMAND external_interface --driver "${CMAKE_CURRENT_LIST_DIR}/examples/nevpt2_v2.json"
)
set_tests_properties("sequant/external_interface/nevpt2_v2/generate" PROPERTIES FIXTURES_SETUP "SEQUANT_TEST_EXTERNAL_INTERFACE_NEVPT2v2")

add_test(
	NAME "sequant/external_interface/nevpt2_v2/verify"
	COMMAND "${CMAKE_COMMAND}" -E compare_files "nevpt2_v2.itfaa" "nevpt2.itfaa.expected"
	WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/examples/"
)
set_tests_properties("sequant/external_interface/nevpt2_v2/verify" PROPERTIES FIXTURES_REQUIRED "SEQUANT_TEST_EXTERNAL_INTERFACE_NEVPT2v2")

add_test(
	NAME "sequant/external_interface/nevpt2_v2_cse/verify"
	COMMAND "${CMAKE_COMMAND}" -E compare_files "nevpt2_v2_cse.itfaa" "nevpt2_cse.itfaa.expected"
	WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/examples/"
)
set_tests_properties("sequant/external_interface/nevpt2_v2_cse/verify" PROPERTIES FIXTURES_REQUIRED "SEQUANT_TEST_EXTERNAL_INTERFACE_NEVPT2v2")

add_test(
	NAME "sequant/external_interface/nevpt2_v2_cse_merged/verify"
	COMMAND "${CMAKE_COMMAND}" -E compare_files "nevpt2_v2_cse_merged.itfaa" "nevpt2_cse_merged.itfaa.expected"
	WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/examples/"
)
set_tests_properties("sequant/external_interface/nevpt2_v2_cse_merged/verify" PROPERTIES FIXTURES_REQUIRED "SEQUANT_TEST_EXTERNAL_INTERFACE_NEVPT2v2")
