include(FindOrFetchCatch2)
include(FindOrFetchDTL)
include(Catch)

add_executable(unit_tests-external-interface
    test_main.cpp
    test_execution_context.cpp
    test_executor.cpp
    test_steps.cpp
    test_expr_steps.cpp
    test_tree_steps.cpp
    test_json_include.cpp
    test_template_expansion.cpp
)
set_target_properties(unit_tests-external-interface PROPERTIES CXX_SCAN_FOR_MODULES OFF)

# external_interface-obj carries the processing pipeline (and PUBLIC-propagates
# SeQuant::export, nlohmann_json, ... and its own include dir); the tests/unit
# include dir provides the shared Catch2 helpers (catch2_sequant.hpp) and the
# export index-space context (test_export.hpp), reused rather than duplicated.
target_link_libraries(unit_tests-external-interface
    PRIVATE
        external_interface-obj
        Catch2::Catch2
        dtl::dtl
)
target_include_directories(unit_tests-external-interface PRIVATE
    "${PROJECT_SOURCE_DIR}/tests/unit")
target_compile_definitions(unit_tests-external-interface PRIVATE
    SEQUANT_EXTINT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

target_set_compiler_flags(unit_tests-external-interface)
target_set_optimization_flags(unit_tests-external-interface)

# As in tests/unit: a file-scope `using namespace sequant` would, under a unity
# build, leak into sources concatenated after it and make a range-v3 header
# parsed there ambiguous (range-v3's ::meta vs. sequant::meta). Parse range-v3
# ahead of every file instead (include guards make this free after the first).
set_property(TARGET unit_tests-external-interface PROPERTY
    UNITY_BUILD_CODE_BEFORE_INCLUDE "#include <range/v3/all.hpp>")

catch_discover_tests(
    unit_tests-external-interface
    TEST_PREFIX "sequant/external_interface/unit/"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    PROPERTIES FIXTURES_REQUIRED "SEQUANT_BUILD_ALL"
)
