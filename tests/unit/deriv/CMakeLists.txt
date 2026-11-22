# Unit tests of SeQuant-deriv (symbolic differentiation), kept in their own executable
add_executable(unit_tests-sequant-deriv ${BUILD_BY_DEFAULT}
    "../test_main.cpp"
    "contexts.cpp"
    "evaluator.cpp"
    "hf_model.cpp"
    "identities.cpp"
    "judge.cpp"
    "test_deriv_hf.cpp"
    "test_deriv_hf_numeric.cpp"
    "test_deriv_identities.cpp"
    "test_deriv_judge.cpp"
    "test_deriv_rules.cpp"
    "test_deriv_series.cpp"
)
set_target_properties(unit_tests-sequant-deriv PROPERTIES CXX_SCAN_FOR_MODULES OFF)
# SeQuant::mbpt: test_main.cpp sets up the mbpt Context; dtl: catch2_sequant.hpp;
# Eigen: hf_model.cpp
target_link_libraries(unit_tests-sequant-deriv PRIVATE
    SeQuant::deriv
    SeQuant::mbpt
    Catch2::Catch2
    dtl::dtl
    Eigen3::Eigen)
target_include_directories(unit_tests-sequant-deriv PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_compile_definitions(unit_tests-sequant-deriv PRIVATE
    SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/..")
set_source_files_properties("../test_main.cpp"
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
# range-v3 names ::meta unqualified, ambiguous with sequant::meta under the
# file-scope `using namespace sequant` that a Unity build leaks into later files
set_property(TARGET unit_tests-sequant-deriv PROPERTY
    UNITY_BUILD_CODE_BEFORE_INCLUDE "#include <range/v3/all.hpp>")
target_set_compiler_flags(unit_tests-sequant-deriv)
target_set_optimization_flags(unit_tests-sequant-deriv)

if (SEQUANT_TESTS)
    catch_discover_tests(
        unit_tests-sequant-deriv
        TEST_PREFIX "sequant/unit/"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        PROPERTIES FIXTURES_REQUIRED "SEQUANT_BUILD_ALL"
    )
else()
    add_test(
        NAME "sequant/unit/deriv"
        COMMAND unit_tests-sequant-deriv
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )

    set(TEST_NAMES "sequant/unit/deriv")
    build_test_as_needed(unit_tests-sequant-deriv "sequant/unit/deriv" TEST_NAMES)
endif()
