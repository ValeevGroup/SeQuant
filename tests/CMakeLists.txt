option(SEQUANT_TESTS "Enable building of test cases" ${BUILD_TESTING})
add_feature_info(TESTS SEQUANT_TESTS "Enable building of test cases")

if (NOT SEQUANT_TESTS)
    set(BUILD_BY_DEFAULT "EXCLUDE_FROM_ALL")
endif()

function(build_test_as_needed TARGET_NAME TEST_PREFIX DEPENDENT_TESTS)
    if (SEQUANT_TESTS)
        # Tests are built by default -> no need to add extra build steps
        return()
    endif()

    add_test("${TEST_PREFIX}/build" "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target "${TARGET_NAME}")

    foreach(current IN LISTS ${DEPENDENT_TESTS})
        set_tests_properties("${current}" PROPERTIES DEPENDS "${TEST_PREFIX}/build")
    endforeach()
endfunction()

add_subdirectory(unit)
add_subdirectory(integration)

# target for running ALL tests
add_custom_target_subproject(sequant check
    USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -R "^sequant")

