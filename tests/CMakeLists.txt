include(FindOrFetchDTL)

option(SEQUANT_TESTS "Enable building of test cases" ${BUILD_TESTING})
add_feature_info(TESTS SEQUANT_TESTS "Enable building of test cases")

if (NOT SEQUANT_TESTS)
    set(BUILD_BY_DEFAULT "EXCLUDE_FROM_ALL")
endif()

option(SEQUANT_SKIP_LONG_TESTS "Whether to skip execution of long-running tests" OFF)

if (SEQUANT_TESTS)
    # Add a "test" that compiles that "all" target, which includes the tests binaries
    add_test("sequant/compile" "${CMAKE_COMMAND}" --build "${PROJECT_BINARY_DIR}")
    set_tests_properties("sequant/compile" PROPERTIES
        FIXTURES_SETUP SEQUANT_BUILD_ALL
        TIMEOUT 0 # disable timeout for the compile step
    )
endif()

function(build_test_as_needed TARGET_NAME TEST_PREFIX DEPENDENT_TESTS)
    if (SEQUANT_TESTS)
        # Tests are built by default -> no need to add extra build steps
        # However, we do depend on the SEQUANT_BUILD_ALL fixture to allow for
        # running ctest without running make first
        foreach(current IN LISTS ${DEPENDENT_TESTS})
            set_tests_properties("${current}" PROPERTIES FIXTURES_REQUIRED "SEQUANT_BUILD_ALL")
        endforeach()

        return()
    endif()

    add_test("${TEST_PREFIX}/build" "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target "${TARGET_NAME}")

    set_tests_properties("${TEST_PREFIX}/build" PROPERTIES FIXTURES_SETUP "SEQUANT_TEST_BUILD_${TARGET_NAME}")

    foreach(current IN LISTS ${DEPENDENT_TESTS})
        set_tests_properties("${current}" PROPERTIES FIXTURES_REQUIRED "SEQUANT_TEST_BUILD_${TARGET_NAME}")
    endforeach()
endfunction()

add_subdirectory(unit)
add_subdirectory(integration)

# target for running ALL tests
add_custom_target_subproject(sequant check
    USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -R "^sequant"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
