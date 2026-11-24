include(FindOrFetchDTL)

option(SEQUANT_TESTS "Enable building of test cases" ${BUILD_TESTING})
add_feature_info(TESTS SEQUANT_TESTS "Enable building of test cases")

if (NOT SEQUANT_TESTS)
    set(BUILD_BY_DEFAULT "EXCLUDE_FROM_ALL")
endif()

option(SEQUANT_SKIP_LONG_TESTS "Whether to skip execution of long-running tests" OFF)

# The eval tests are few, very large translation units, each a deep stack of
# Catch2, eval and (for the TiledArray backend) TiledArray templates: one
# compile of the larger tests/unit/test_eval_ta*.cpp runs for minutes and takes
# most of the memory of a small machine, and two of them side by side do not fit a
# 16 GB CI runner. Setting this to N lets at most N of these compiles run at a
# time, whatever the build's -j; everything else still builds in parallel.
# Honored by the Ninja generators only (JOB_POOLS).
set(SEQUANT_TESTS_HEAVY_COMPILE_JOBS "" CACHE STRING
    "Maximum number of concurrent compiles of the heaviest test sources (empty: no limit; Ninja only)")
if (SEQUANT_TESTS_HEAVY_COMPILE_JOBS)
    set_property(GLOBAL APPEND PROPERTY JOB_POOLS
        sequant_heavy_compile=${SEQUANT_TESTS_HEAVY_COMPILE_JOBS})
endif()

# puts the compiles of TARGET_NAME into the pool of heavy compiles, if any
function(sequant_heavy_compile TARGET_NAME)
    if (SEQUANT_TESTS_HEAVY_COMPILE_JOBS)
        set_target_properties(${TARGET_NAME} PROPERTIES
            JOB_POOL_COMPILE sequant_heavy_compile)
    endif()
endfunction()

# Generator expression so that multi-config generators decide per configuration
set(SEQUANT_SKIP_LONG_TESTS_GENEX "$<OR:$<BOOL:${SEQUANT_SKIP_LONG_TESTS}>,$<CONFIG:Debug>>")

if (SEQUANT_TESTS)
    # Add a "test" that builds that "all" target, which includes the tests binaries
    add_test("sequant/build" "${CMAKE_COMMAND}" --build "${PROJECT_BINARY_DIR}")
    set_tests_properties("sequant/build" PROPERTIES
        FIXTURES_SETUP SEQUANT_BUILD_ALL
        TIMEOUT 0 # disable timeout for the build step
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
