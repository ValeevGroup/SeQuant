set(TEST_CASES
    "antisymmetrizer.cpp"
)
if (NOT SEQUANT_SKIP_LONG_TESTS)
    list(APPEND TEST_CASES
        # Single-Reference Coupled-Cluster equation generation (spin-orbital and spin-free)
        "srcc.cpp ->  |3 t std so|3 t csv so|3 lambda std so|3 lambda csv so|2 t std sf|2 t csv sf|2 lambda std sf|2 lambda csv sf"
        # Single-Reference open-shell equation generation (spin-traced)
        "osstcc.cpp"
        # Equation-of-motion Coupled-Cluster
        "eomcc.cpp"
            "os_eom.cpp"
    )

    if (TARGET Eigen3::Eigen)
        # these examples require Eigen for full functionality
        list(APPEND TEST_CASES
            # Single-Reference closed-shell Coupled-Cluster equation generation
            "stcc.cpp -> |2 v1 naive|2 v1|2 v2|3 v1|3 v2|4 v1|4 v2"
        )
    endif()
else()
    list(APPEND TEST_CASES
        # Single-Reference Coupled-Cluster equation generation (reduced test set)
        "srcc.cpp ->  |2 t csv sf"
    )
    if (TARGET Eigen3::Eigen)
        # these examples require Eigen for full functionality
        list(APPEND TEST_CASES
                # Single-Reference closed-shell Coupled-Cluster equation generation
                "stcc.cpp -> |2 v1 naive|3 v1|3 v2"
        )
    endif()
endif()

foreach(current IN LISTS TEST_CASES)
    if (NOT "${current}" MATCHES "([a-zA-z0-9_-]+\\.cpp) *(-> *(.+))?")
        message(FATAL_ERROR "Test case string has invalid format: '${current}'")
    endif()
    set(test_source "${CMAKE_MATCH_1}")
    set(test_variants "${CMAKE_MATCH_3}")

    string(REGEX REPLACE "\\.cpp$" "" test_name "${test_source}")

    add_executable(${test_name} ${BUILD_BY_DEFAULT} ${test_source})
    target_link_libraries(${test_name} PRIVATE SeQuant)

    if (SEQUANT_SKIP_LONG_TESTS)
        target_compile_definitions(${test_name} PRIVATE SEQUANT_SKIP_LONG_TESTS=1)
    endif()

    target_set_warning_flags("${test_name}")

    set(variant_names "")

    if (NOT test_variants)
        add_test(
            NAME "sequant/integration/${test_name}"
            COMMAND ${test_name}
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        )
        list(APPEND variant_names "sequant/integration/${test_name}")
    else()
        # By inserting semicolons, we turn test_variants into a list
        string(REPLACE "|" ";" test_variants "${test_variants}")

        set(counter 1)
        foreach(arguments IN LISTS test_variants)
            add_test(
                NAME "sequant/integration/${test_name}/${counter}"
                COMMAND ${test_name} ${arguments}
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            )
            list(APPEND variant_names "sequant/integration/${test_name}/${counter}")

            math(EXPR counter "${counter} + 1")
        endforeach()
    endif()

    build_test_as_needed(${test_name} "sequant/integration/${test_name}" variant_names)
endforeach()

add_subdirectory(eval)
