# Make sure all example snippets are compiled to avoid documented examples becoming out of date

option(SEQUANT_COMPILE_DOC_EXAMPLES "Whether to compile all code examples from the documentation" ${PROJECT_IS_TOP_LEVEL})

if (NOT SEQUANT_COMPILE_DOC_EXAMPLES)
    return()
endif()

file(GLOB_RECURSE EXAMPLE_CODES LIST_DIRECTORIES false CONFIGURE_DEPENDS "*.cpp")

# Create a list to track runnable examples
set(RUNNABLE_DOC_EXAMPLES "")

set(counter 1)
foreach(current IN LISTS EXAMPLE_CODES)
    file(READ "${current}" contents)

    set(current_name "sequant_doc_example_${counter}")

    if ("${contents}" MATCHES " main\\(")
        add_executable(${current_name} ${current})

        # check if running is disabled, we do not want to run very expensive examples
        if (NOT "${contents}" MATCHES "// SEQUANT_DOC_EXAMPLE: DO NOT RUN")
            list(APPEND RUNNABLE_DOC_EXAMPLES ${current_name})
        endif ()
    else()
        add_library(${current_name} ${current})
    endif()

    target_link_libraries(${current_name} PRIVATE SeQuant)

    # Make sure the user sees as little from these as possible
    set_target_properties(${current_name}
        PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )

    target_set_warning_flags("${current_name}")

    math(EXPR counter "${counter} + 1")
endforeach()

# Testing runnable examples
if (SEQUANT_TESTS AND RUNNABLE_DOC_EXAMPLES)
    foreach (example IN LISTS RUNNABLE_DOC_EXAMPLES)
        add_test(NAME "sequant/doc-examples/${example}/build"
                 COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target "${example}")
        set_tests_properties("sequant/doc-examples/${example}/build" PROPERTIES FIXTURES_SETUP SEQUANT_FIXTURE_${example})

        # these tests should fail if any runtime issue is encountered
        add_test(
                NAME "sequant/doc-examples/${example}/run"
                COMMAND ${example}
        )
        set_tests_properties("sequant/doc-examples/${example}/run" PROPERTIES
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                FIXTURES_REQUIRED SEQUANT_FIXTURE_${example})

    endforeach ()
endif ()
