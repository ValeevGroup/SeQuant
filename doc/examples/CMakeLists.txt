# Make sure all example snippets are compiled to avoid documented examples becoming out of date

option(SEQUANT_COMPILE_DOC_EXAMPLES "Whether to compile all code examples from the documentation" ${PROJECT_IS_TOP_LEVEL})

if (NOT SEQUANT_COMPILE_DOC_EXAMPLES)
    return()
endif()

file(GLOB_RECURSE EXAMPLE_CODES LIST_DIRECTORIES false CONFIGURE_DEPENDS "*.cpp")

set(counter 1)
foreach(current IN LISTS EXAMPLE_CODES)
    file(READ "${current}" contents)

    set(current_name "sequant_doc_example_${counter}")

    if ("${contents}" MATCHES " main\\(")
        add_executable(${current_name} ${current})
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

    math(EXPR counter "${counter} + 1")
endforeach()
