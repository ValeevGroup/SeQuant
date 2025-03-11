# Main documentation targets
add_custom_target_subproject(doc sequant)

# Find Python for Sphinx dependencies
find_package(Python 3.8 COMPONENTS Interpreter REQUIRED)

# Find Doxygen for API documentation
find_package(Doxygen 1.9.3)
if (DOXYGEN_FOUND)
    set(DOXYGEN_OUTPUT_DIR ${PROJECT_BINARY_DIR}/doc/doxygen)
    file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

    # Configure Doxyfile
    configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in
            ${DOXYGEN_OUTPUT_DIR}/Doxyfile @ONLY IMMEDIATE
    )

    # Set DOT path if available
    if (DOXYGEN_DOT_EXECUTABLE)
        get_filename_component(DOXYGEN_DOT_DIRECTORY ${DOXYGEN_DOT_EXECUTABLE} DIRECTORY)
        file(TO_NATIVE_PATH ${DOXYGEN_DOT_DIRECTORY} DOXYGEN_DOT_PATH)
    endif ()

    # Doxygen build target
    add_custom_target_subproject(sequant doxygen
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM USES_TERMINAL
            SOURCES ${DOXYGEN_OUTPUT_DIR}/Doxyfile
    )

    # Link Doxygen to main documentation target
    add_dependencies(sequant-doc doxygen-sequant)
    set_property(GLOBAL PROPERTY SEQUANT_HAS_DOXYGEN TRUE)
else ()
    set_property(GLOBAL PROPERTY SEQUANT_HAS_DOXYGEN FALSE)
    message(STATUS "Doxygen not found - Skipping API documentation")
endif ()

# Find Sphinx for documentation generation
find_package(Sphinx)
if (SPHINX_FOUND)
    set(SPHINX_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
    set(SPHINX_BUILD_DIR ${PROJECT_BINARY_DIR}/doc/sphinx)
    set(SPHINX_HTML_DIR ${SPHINX_BUILD_DIR}/html)

    # Ensure output directories exist
    file(MAKE_DIRECTORY ${SPHINX_BUILD_DIR} ${SPHINX_HTML_DIR})

    # Check required Python modules from requirements.txt
    file(STRINGS ${SPHINX_SOURCE_DIR}/requirements.txt SPHINX_RAW_MODULES)

    set(SPHINX_PYTHON_MODULES "")
    foreach (module ${SPHINX_RAW_MODULES})
        string(REGEX REPLACE "#.*" "" module_clean ${module})  # Remove comments
        string(STRIP ${module_clean} module_clean)  # Remove whitespace
        if (NOT "${module_clean}" STREQUAL "")
            list(APPEND SPHINX_PYTHON_MODULES ${module_clean})
        endif ()
    endforeach ()

    if (SPHINX_MISSING_MODULES)
        message(FATAL_ERROR
                "The following Python modules required for Sphinx documentation are missing:\n"
                "  ${SPHINX_MISSING_MODULES}\n"
                "Run the following command to install them:\n"
                "  '${Python_EXECUTABLE} -m pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt'\n"
        )
    endif ()

    # copy index.rst to the build directory
    file(COPY ${SPHINX_SOURCE_DIR}/index.rst DESTINATION ${SPHINX_BUILD_DIR})

    # Configure Sphinx conf.py
    configure_file(
            ${SPHINX_SOURCE_DIR}/conf.py.in
            ${SPHINX_BUILD_DIR}/conf.py @ONLY IMMEDIATE
    )

    # Copy documentation sources, excluding _static directory
    file(GLOB DIRS_TO_COPY RELATIVE ${SPHINX_SOURCE_DIR} ${SPHINX_SOURCE_DIR}/*)
    list(REMOVE_ITEM DIRS_TO_COPY "_static")

    foreach (dir ${DIRS_TO_COPY})
        if (IS_DIRECTORY ${SPHINX_SOURCE_DIR}/${dir})
            file(COPY ${SPHINX_SOURCE_DIR}/${dir} DESTINATION ${SPHINX_BUILD_DIR})
        endif ()
    endforeach ()

    # Handle legacy source directory
    if (EXISTS ${SPHINX_SOURCE_DIR}/source)
        file(COPY ${SPHINX_SOURCE_DIR}/source DESTINATION ${SPHINX_BUILD_DIR})
    endif ()

    # Determine if Doxygen is available
    get_property(HAS_DOXYGEN GLOBAL PROPERTY SEQUANT_HAS_DOXYGEN)

    # Sphinx build target
    add_custom_target_subproject(sequant sphinx
            COMMAND ${SPHINX_EXECUTABLE} -E -b html -j auto -a
            -Dbreathe_projects.SeQuant=${PROJECT_BINARY_DIR}/doc/doxygen/xml
            ${SPHINX_BUILD_DIR} ${SPHINX_HTML_DIR}
            WORKING_DIRECTORY ${SPHINX_BUILD_DIR}
            COMMENT "Building HTML documentation with Sphinx"
            VERBATIM USES_TERMINAL
    )

    # Make Sphinx depend on Doxygen if available
    if (HAS_DOXYGEN)
        add_dependencies(sphinx-sequant doxygen-sequant)
    endif ()

    # Add Sphinx to main documentation target
    add_dependencies(sequant-doc sphinx-sequant)
else ()
    message(STATUS "Sphinx not found - Skipping documentation generation")
endif ()

# Alias for HTML generation
add_custom_target_subproject(sequant html DEPENDS sequant-doc)
