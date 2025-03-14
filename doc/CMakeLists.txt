# Main documentation targets
add_custom_target_subproject(sequant doc)

# Find Python for Sphinx dependencies
find_package(Python 3.9 COMPONENTS Interpreter REQUIRED)

# Find Doxygen for API documentation
find_package(Doxygen 1.9.2 REQUIRED)

set(DOXYGEN_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.doxygen")
set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/.doxygen")

file(REAL_PATH "${DOXYGEN_SOURCE_DIR}" DOXYGEN_SOURCE_DIR)
file(REAL_PATH "${DOXYGEN_OUTPUT_DIR}" DOXYGEN_OUTPUT_DIR)

file(MAKE_DIRECTORY "${DOXYGEN_OUTPUT_DIR}")

# Configure Doxyfile
configure_file(
	${DOXYGEN_SOURCE_DIR}/Doxyfile.in
    ${DOXYGEN_OUTPUT_DIR}/Doxyfile @ONLY
)

doxygen_add_docs(
    doxygen-sequant
    COMMENT "Generating API documentation via Doxygen"
    CONFIG_FILE "${DOXYGEN_OUTPUT_DIR}/Doxyfile"
)

# Ensure the config file installed via sphinx-cmake pip package can be found
if ("${Python_EXECUTABLE}" MATCHES "(.+)/bin/python[^/]*")
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_MATCH_1}")
endif()

# Find Sphinx for documentation generation
find_package(Sphinx REQUIRED)

set(SPHINX_DOCS_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(SPHINX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.sphinx")
set(SPHINX_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/.sphinx")
set(SPHINX_HTML_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)

# Directories need to exist in order for REAL_PATH to work properly
file(MAKE_DIRECTORY "${SPHINX_BUILD_DIR}")
file(MAKE_DIRECTORY "${SPHINX_HTML_DIR}")

file(REAL_PATH "${SPHINX_DOCS_DIR}" SPHINX_DOCS_DIR)
file(REAL_PATH "${SPHINX_SOURCE_DIR}" SPHINX_SOURCE_DIR)
file(REAL_PATH "${SPHINX_BUILD_DIR}" SPHINX_BUILD_DIR)
file(REAL_PATH "${SPHINX_HTML_DIR}" SPHINX_HTML_DIR)

# Ensure output directories exist
file(MAKE_DIRECTORY ${SPHINX_BUILD_DIR} ${SPHINX_HTML_DIR})

# Check required Python modules from requirements.txt
file(STRINGS ${SPHINX_SOURCE_DIR}/requirements.txt SPHINX_RAW_MODULES)

set(SPHINX_PYTHON_MODULES "")
set(MODULE_SKIP_CHECK_REGEX "^(sphinx-cmake|sphinx-design)$")
foreach (module ${SPHINX_RAW_MODULES})
    string(REGEX REPLACE "#.*" "" module_clean ${module})  # Remove comments
    string(STRIP ${module_clean} module_clean)  # Remove whitespace
	if (NOT "${module_clean}" STREQUAL "" AND NOT "${module_clean}" MATCHES "${MODULE_SKIP_CHECK_REGEX}")
        list(APPEND SPHINX_PYTHON_MODULES ${module_clean})
    endif ()
endforeach ()

set(SPHINX_MISSING_MODULES "")
foreach(module ${SPHINX_PYTHON_MODULES})
    execute_process(
            COMMAND ${Python_EXECUTABLE} -c "import ${module}"
            RESULT_VARIABLE module_result
            ERROR_QUIET
    )
    if(NOT module_result EQUAL 0)
        list(APPEND SPHINX_MISSING_MODULES ${module})
    endif()
endforeach()

if (SPHINX_MISSING_MODULES)
    message(FATAL_ERROR
            "The following Python modules required for Sphinx documentation are missing:\n"
            "  ${SPHINX_MISSING_MODULES}\n"
            "Run the following command to install them:\n"
            "  '${Python_EXECUTABLE} -m pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/requirements.txt'\n"
    )
endif ()

# Configure Sphinx conf.py
string(TIMESTAMP CURRENT_YEAR "%Y")
configure_file(
        ${SPHINX_SOURCE_DIR}/conf.py.in
        ${SPHINX_BUILD_DIR}/conf.py @ONLY
)

sphinx_add_docs(
    sphinx-sequant
    COMMENT "Generating HTML documentation via Sphinx"
    WORKING_DIRECTORY "${SPHINX_BUILD_DIR}"
    BUILDER html
	SOURCE_DIRECTORY "${SPHINX_DOCS_DIR}"
	CONFIG_DIRECTORY "${SPHINX_BUILD_DIR}"
    OUTPUT_DIRECTORY "${SPHINX_HTML_DIR}"
	DEFINE breathe_projects.SeQuant=${DOXYGEN_OUTPUT_DIR}/xml
    DEPENDS doxygen-sequant
    FRESH_ENV
)

# Add Sphinx to main documentation target
add_dependencies(doc-sequant sphinx-sequant)

# `html` is a standard GNU  target: https://www.gnu.org/prep/standards/html_node/Standard-Targets.html#Standard-Targets
add_custom_target_subproject(sequant html DEPENDS doc-sequant)
