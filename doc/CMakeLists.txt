# Main documentation targets
add_custom_target_subproject(sequant doc)

# Find Python for Sphinx dependencies
find_package(Python 3.9 COMPONENTS Interpreter REQUIRED)

# Find Doxygen for API documentation
find_package(Doxygen 1.9.2 REQUIRED)

set(DOXYGEN_OUTPUT_DIR ${PROJECT_BINARY_DIR}/doc/doxygen)
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

# Check for dot support
if(DOXYGEN_HAVE_DOT)
	set(DOXYGEN_DOT_FOUND YES)
endif()

# Configure Doxyfile
configure_file(
		${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in
		${DOXYGEN_OUTPUT_DIR}/Doxyfile @ONLY
)

# Doxygen build target
add_custom_target_subproject(sequant doxygen
		COMMAND Doxygen::doxygen ${DOXYGEN_OUTPUT_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		COMMENT "Generating API documentation with Doxygen"
		VERBATIM USES_TERMINAL
		SOURCES ${DOXYGEN_OUTPUT_DIR}/Doxyfile
)

# Link Doxygen to main documentation target
add_dependencies(doc-sequant doxygen-sequant)

# Find Sphinx for documentation generation
find_package(Sphinx REQUIRED)

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

# copy index.rst to the build directory
file(COPY ${SPHINX_SOURCE_DIR}/index.rst DESTINATION ${SPHINX_BUILD_DIR})

# Configure Sphinx conf.py
string(TIMESTAMP CURRENT_YEAR "%Y")
configure_file(
		${SPHINX_SOURCE_DIR}/conf.py.in
		${SPHINX_BUILD_DIR}/conf.py @ONLY
)

# Copy documentation sources, excluding _static directory
set(SOURCE_DIRS user internal api)

foreach (dir ${SOURCE_DIRS})
	if (IS_DIRECTORY ${SPHINX_SOURCE_DIR}/${dir})
		file(COPY ${SPHINX_SOURCE_DIR}/${dir} DESTINATION ${SPHINX_BUILD_DIR})
	endif ()
endforeach ()

# Sphinx build target
add_custom_target_subproject(sequant sphinx
		COMMAND ${SPHINX_EXECUTABLE} -E -b html -j auto -a
		-Dbreathe_projects.SeQuant=${PROJECT_BINARY_DIR}/doc/doxygen/xml
		${SPHINX_BUILD_DIR} ${SPHINX_HTML_DIR}
		WORKING_DIRECTORY ${SPHINX_BUILD_DIR}
		COMMENT "Building HTML documentation with Sphinx"
		VERBATIM USES_TERMINAL
)
add_dependencies(sphinx-sequant doxygen-sequant)

# Add Sphinx to main documentation target
add_dependencies(doc-sequant sphinx-sequant)

# `html` is a standard GNU  target: https://www.gnu.org/prep/standards/html_node/Standard-Targets.html#Standard-Targets
add_custom_target_subproject(sequant html DEPENDS doc-sequant)
