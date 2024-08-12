# Sphinx
find_package(Sphinx)
if (SPHINX_FOUND)

    # Check if required python modules are available
    file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt SPHINX_PYTHON_MODULES)
    foreach (module ${SPHINX_PYTHON_MODULES})
        execute_process(COMMAND python3 -c "import ${module}"
                RESULT_VARIABLE PYMODULE_CHECK_EXITCODE ERROR_QUIET)
        if (NOT PYMODULE_CHECK_EXITCODE EQUAL 0)
            message(WARNING "Docs: Required python module ${module} is not found")
        endif()
    endforeach()

    # Configure conf.py and copy source files to binary tree, otherwise a lot files will be generated in source tree
    configure_file(conf.py.in conf.py @ONLY IMMEDIATE)
    file(COPY index.rst DESTINATION ${PROJECT_BINARY_DIR}/doc/sphinx)
    # Copy source files to Sphinx root, because everything needs to be relative to conf.py
    file(COPY source DESTINATION ${PROJECT_BINARY_DIR}/doc/sphinx)
    file(COPY ${PROJECT_SOURCE_DIR}/INSTALL.md DESTINATION ${PROJECT_BINARY_DIR}/doc/sphinx/source)

    set(SPHINX_SOURCE ${PROJECT_BINARY_DIR}/doc/sphinx)
    set(SPHINX_BUILD ${PROJECT_BINARY_DIR}/doc/sphinx/html)
    file(MAKE_DIRECTORY ${SPHINX_BUILD})

    add_custom_target_subproject(sequant sphinx USES_TERMINAL COMMAND ${SPHINX_EXECUTABLE} -b html
            -Dbreathe_projects.SeQuant=${PROJECT_BINARY_DIR}/doc/doxygen/xml ${SPHINX_SOURCE} ${SPHINX_BUILD}
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc DEPENDS html-sequant)
else ()
    message(WARNING "Sphinx not found, not building Sphinx documentation.")
endif ()
