add_subdirectory(doxygen)

# Sphinx

#find_package(Sphinx)
#if (SPHINX_FOUND)
#  message(STATUS "Sphinx-build found: ${SPHINX_EXECUTABLE}")
#  set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
#  set(SPHINX_BUILD ${PROJECT_BINARY_DIR}/doc/sphinx)
#  file(MAKE_DIRECTORY ${SPHINX_BUILD})
#
#  add_custom_target_subproject(sequant sphinx USES_TERMINAL COMMAND ${SPHINX_EXECUTABLE} -b html
#          -Dbreathe_projects.SeQuant=${PROJECT_BINARY_DIR}/doc/xml ${SPHINX_SOURCE} ${SPHINX_BUILD}
#          WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc DEPENDS html-sequant)
#else ()
#  message(WARNING "Sphinx not found, not building Sphinx documentation.")
#endif ()
