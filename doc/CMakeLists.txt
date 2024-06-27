find_package(Doxygen 1.9.3)
if(DOXYGEN_FOUND)
  configure_file(Doxyfile.in Doxyfile @ONLY IMMEDIATE)
  add_custom_target_subproject(sequant html USES_TERMINAL
    COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/Doxyfile
    SOURCES ${PROJECT_BINARY_DIR}/doc/Doxyfile)
  add_custom_target_subproject(sequant doc DEPENDS html-sequant)

  # if DOT found obtain DOXYGEN_DOT_PATH
  if (DOXYGEN_DOT_EXECUTABLE)
    get_filename_component(DOXYGEN_DOT_DIRECTORY ${DOXYGEN_DOT_EXECUTABLE} DIRECTORY)
    file(TO_NATIVE_PATH ${DOXYGEN_DOT_DIRECTORY} DOXYGEN_DOT_PATH)
  endif()
endif()

# Sphinx
find_package(Sphinx REQUIRED)
if (SPHINX_FOUND)
  message(STATUS "Sphinx-build found: ${SPHINX_EXECUTABLE}")
  set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
  set(SPHINX_BUILD ${PROJECT_BINARY_DIR}/doc/sphinx)
  file(MAKE_DIRECTORY ${SPHINX_BUILD})

  add_custom_target_subproject(sequant sphinx USES_TERMINAL COMMAND ${SPHINX_EXECUTABLE} -b html -Dbreathe_projects.SeQuant=${PROJECT_BINARY_DIR}/doc/xml ${SPHINX_SOURCE} ${SPHINX_BUILD} WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc DEPENDS html-sequant)
else ()
  message(FATAL_ERROR "Sphinx not found")
endif ()
