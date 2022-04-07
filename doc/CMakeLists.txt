find_package(Doxygen 1.9.3)
if(DOXYGEN_FOUND)
  configure_file(Doxyfile.in Doxyfile @ONLY IMMEDIATE)
  add_custom_target_subproject(sequant html USES_TERMINAL
    COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/Doxyfile
    SOURCES ${PROJECT_BINARY_DIR}/doc/Doxyfile)
  add_custom_target_subproject(sequant doc DEPENDS html-sequant)
endif()
