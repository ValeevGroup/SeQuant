find_package(Doxygen 1.8.12)
if(DOXYGEN_FOUND)
  configure_file(Doxyfile.in Doxyfile @ONLY IMMEDIATE)
  add_custom_target (html USES_TERMINAL
    COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/Doxyfile
    SOURCES ${PROJECT_BINARY_DIR}/doc/Doxyfile)
endif()
