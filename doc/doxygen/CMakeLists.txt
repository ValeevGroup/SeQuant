find_package(Doxygen 1.9.3)
if(DOXYGEN_FOUND)
    configure_file(Doxyfile.in Doxyfile @ONLY IMMEDIATE)
    add_custom_target_subproject(sequant html USES_TERMINAL
            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile
            SOURCES ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile)
    add_custom_target_subproject(sequant doc DEPENDS html-sequant)

    # if DOT found obtain DOXYGEN_DOT_PATH
    if (DOXYGEN_DOT_EXECUTABLE)
        get_filename_component(DOXYGEN_DOT_DIRECTORY ${DOXYGEN_DOT_EXECUTABLE} DIRECTORY)
        file(TO_NATIVE_PATH ${DOXYGEN_DOT_DIRECTORY} DOXYGEN_DOT_PATH)
    endif()
endif()
