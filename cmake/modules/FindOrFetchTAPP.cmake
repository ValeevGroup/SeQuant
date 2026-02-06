# try find_package
if (NOT TARGET tapp-reference)
  include(FindPackageRegimport)
  find_package_regimport(tapp-reference QUIET CONFIG)
  if (TARGET tapp-reference)
    message(STATUS "Found tapp-reference CONFIG at ${tapp-reference_CONFIG}")
  endif (TARGET tapp-reference)
endif (NOT TARGET tapp-reference)

# if not found, build via FetchContent
if (NOT TARGET tapp-reference)

  include(FetchContent)
  FetchContent_Declare(
      tapp-reference
      GIT_REPOSITORY      https://github.com/TAPPorg/reference-implementation.git
      GIT_TAG             ${SEQUANT_TRACKED_TAPP_TAG}
      EXCLUDE_FROM_ALL
      SYSTEM
  )
  FetchContent_MakeAvailable(tapp-reference)
  FetchContent_GetProperties(tapp-reference
      SOURCE_DIR TAPP_SOURCE_DIR
      BINARY_DIR TAPP_BINARY_DIR
  )

  set(SEQUANT_TAPP_BUILT_FROM_SOURCE TRUE)

  # The TAPP reference implementation does not provide install/export config,
  # so we patch its targets and install them as part of SeQuant.

  # The TAPP upstream CMakeLists.txt uses bare source paths for include
  # directories and INTERFACE_SOURCES which CMake rejects during install/export.
  # Patch them to use proper BUILD_INTERFACE/INSTALL_INTERFACE generator exprs.

  # Fix tapp-reference include directories
  set_property(TARGET tapp-reference PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  target_include_directories(tapp-reference PUBLIC
      $<BUILD_INTERFACE:${TAPP_SOURCE_DIR}/reference_implementation/include>
      $<INSTALL_INTERFACE:${SEQUANT_INSTALL_INCLUDEDIR}>
  )

  # Fix tapp-api: clear INTERFACE_SOURCES (bare source paths not exportable;
  # these are just for IDE display, headers are found via include dirs)
  set_property(TARGET tapp-api PROPERTY INTERFACE_SOURCES)

  # Install both tapp-api (INTERFACE) and tapp-reference (SHARED) in the
  # sequant export set so that downstream find_package(SeQuant) works.
  install(TARGETS tapp-api
          EXPORT sequant
          COMPONENT sequant)
  install(TARGETS tapp-reference
          EXPORT sequant
          COMPONENT sequant
          LIBRARY DESTINATION ${SEQUANT_INSTALL_LIBDIR})

  # Install headers from both the API and the reference implementation
  if (EXISTS "${TAPP_SOURCE_DIR}/api/include")
    install(DIRECTORY "${TAPP_SOURCE_DIR}/api/include/"
            COMPONENT sequant
            DESTINATION "${SEQUANT_INSTALL_INCLUDEDIR}")
  endif()
  if (EXISTS "${TAPP_SOURCE_DIR}/reference_implementation/include")
    install(DIRECTORY "${TAPP_SOURCE_DIR}/reference_implementation/include/"
            COMPONENT sequant
            DESTINATION "${SEQUANT_INSTALL_INCLUDEDIR}")
  endif()

endif(NOT TARGET tapp-reference)

# postcond check
if (NOT TARGET tapp-reference)
  message(FATAL_ERROR "FindOrFetchTAPP could not make tapp-reference target available")
endif(NOT TARGET tapp-reference)
