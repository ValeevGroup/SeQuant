# try find_package
if (NOT TARGET tapp::reference)
  include(FindPackageRegimport)
  find_package_regimport(tapp QUIET CONFIG COMPONENTS reference)
  if (TARGET tapp::reference)
    message(STATUS "Found tapp CONFIG at ${tapp_CONFIG}")
    set(SEQUANT_FETCHED_TAPP FALSE)
  elseif (NOT SEQUANT_FETCH_MISSING_DEPENDENCIES)
    message(FATAL_ERROR "TAPP (tapp, component reference) not found: set tapp_DIR to the directory containing tapp-config.cmake, or set SEQUANT_FETCH_MISSING_DEPENDENCIES=ON")
  endif()
endif()

# if not found, build via FetchContent
if (NOT TARGET tapp::reference)
  # a fetched TAPP is installed only on request
  sequant_backend_installable(TAPP _install)
  set(_exclude)
  if (NOT _install)
    set(_exclude EXCLUDE_FROM_ALL)
  endif()

  include(FetchContent)
  FetchContent_Declare(
      tapp
      GIT_REPOSITORY      https://github.com/TAPPorg/reference-implementation.git
      GIT_TAG             ${SEQUANT_TRACKED_TAPP_TAG}
      SYSTEM
      ${_exclude}
  )
  FetchContent_MakeAvailable(tapp)

  set(SEQUANT_FETCHED_TAPP TRUE)
  if (_install)
    get_directory_property(_tapp_version DIRECTORY ${tapp_SOURCE_DIR} DEFINITION TAPP_VERSION)
    # set tapp_CONFIG to the install location so that sequant-config.cmake knows where to find it
    set(tapp_CONFIG ${CMAKE_INSTALL_PREFIX}/${TAPP_INSTALL_CMAKEDIR}/tapp-config.cmake)
    sequant_record_installed_dependency(tapp VERSION ${_tapp_version})
  else()
    unset(tapp_CONFIG)
  endif()

endif()

# postcond check
if (NOT TARGET tapp::reference)
  message(FATAL_ERROR "FindOrFetchTAPP could not make tapp::reference target available")
endif()
