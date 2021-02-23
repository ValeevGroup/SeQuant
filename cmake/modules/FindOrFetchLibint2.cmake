if (NOT TARGET Libint2::cxx)
  find_package(Libint2 2.7.0 CONFIG QUIET)
endif(NOT TARGET Libint2::cxx)

if (TARGET Libint2::libint2_cxx)
  message(STATUS "Found Libint2 CONFIG at ${Libint2_CONFIG}")
else(TARGET Libint2::libint2_cxx)

  if (NOT DEFINED LIBINT2_URL)
    set(LIBINT2_URL https://github.com/evaleev/libint/releases/download/v2.7.0-beta.6/libint-2.7.0-beta.6-test-mpqc4.tgz)
  endif(NOT DEFINED LIBINT2_URL)
  message(STATUS "Will obtain Libint2 from ${LIBINT2_URL}")
  set(URL_SPEC URL "${LIBINT2_URL}")
  if (DEFINED LIBINT2_URL_HASH)
    list(APPEND URL_SPEC URL_HASH "${LIBINT2_URL_HASH}")
  endif(DEFINED LIBINT2_URL_HASH)

  include(FetchContent)
  FetchContent_Declare(
      Libint2
      ${URL_SPEC}
  )
  FetchContent_MakeAvailable(Libint2)
  FetchContent_GetProperties(Libint2
      SOURCE_DIR Libint2_SOURCE_DIR
      BINARY_DIR Libint2_BINARY_DIR
      )
  set_property(DIRECTORY ${Libint2_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL TRUE)

  # ALWAYS build libint using unity build by setting UNITY_BUILD for libint2_obj
  # this saves time and works around issues with Libint library being too big for Ninja on MacOS
  if (NOT CMAKE_UNITY_BUILD)
    set_target_properties(libint2_obj PROPERTIES UNITY_BUILD ON)
    message(STATUS "Will unity-build Libint2")
  endif()

  # use subproject targets as if they were in exported namespace ...
  if (TARGET libint2_cxx AND NOT TARGET Libint2::libint2_cxx)
    add_library(Libint2::libint2_cxx ALIAS libint2_cxx)

    install(TARGETS libint2 EXPORT mpqc
        COMPONENT libint2
        LIBRARY DESTINATION "${MPQC_INSTALL_LIBDIR}")
    install(TARGETS libint2_cxx EXPORT mpqc
        COMPONENT libint2
        LIBRARY DESTINATION "${MPQC_INSTALL_LIBDIR}")
    install(TARGETS libint2_Eigen EXPORT mpqc
        COMPONENT libint2
        LIBRARY DESTINATION "${MPQC_INSTALL_LIBDIR}")
  endif(TARGET libint2_cxx AND NOT TARGET Libint2::libint2_cxx)

  # this is where libint2-config.cmake will end up
  # must be in sync with the "install(FILES ...libint2-config.cmake" statement in https://github.com/evaleev/libint/blob/master/export/cmake/CMakeLists.txt.export
  set(Libint2_CONFIG "${CMAKE_INSTALL_PREFIX}/${LIBINT2_INSTALL_CMAKEDIR}" CACHE INTERNAL "The location of installed libint2-config.cmake file")

endif(TARGET Libint2::libint2_cxx)

# postcond check
if (NOT TARGET Libint2::libint2_cxx)
  message(FATAL_ERROR "FindOrFetchLibint2 could not make Libint2::libint2_cxx target available")
endif(NOT TARGET Libint2::libint2_cxx)
