if (NOT TARGET tiledarray)
  find_package(TiledArray CONFIG COMPONENTS tiledarray)
  if(TiledArray_INSTALL_DIR)
    set(TiledArray_DIR ${TiledArray_INSTALL_DIR}/lib/cmake/tiledarray)
  endif()
  find_package(TiledArray CONFIG QUIET COMPONENTS tiledarray)
endif (NOT TARGET tiledarray)

set(TA_PYTHON OFF)

if (TARGET tiledarray)
  message(STATUS "Found TiledArray CONFIG at ${TiledArray_CONFIG}")

else (TARGET tiledarray)

  set(ENABLE_DQ_PREBUF OFF CACHE BOOL "Whether to enable prebuffering in madness::DQueue")

  include(FetchContent)
  FetchContent_Declare(
      TILEDARRAY
      GIT_REPOSITORY      https://github.com/ValeevGroup/tiledarray.git
      GIT_TAG             ${SEQUANT_TRACKED_TILEDARRAY_TAG}
  )
  FetchContent_MakeAvailable(TILEDARRAY)
  FetchContent_GetProperties(TILEDARRAY
      SOURCE_DIR TILEDARRAY_SOURCE_DIR
      BINARY_DIR TILEDARRAY_BINARY_DIR
      )
  # TA includes dependencies that are built manually, not using FetchContent, hence make sure we build them before building any SeQuant-dependent code
  add_dependencies(tiledarray External-tiledarray)

  include("${TILEDARRAY_BINARY_DIR}/cmake/modules/ReimportTargets.cmake")
  if (NOT TARGET MADworld)
    message(FATAL_ERROR "did not find re-imported target MADworld")
  endif(NOT TARGET MADworld)

  # this is where tiledarray-config.cmake will end up
  # must be in sync with the "install(FILES ...tiledarray-config.cmake" statement in https://github.com/ValeevGroup/tiledarray/blob/${MPQC_TRACKED_TILEDARRAY_TAG}/CMakeLists.txt
  set(TiledArray_CONFIG "${CMAKE_INSTALL_PREFIX}/${TILEDARRAY_INSTALL_CMAKEDIR}" CACHE INTERNAL "The location of installed tiledarray-config.cmake file")
endif(TARGET tiledarray)
