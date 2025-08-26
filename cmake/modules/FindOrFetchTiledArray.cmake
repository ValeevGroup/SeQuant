if (NOT TARGET tiledarray)
    include(FetchContent)

    set(TA_PYTHON OFF)
    set(ENABLE_DQ_PREBUF OFF CACHE BOOL "Whether to enable prebuffering in madness::DQueue" FORCE)

    FetchContent_Declare(
        TiledArray
        GIT_REPOSITORY "https://github.com/ValeevGroup/tiledarray.git"
        GIT_TAG "${SEQUANT_TRACKED_TILEDARRAY_TAG}"
        GIT_SHALLOW
        EXCLUDE_FROM_ALL
        SYSTEM
        FIND_PACKAGE_ARGS NAMES TiledArray COMPONENTS tiledarray
    )

    FetchContent_MakeAvailable(TiledArray)

    if (NOT DEFINED TiledArray_CONFIG)
        # this is where tiledarray-config.cmake will end up
        # must be in sync with the "install(FILES ...tiledarray-config.cmake" statement in https://github.com/ValeevGroup/tiledarray/blob/${MPQC_TRACKED_TILEDARRAY_TAG}/CMakeLists.txt
        set(TiledArray_CONFIG "${CMAKE_INSTALL_PREFIX}/${TILEDARRAY_INSTALL_CMAKEDIR}")
    endif()
endif()

# postcond check
if (NOT TARGET tiledarray)
    message(FATAL_ERROR "FindOrFetchCatch2 could not make TARGET tiledarray available")
endif()
