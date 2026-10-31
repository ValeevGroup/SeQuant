if (NOT TARGET tiledarray)
    include(FetchContent)

    set(TA_PYTHON OFF)
    set(ENABLE_DQ_PREBUF OFF CACHE BOOL "Whether to enable prebuffering in madness::DQueue" FORCE)

    # Seed TiledArray's assertion policy from SeQuant's, so that the two agree
    # by default; the THROW/ABORT/IGNORE values map one-to-one. Like any cached
    # option this applies on the first configure of a build directory only: an
    # explicit -DTA_ASSERT_POLICY=... (from the user or a parent project) is
    # already in the cache and wins, and a later change of
    # SEQUANT_ASSERT_BEHAVIOR does not re-seed it (set TA_ASSERT_POLICY
    # explicitly, or use a fresh build directory). If TiledArray ends up being
    # found installed instead, the seeded entry is simply unused.
    # BTAS is TiledArray's business: it builds BTAS itself and seeds
    # BTAS_ASSERT_POLICY from TA_ASSERT_POLICY the same way.
    if (NOT DEFINED CACHE{TA_ASSERT_POLICY})
        set(TA_ASSERT_POLICY TA_ASSERT_${SEQUANT_ASSERT_BEHAVIOR} CACHE STRING "Controls the behavior of TA_ASSERT (seeded from SEQUANT_ASSERT_BEHAVIOR)")
    endif()

    sequant_find_unless_fetching(TiledArray COMPONENTS tiledarray)

    # a fetched TiledArray (with MADNESS and BTAS) is installed only on request;
    # the nlohmann_json it reuses is made installable by FindOrFetchJSON
    sequant_backend_installable(TiledArray _install)
    set(_exclude)
    if (NOT _install)
        set(_exclude EXCLUDE_FROM_ALL)
    endif()

    FetchContent_Declare(
        TiledArray
        GIT_REPOSITORY "https://github.com/ValeevGroup/tiledarray.git"
        GIT_TAG "${SEQUANT_TRACKED_TILEDARRAY_TAG}"
        GIT_SHALLOW
        SYSTEM
        ${_exclude}
        FIND_PACKAGE_ARGS NAMES TiledArray COMPONENTS tiledarray
    )

    FetchContent_MakeAvailable(TiledArray)

    if (TiledArray_FOUND)
        set(SEQUANT_FETCHED_TILEDARRAY FALSE)
    else()
        set(SEQUANT_FETCHED_TILEDARRAY TRUE)
        if (_install)
            get_directory_property(_ta_version DIRECTORY ${tiledarray_SOURCE_DIR} DEFINITION TILEDARRAY_VERSION)
            # this is where tiledarray-config.cmake will end up
            # must be in sync with the "install(FILES ...tiledarray-config.cmake" statement in https://github.com/ValeevGroup/tiledarray/blob/${MPQC_TRACKED_TILEDARRAY_TAG}/CMakeLists.txt
            set(TiledArray_CONFIG "${CMAKE_INSTALL_PREFIX}/${TILEDARRAY_INSTALL_CMAKEDIR}/tiledarray-config.cmake")
            sequant_record_installed_dependency(TiledArray VERSION ${_ta_version})
        else()
            unset(TiledArray_CONFIG)
        endif()
    endif()
endif()

# postcond check
if (NOT TARGET tiledarray)
	message(FATAL_ERROR "FindOrFetchTiledArray could not make TARGET tiledarray available")
endif()
