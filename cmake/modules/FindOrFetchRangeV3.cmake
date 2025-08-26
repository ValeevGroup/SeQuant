if (NOT TARGET range-v3::range-v3)
    include(FetchContent)

    FetchContent_Declare(
        RangeV3
        GIT_REPOSITORY "https://github.com/ericniebler/range-v3.git"
        GIT_TAG "${SEQUANT_TRACKED_RANGEV3_TAG}"
        GIT_SHALLOW
        EXCLUDE_FROM_ALL
        SYSTEM
        FIND_PACKAGE_ARGS NAMES range-v3
    )

    FetchContent_MakeAvailable(RangeV3)

    if (NOT DEFINED range-v3_CONFIG)
        set(range-v3_CONFIG "${CMAKE_INSTALL_PREFIX}/lib/cmake/range-v3/range-v3-config.cmake")
    endif()
endif()

# postcond check
if (NOT TARGET range-v3::range-v3)
  message(FATAL_ERROR "FindOrFetchRangeV3 could not make range-v3::range-v3 target available")
endif()
