if (NOT TARGET range-v3::range-v3)
    include(FetchContent)

    # range-v3's version file is ExactVersion, so look for it without a version
    sequant_find_unless_fetching(range-v3)

    # range-v3 has no option to disable its install rules
    set(_exclude)
    if (NOT SEQUANT_INSTALL_FETCHED_DEPENDENCIES)
        set(_exclude EXCLUDE_FROM_ALL)
    endif()

    FetchContent_Declare(
        range-v3
        GIT_REPOSITORY "https://github.com/ericniebler/range-v3.git"
        GIT_TAG "${SEQUANT_TRACKED_RANGEV3_TAG}"
        GIT_SHALLOW
        SYSTEM
        ${_exclude}
        FIND_PACKAGE_ARGS
    )

    FetchContent_MakeAvailable(range-v3)

    if (range-v3_FOUND)
        set(_rangev3_version "${range-v3_VERSION}")
        set(SEQUANT_FETCHED_RANGEV3 FALSE)
    else()
        get_directory_property(_rangev3_version DIRECTORY ${range-v3_SOURCE_DIR} DEFINITION RANGE_V3_VERSION)
        set(SEQUANT_FETCHED_RANGEV3 TRUE)
    endif()
    if (_rangev3_version VERSION_LESS SEQUANT_OLDEST_RANGEV3_VERSION)
        message(FATAL_ERROR "Found range-v3 ${_rangev3_version}, but SeQuant requires range-v3 ${SEQUANT_OLDEST_RANGEV3_VERSION} or later")
    endif()

    if (NOT SEQUANT_FETCHED_RANGEV3)
        # range-v3_CONFIG is the found config
    elseif (SEQUANT_INSTALL_FETCHED_DEPENDENCIES)
        set(range-v3_CONFIG "${CMAKE_INSTALL_PREFIX}/lib/cmake/range-v3/range-v3-config.cmake")
        sequant_record_installed_dependency(range-v3 VERSION ${_rangev3_version})
    else()
        unset(range-v3_CONFIG)
    endif()
endif()

# postcond check
if (NOT TARGET range-v3::range-v3)
  message(FATAL_ERROR "FindOrFetchRangeV3 could not make range-v3::range-v3 target available")
endif()
