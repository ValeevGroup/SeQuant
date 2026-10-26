if (NOT TARGET nlohmann_json::nlohmann_json)
    include(FetchContent)

    sequant_find_unless_fetching(nlohmann_json ${SEQUANT_OLDEST_JSON_VERSION})

    # a fetched nlohmann_json is installed with the other fetched dependencies; this also
    # covers a fetched and installed TiledArray, whose MADNESS reuses this target and
    # links its installed libraries against it. The export set must exist even when
    # nothing is installed (MADNESS's install(EXPORT) requires it), hence EXCLUDE_FROM_ALL
    # rather than JSON_Install=OFF
    set(_exclude)
    if (NOT SEQUANT_INSTALL_FETCHED_DEPENDENCIES)
        set(_exclude EXCLUDE_FROM_ALL)
    endif()

    FetchContent_Declare(
        nlohmann_json
        URL "https://github.com/nlohmann/json/releases/download/${SEQUANT_TRACKED_JSON_TAG}/json.tar.xz"
        GIT_SHALLOW
        SYSTEM
        ${_exclude}
        FIND_PACKAGE_ARGS ${SEQUANT_OLDEST_JSON_VERSION} NAMES nlohmann_json
    )

    set(JSON_BuildTests OFF CACHE INTERNAL "")
    set(JSON_Install ON CACHE INTERNAL "")

    FetchContent_MakeAvailable(nlohmann_json)

    if (nlohmann_json_FOUND)
        set(SEQUANT_FETCHED_JSON FALSE)
    else()
        set(SEQUANT_FETCHED_JSON TRUE)
        if (SEQUANT_INSTALL_FETCHED_DEPENDENCIES)
            get_directory_property(_json_version DIRECTORY ${nlohmann_json_SOURCE_DIR} DEFINITION nlohmann_json_VERSION)
            sequant_record_installed_dependency(nlohmann_json VERSION ${_json_version})
        endif()
    endif()
endif()

# postcond check
if (NOT TARGET nlohmann_json::nlohmann_json)
    message(FATAL_ERROR "FindOrFetchJSON could not make TARGET nlohmann_json::nlohmann_json available")
endif()
