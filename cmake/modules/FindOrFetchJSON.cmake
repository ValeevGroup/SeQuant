if (NOT TARGET nlohmann_json::nlohmann_json)
    include(FetchContent)

    sequant_find_unless_fetching(nlohmann_json ${SEQUANT_OLDEST_JSON_VERSION})

    FetchContent_Declare(
        nlohmann_json
        URL "https://github.com/nlohmann/json/releases/download/${SEQUANT_TRACKED_JSON_TAG}/json.tar.xz"
        GIT_SHALLOW
        SYSTEM
        FIND_PACKAGE_ARGS ${SEQUANT_OLDEST_JSON_VERSION} NAMES nlohmann_json
    )

    # a fetched nlohmann_json is installed with the other fetched dependencies; this also
    # covers a fetched and installed TiledArray, whose MADNESS reuses this target and
    # links its installed libraries against it
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    set(JSON_Install ${SEQUANT_INSTALL_FETCHED_DEPENDENCIES} CACHE INTERNAL "")

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
