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

    # a fetched nlohmann_json must also be installed when TiledArray will be
    # fetched and installed: MADNESS reuses this target and links its installed
    # libraries against it; whether TiledArray will be found instead is not
    # known yet, so at worst nlohmann_json is installed needlessly
    sequant_backend_installable(TiledArray _ta_install)
    if (SEQUANT_INSTALL_FETCHED_DEPENDENCIES OR (SEQUANT_TILEDARRAY AND _ta_install AND NOT TARGET tiledarray))
        set(_json_install ON)
    else()
        set(_json_install OFF)
    endif()
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    set(JSON_Install ${_json_install} CACHE INTERNAL "")

    FetchContent_MakeAvailable(nlohmann_json)

    if (nlohmann_json_FOUND)
        set(SEQUANT_FETCHED_JSON FALSE)
    else()
        set(SEQUANT_FETCHED_JSON TRUE)
        if (_json_install)
            get_directory_property(_json_version DIRECTORY ${nlohmann_json_SOURCE_DIR} DEFINITION nlohmann_json_VERSION)
            sequant_record_installed_dependency(nlohmann_json VERSION ${_json_version})
        endif()
    endif()
endif()

# postcond check
if (NOT TARGET nlohmann_json::nlohmann_json)
    message(FATAL_ERROR "FindOrFetchJSON could not make TARGET nlohmann_json::nlohmann_json available")
endif()
