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

    # a fetched nlohmann_json must be installable and exported when a
    # dependency fetched later (MADNESS via TiledArray) reuses this target and
    # links its installed libraries against it; until that decision is made
    # (it turns JSON_Install ON again and records nlohmann_json when TiledArray
    # is installed) it follows SEQUANT_INSTALL_FETCHED_DEPENDENCIES
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    set(JSON_Install ${SEQUANT_INSTALL_FETCHED_DEPENDENCIES} CACHE INTERNAL "")

    FetchContent_MakeAvailable(nlohmann_json)

    if (nlohmann_json_FOUND)
        set(SEQUANT_FETCHED_JSON FALSE)
    else()
        set(SEQUANT_FETCHED_JSON TRUE)
    endif()
endif()

# postcond check
if (NOT TARGET nlohmann_json::nlohmann_json)
    message(FATAL_ERROR "FindOrFetchJSON could not make TARGET nlohmann_json::nlohmann_json available")
endif()
