if (NOT TARGET nlohmann_json::nlohmann_json)
    include(FetchContent)

    FetchContent_Declare(
        nlohmann_json
        URL "https://github.com/nlohmann/json/releases/download/${SEQUANT_TRACKED_JSON_TAG}/json.tar.xz"
        GIT_SHALLOW
        SYSTEM
        FIND_PACKAGE_ARGS ${SEQUANT_OLDEST_JSON_VERSION} NAMES nlohmann_json
    )

    # a fetched nlohmann_json must be installable and exported: dependencies
    # fetched later (MADNESS via TiledArray) reuse this target when they find
    # it already defined and link their installed libraries against it
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    set(JSON_Install ON CACHE INTERNAL "")

    FetchContent_MakeAvailable(nlohmann_json)
endif()

# postcond check
if (NOT TARGET nlohmann_json::nlohmann_json)
    message(FATAL_ERROR "FindOrFetchJSON could not make TARGET nlohmann_json::nlohmann_json available")
endif()
