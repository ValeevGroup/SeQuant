if (NOT TARGET CLI11::CLI11)
    include(FetchContent)

    FetchContent_Declare(
        CLI11
        GIT_REPOSITORY "https://github.com/CLIUtils/CLI11.git"
        GIT_TAG "${SEQUANT_TRACKED_CLI11_TAG}"
        GIT_SHALLOW
        EXCLUDE_FROM_ALL
        SYSTEM
        FIND_PACKAGE_ARGS ${SEQUANT_OLDEST_CLI11_VERSION} NAMES CLI11
    )

    FetchContent_MakeAvailable(CLI11)
endif()

# postcond check
if (NOT TARGET CLI11::CLI11)
    message(FATAL_ERROR "FindOrFetchCLI11 could not make TARGET CLI11::CLI11 available")
endif()
