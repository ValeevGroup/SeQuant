if (NOT TARGET Catch2::Catch2)
    include(FetchContent)

    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY "https://github.com/catchorg/Catch2.git"
        GIT_TAG "${SEQUANT_TRACKED_CATCH2_TAG}"
        GIT_SHALLOW
        EXCLUDE_FROM_ALL
        SYSTEM
        FIND_PACKAGE_ARGS ${SEQUANT_OLDEST_CATCH2_VERSION} NAMES Catch2
    )

    FetchContent_MakeAvailable(Catch2)
endif()

# postcond check
if (NOT TARGET Catch2::Catch2)
    message(FATAL_ERROR "FindOrFetchCatch2 could not make TARGET Catch2 available")
endif()
