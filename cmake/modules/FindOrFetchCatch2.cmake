include(FetchContent)

if (NOT TARGET Catch2::Catch2)
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

# A fetched Catch2 puts its extras/ dir (Catch.cmake, for catch_discover_tests)
# on CMAKE_MODULE_PATH only in the directory scope that fetched it, and its
# targets are global, so a sibling directory that includes this module skips the
# fetch and would fail include(Catch). When find_package satisfied the fetch
# instead, catch2_SOURCE_DIR is unset and each including directory's own
# find_package call already extends CMAKE_MODULE_PATH.
FetchContent_GetProperties(Catch2)
if (catch2_SOURCE_DIR AND EXISTS "${catch2_SOURCE_DIR}/extras/Catch.cmake"
    AND NOT "${catch2_SOURCE_DIR}/extras" IN_LIST CMAKE_MODULE_PATH)
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
endif()
