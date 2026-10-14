if (NOT TARGET libperm)
    include(FetchContent)

    # libPerm's own tests are never run by check-sequant (ctest -R "^sequant"),
    # and gtest_discover_tests runs the test binary at build time with a 5 s
    # timeout that a loaded CI runner can exceed
    set(LIBPERM_TESTS OFF CACHE BOOL "Whether to build libPerm's own tests" FORCE)

    FetchContent_Declare(
        libPerm
        GIT_REPOSITORY "https://github.com/Krzmbrzl/libPerm.git"
        GIT_TAG "${SEQUANT_TRACKED_LIBPERM_TAG}"
        GIT_SHALLOW
        SYSTEM
    )

    FetchContent_MakeAvailable(libPerm)
endif()

# postcond check
if (NOT TARGET libperm)
    message(FATAL_ERROR "FindOrFetchLibPerm could not make TARGET libperm available")
endif()
