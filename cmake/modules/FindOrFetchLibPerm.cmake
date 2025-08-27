if (NOT TARGET libperm)
    include(FetchContent)

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
