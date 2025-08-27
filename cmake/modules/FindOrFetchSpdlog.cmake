if (NOT TARGET spdlog::spdlog)
    include(FetchContent)

	set(SPDLOG_USE_STD_FORMAT ON)

    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY "https://github.com/gabime/spdlog.git"
        GIT_TAG "${SEQUANT_TRACKED_SPDLOG_TAG}"
        GIT_SHALLOW
        EXCLUDE_FROM_ALL
        SYSTEM
        FIND_PACKAGE_ARGS NAMES spdlog
    )

    FetchContent_MakeAvailable(spdlog)
endif()

# postcond check
if (NOT TARGET spdlog::spdlog)
    message(FATAL_ERROR "FindOrFetchSpdlog could not make TARGET spdlog::spdlog available")
endif()
