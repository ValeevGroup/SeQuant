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
		# The CMake option is not honored if spdlog is found via CMake config file
		# However, our code doesn't compile when using libfmt as we're specializing templates for std::format
		# see https://github.com/gabime/spdlog/issues/3467
		#FIND_PACKAGE_ARGS NAMES spdlog
    )

    FetchContent_MakeAvailable(spdlog)
endif()

# postcond check
if (NOT TARGET spdlog::spdlog)
    message(FATAL_ERROR "FindOrFetchSpdlog could not make TARGET spdlog::spdlog available")
endif()
