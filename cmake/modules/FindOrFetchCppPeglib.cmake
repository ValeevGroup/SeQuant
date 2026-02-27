if (NOT TARGET peglib)
    include(FetchContent)

    FetchContent_Declare(
        peglib
        GIT_REPOSITORY "https://github.com/yhirose/cpp-peglib.git"
		GIT_TAG "${SEQUANT_TRACKED_CPPPEGLIB_TAG}"
        GIT_SHALLOW
    )

    FetchContent_MakeAvailable(peglib)
endif()

# postcond check
if (NOT TARGET peglib)
	message(FATAL_ERROR "FindOrFetchCppPeglib could not make TARGET peglib available")
endif()
