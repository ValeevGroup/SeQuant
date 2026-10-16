if (NOT TARGET peglib)
    include(FetchContent)

    FetchContent_Declare(
        peglib
        GIT_REPOSITORY "https://github.com/yhirose/cpp-peglib.git"
		GIT_TAG "${SEQUANT_TRACKED_CPPPEGLIB_TAG}"
        GIT_SHALLOW TRUE
    )

    # Only the library target is wanted; cpp-peglib's own tests and benchmark
    # default to ON and the latter defines a `benchmark` target that collides
    # with Google Benchmark's.
    set(BUILD_TESTS OFF)
    set(PEGLIB_BUILD_BENCHMARK OFF)
    FetchContent_MakeAvailable(peglib)
    unset(BUILD_TESTS)
endif()

# postcond check
if (NOT TARGET peglib)
	message(FATAL_ERROR "FindOrFetchCppPeglib could not make TARGET peglib available")
endif()
