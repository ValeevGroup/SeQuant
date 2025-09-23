if (NOT TARGET benchmark::benchmark)
    include(FetchContent)

    set(BENCHMARK_ENABLE_TESTING OFF CACHE INTERNAL "" FORCE)

    FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY "https://github.com/google/benchmark.git"
        GIT_TAG "${SEQUANT_TRACKED_GOOGLEBENCHMARK_TAG}"
        GIT_SHALLOW
        EXCLUDE_FROM_ALL
        SYSTEM
        FIND_PACKAGE_ARGS ${SEQUANT_OLDEST_GOOGLEBENCHMARK_VERSION} NAMES benchmark
    )

    FetchContent_MakeAvailable(googlebenchmark)
endif()

# postcond check
if (NOT TARGET benchmark::benchmark)
    message(FATAL_ERROR "FindOrFetchGoogleBenchmark could not make TARGET benchmark::benchmark available")
endif()
