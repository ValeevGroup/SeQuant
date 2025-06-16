find_package(benchmark)

if (NOT TARGET benchmark::benchmark)
    include(${vg_cmake_kit_SOURCE_DIR}/modules/VRGFindOrFetchPackage.cmake)

	set(BENCHMARK_ENABLE_TESTING OFF CACHE INTERNAL "" FORCE)

    VRGFindOrFetchPackage(benchmark "https://github.com/google/benchmark.git" "${SEQUANT_TRACKED_GOOGLEBENCHMARK_TAG}"
            ADD_SUBDIR
			ADD_SUBDIR_EXCLUDE_FROM_ALL
            CONFIG_SUBDIR
    )
endif()

# postcond check
if (NOT TARGET benchmark::benchmark)
	message(FATAL_ERROR "FindOrFetchGoogleBenchmark could not make TARGET benchmark::benchmark available")
endif()
