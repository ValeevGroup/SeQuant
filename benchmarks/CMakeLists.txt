include(FindOrFetchGoogleBenchmark)

add_executable(sequant_benchmarks
        "canonicalize.cpp"
        "coupled_cluster.cpp"
        "main.cpp"
        "simplify.cpp"
        "spintrace.cpp"
        "tensor_block_compare.cpp"
        "tensor_network.cpp"
        "wick.cpp"
)

target_set_warning_flags(sequant_benchmarks)
target_set_optimization_flags(sequant_benchmarks)

set_target_properties(sequant_benchmarks PROPERTIES CXX_SCAN_FOR_MODULES OFF)

option(SEQUANT_QUICK_BENCHMARKS "Whether to exclude computational expensive benchmarks" OFF)

# Running benchmarks in Debug mode doesn't make a lot of sense but in case someone does,
# we at least want to exclude benchmarks that are likely to take forever in this setup
if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR SEQUANT_QUICK_BENCHMARKS)
	target_compile_definitions(sequant_benchmarks
		PRIVATE
			SEQUANT_BENCH_MAX_CC_RANK=3
			SEQUANT_BENCH_MAX_WICK_INPUTS=2
	)
endif()

# The benchmark implementations contain a bunch of same-name symbol definition that
# ought to be local to the individual translation unit. Hence, we must not combine
# all of them into a single translation unit.
set_target_properties(sequant_benchmarks PROPERTIES UNITY_BUILD OFF)

target_link_libraries(sequant_benchmarks
        PRIVATE
        benchmark::benchmark
        SeQuant::SeQuant
)
