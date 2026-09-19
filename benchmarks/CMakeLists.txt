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

target_set_compiler_flags(sequant_benchmarks)
target_set_optimization_flags(sequant_benchmarks)

set_target_properties(sequant_benchmarks PROPERTIES CXX_SCAN_FOR_MODULES OFF)

# Running benchmarks in Debug mode doesn't make a lot of sense but in case someone does,
# we at least want to exclude benchmarks that are likely to take forever in this setup.
# Use a generator expression rather than comparing CMAKE_BUILD_TYPE: for a multi-config
# generator (Visual Studio, Ninja Multi-Config) that variable is always empty at configure
# time, since the actual configuration is chosen later, per build invocation.
target_compile_definitions(sequant_benchmarks
    PRIVATE
        $<$<OR:$<CONFIG:Debug>,$<BOOL:${SEQUANT_QUICK_BENCHMARKS}>>:SEQUANT_BENCH_MAX_CC_RANK=3>
        $<$<OR:$<CONFIG:Debug>,$<BOOL:${SEQUANT_QUICK_BENCHMARKS}>>:SEQUANT_BENCH_MAX_WICK_INPUTS=2>
)

# The benchmark implementations contain a bunch of same-name symbol definition that
# ought to be local to the individual translation unit. Hence, we must not combine
# all of them into a single translation unit.
set_target_properties(sequant_benchmarks PROPERTIES UNITY_BUILD OFF)

target_link_libraries(sequant_benchmarks
        PRIVATE
        benchmark::benchmark
        SeQuant::SeQuant
)
