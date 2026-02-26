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

set_target_properties(sequant_benchmarks PROPERTIES CXX_SCAN_FOR_MODULES OFF)

# The benchmark implementations contain a bunch of same-name symbol definition that
# ought to be local to the individual translation unit. Hence, we must not combine
# all of them into a single translation unit.
set_target_properties(sequant_benchmarks PROPERTIES UNITY_BUILD OFF)

target_link_libraries(sequant_benchmarks
        PRIVATE
        benchmark::benchmark
        SeQuant::SeQuant
)
