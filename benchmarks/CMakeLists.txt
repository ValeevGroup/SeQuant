include(FindOrFetchGoogleBenchmark)

add_executable(sequant_benchmarks
	"canonicalize.cpp"
	"main.cpp"
	"simplify.cpp"
	"spintrace.cpp"
	"tensor_network.cpp"
	"wick.cpp"
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
