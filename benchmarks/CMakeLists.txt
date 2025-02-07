include(FindOrFetchGoogleBenchmark)

add_executable(sequant_benchmarks
	"main.cpp"
	"wick.cpp"
)

target_link_libraries(sequant_benchmarks
	PRIVATE
		benchmark::benchmark
		SeQuant::SeQuant
)
