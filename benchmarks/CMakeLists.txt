include(FetchContent)

FetchContent_Declare(
	googlebenchmark
	GIT_REPOSITORY https://github.com/google/benchmark
	GIT_TAG        v1.8.3
	GIT_SHALLOW    ON
)

set(BENCHMARK_ENABLE_TESTING OFF CACHE INTERNAL "")
set(BENCHMARK_ENABLE_LTO OFF CACHE INTERNAL "${CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE}")
set(BENCHMARK_INSTALL OFF CACHE INTERNAL "")
set(BENCHMARK_INSTALL_DOCS OFF CACHE INTERNAL "")

FetchContent_MakeAvailable(googlebenchmark)

add_executable(wick_benchmark wick_benchmark.cpp)
target_link_libraries(wick_benchmark PRIVATE benchmark::benchmark SeQuant::SeQuant)
