include(FetchContent)

FetchContent_Declare(
	cli11
	GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
	GIT_TAG        v2.4.1
	GIT_SHALLOW    ON
)
FetchContent_Declare(
	json
	URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
)

FetchContent_MakeAvailable(cli11 json)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

add_executable(external_interface
	external_interface.cpp
	processing.cpp
	utils.cpp
)

target_link_libraries(external_interface
	PRIVATE
		SeQuant::SeQuant
		nlohmann_json::nlohmann_json
		CLI11::CLI11
		Boost::headers
)
