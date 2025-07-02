include(FetchContent)

FetchContent_Declare(
	DTL
	GIT_REPOSITORY https://github.com/cubicdaiya/dtl
	GIT_TAG v1.21
	GIT_SHALLOW TRUE
	EXCLUDE_FROM_ALL
)

FetchContent_MakeAvailable(DTL)

add_library(sequant_dtl INTERFACE)
target_include_directories(sequant_dtl INTERFACE "${dtl_SOURCE_DIR}")
add_library(dtl::dtl ALIAS sequant_dtl)
