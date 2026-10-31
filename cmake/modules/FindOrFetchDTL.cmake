include(FetchContent)

if (NOT TARGET dtl::dtl)
    if (SEQUANT_FETCH_MISSING_DEPENDENCIES)
        FetchContent_Declare(
            DTL
            GIT_REPOSITORY https://github.com/cubicdaiya/dtl
            GIT_TAG v1.21
            GIT_SHALLOW TRUE
            EXCLUDE_FROM_ALL
        )

        FetchContent_MakeAvailable(DTL)
        set(_dtl_include_dir "${dtl_SOURCE_DIR}")
    else()
        # dtl has no CMake package; tests/ needs it even if SEQUANT_TESTS=OFF
        find_path(DTL_INCLUDE_DIR NAMES dtl/dtl.hpp DOC "directory containing dtl/dtl.hpp")
        if (NOT DTL_INCLUDE_DIR OR NOT EXISTS "${DTL_INCLUDE_DIR}/dtl/dtl.hpp")
            message(FATAL_ERROR "dtl (https://github.com/cubicdaiya/dtl) not found: set DTL_INCLUDE_DIR to the directory containing dtl/dtl.hpp, or set SEQUANT_FETCH_MISSING_DEPENDENCIES=ON")
        endif()
        set(_dtl_include_dir "${DTL_INCLUDE_DIR}")
    endif()

    add_library(sequant_dtl INTERFACE)
    target_include_directories(sequant_dtl SYSTEM INTERFACE "${_dtl_include_dir}")
    add_library(dtl::dtl ALIAS sequant_dtl)
endif()

# postcond check
if (NOT TARGET dtl::dtl)
	message(FATAL_ERROR "FindOrFetchDTL could not make TARGET dtl::dtl available")
endif()
