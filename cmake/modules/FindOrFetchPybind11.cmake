if (NOT TARGET pybind11::module)
    include(FetchContent)

    FetchContent_Declare(
        Pybind11
        GIT_REPOSITORY "https://github.com/pybind/pybind11.git"
		GIT_TAG "${SEQUANT_TRACKED_PYBIND11_TAG}"
        GIT_SHALLOW
        EXCLUDE_FROM_ALL
        SYSTEM
        FIND_PACKAGE_ARGS ${SEQUANT_OLDEST_PYBIND11_VERSION} NAMES pybind11
    )

    FetchContent_MakeAvailable(Pybind11)
endif()

# postcond check
if (NOT TARGET pybind11::module)
    message(FATAL_ERROR "FindOrFetchPybind11 could not make TARGET pybind11::module available")
endif()
