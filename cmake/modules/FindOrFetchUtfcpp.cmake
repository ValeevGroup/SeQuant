if (NOT TARGET utfcpp)
    include(FetchContent)

    FetchContent_Declare(
        utfcpp
        GIT_REPOSITORY "https://github.com/nemtrif/utfcpp.git"
        GIT_TAG "v4.0.6"
        GIT_SHALLOW
        # Setting SOURCE_SUBDIR to a non-existing directory is the suggested workaround
        # to prevent FetchContent_MakeAvailable to add_subdirectory until
        # https://gitlab.kitware.com/cmake/cmake/-/issues/26220 gets implemented
        # See also https://discourse.cmake.org/t/prevent-fetchcontent-makeavailable-to-execute-cmakelists-txt/12704
        SOURCE_SUBDIR "_Don't use add_subdirectory_"
    )

    FetchContent_MakeAvailable(utfcpp)

    add_library(utfcpp INTERFACE)
    target_include_directories(utfcpp SYSTEM INTERFACE "${utfcpp_SOURCE_DIR}/..")
endif()

# postcond check
if (NOT TARGET utfcpp)
    message(FATAL_ERROR "FindOrFetchUtfcpp could not make TARGET utfcpp available")
endif()
