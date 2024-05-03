if (NOT TARGET Catch2::Catch2)

    include(${vg_cmake_kit_SOURCE_DIR}/modules/VRGFindOrFetchPackage.cmake)
    VRGFindOrFetchPackage(Catch2 "https://github.com/catchorg/Catch2.git" "${SEQUANT_TRACKED_CATCH2_TAG}"
            ADD_SUBDIR
            CONFIG_SUBDIR
    )
    if (TARGET Catch2 AND NOT TARGET Catch2::Catch2)
        add_library(Catch2::Catch2 ALIAS Catch2)
    endif()
endif()

# postcond check
if (NOT TARGET Catch2::Catch2)
    message(FATAL_ERROR "FindOrFetchCatch2 could not make TARGET Catch2 available")
endif(NOT TARGET Catch2::Catch2)
