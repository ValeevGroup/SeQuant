if (NOT TARGET libperm)
    include(${vg_cmake_kit_SOURCE_DIR}/modules/VRGFindOrFetchPackage.cmake)
    VRGFindOrFetchPackage(libperm "https://github.com/Krzmbrzl/libPerm.git" "${SEQUANT_TRACKED_LIBPERM_TAG}"
            ADD_SUBDIR
			ADD_SUBDIR_EXCLUDE_FROM_ALL
            CONFIG_SUBDIR
    )
endif()

# postcond check
if (NOT TARGET libperm)
    message(FATAL_ERROR "FindOrFetchLibPerm could not make TARGET libperm available")
endif(NOT TARGET libperm)
