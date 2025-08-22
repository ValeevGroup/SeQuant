if (NOT TARGET polymorphic_variant::polymorphic_variant)
    include(${vg_cmake_kit_SOURCE_DIR}/modules/VRGFindOrFetchPackage.cmake)
	VRGFindOrFetchPackage(polymorphic_variant "https://github.com/Krzmbrzl/polymorphic_variant.git" "${SEQUANT_TRACKED_POLYMORPHICVARIANT_TAG}"
            ADD_SUBDIR
            CONFIG_SUBDIR
    )
endif()

# postcond check
if (NOT TARGET polymorphic_variant::polymorphic_variant)
    message(FATAL_ERROR "FindOrFetchLibPerm could not make TARGET polymorphic_variant::polymorphic_variant available")
endif()
