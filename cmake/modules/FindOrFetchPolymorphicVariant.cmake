if (NOT TARGET polymorphic_variant::polymorphic_variant)
    include(FetchContent)

    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY "https://github.com/Krzmbrzl/polymorphic_variant.git"
        GIT_TAG "${SEQUANT_TRACKED_POLYMORPHICVARIANT_TAG}"
        GIT_SHALLOW
        EXCLUDE_FROM_ALL
        SYSTEM
    )
endif()

# postcond check
if (NOT TARGET polymorphic_variant::polymorphic_variant)
    message(FATAL_ERROR "FindOrFetchLibPerm could not make TARGET polymorphic_variant::polymorphic_variant available")
endif()
