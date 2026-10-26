if (NOT TARGET polymorphic_variant::polymorphic_variant)
    include(FetchContent)

    sequant_find_unless_fetching(polymorphic_variant ${SEQUANT_OLDEST_POLYMORPHICVARIANT_VERSION})

    # polymorphic_variant has no option to disable its install rules
    set(_exclude)
    if (NOT SEQUANT_INSTALL_FETCHED_DEPENDENCIES)
        set(_exclude EXCLUDE_FROM_ALL)
    endif()

    FetchContent_Declare(
        polymorphic_variant
        GIT_REPOSITORY "https://github.com/Krzmbrzl/polymorphic_variant.git"
        GIT_TAG "${SEQUANT_TRACKED_POLYMORPHICVARIANT_TAG}"
        GIT_SHALLOW
        SYSTEM
        ${_exclude}
        FIND_PACKAGE_ARGS ${SEQUANT_OLDEST_POLYMORPHICVARIANT_VERSION} NAMES polymorphic_variant
    )

    FetchContent_MakeAvailable(polymorphic_variant)

    if (polymorphic_variant_FOUND)
        set(SEQUANT_FETCHED_POLYMORPHICVARIANT FALSE)
    else()
        set(SEQUANT_FETCHED_POLYMORPHICVARIANT TRUE)
        if (SEQUANT_INSTALL_FETCHED_DEPENDENCIES)
            get_directory_property(_pv_version DIRECTORY ${polymorphic_variant_SOURCE_DIR} DEFINITION polymorphic_variant_VERSION)
            get_directory_property(_pv_cmakepath DIRECTORY ${polymorphic_variant_SOURCE_DIR} DEFINITION PV_INSTALL_CMAKEPATH)
            set(polymorphic_variant_CONFIG "${CMAKE_INSTALL_PREFIX}/${_pv_cmakepath}/polymorphic_variant-config.cmake")
            sequant_record_installed_dependency(polymorphic_variant VERSION ${_pv_version})
        else()
            unset(polymorphic_variant_CONFIG)
        endif()
    endif()
endif()

# postcond check
if (NOT TARGET polymorphic_variant::polymorphic_variant)
    message(FATAL_ERROR "FindOrFetchPolymorphicVariant could not make TARGET polymorphic_variant::polymorphic_variant available")
endif()
