if (NOT TARGET libperm::libperm)
    include(FetchContent)

    # libPerm's own tests are never run by check-sequant (ctest -R "^sequant"),
    # and gtest_discover_tests runs the test binary at build time with a 5 s
    # timeout that a loaded CI runner can exceed
    set(LIBPERM_TESTS OFF CACHE BOOL "Whether to build libPerm's own tests" FORCE)

    sequant_find_unless_fetching(libperm ${SEQUANT_OLDEST_LIBPERM_VERSION})

    # keep libPerm's export set (the static SeQuant libraries link it privately), but leave it out of the install
    set(_exclude)
    if (NOT SEQUANT_INSTALL_FETCHED_DEPENDENCIES)
        set(_exclude EXCLUDE_FROM_ALL)
    endif()

    FetchContent_Declare(
        libperm
        GIT_REPOSITORY "https://github.com/Krzmbrzl/libPerm.git"
        GIT_TAG "${SEQUANT_TRACKED_LIBPERM_TAG}"
        GIT_SHALLOW
        SYSTEM
        ${_exclude}
        FIND_PACKAGE_ARGS ${SEQUANT_OLDEST_LIBPERM_VERSION}
    )

    FetchContent_MakeAvailable(libperm)

    if (libperm_FOUND)
        set(SEQUANT_FETCHED_LIBPERM FALSE)
    else()
        set(SEQUANT_FETCHED_LIBPERM TRUE)
        if (SEQUANT_INSTALL_FETCHED_DEPENDENCIES)
            # libPerm's project name has a trailing comma, so its <name>_VERSION is not usable
            get_directory_property(_libperm_version DIRECTORY ${libperm_SOURCE_DIR} DEFINITION PROJECT_VERSION)
            set(libperm_CONFIG "${CMAKE_INSTALL_PREFIX}/${LIBPERM_INSTALL_CMAKEDIR}/libperm-config.cmake")
            sequant_record_installed_dependency(libperm VERSION ${_libperm_version})
        else()
            unset(libperm_CONFIG)
        endif()
    endif()
endif()

# postcond check
if (NOT TARGET libperm::libperm)
    message(FATAL_ERROR "FindOrFetchLibPerm could not make TARGET libperm::libperm available")
endif()
