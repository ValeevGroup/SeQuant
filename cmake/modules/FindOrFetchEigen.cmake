include(FetchContent)

if (NOT TARGET Eigen3::Eigen)
    sequant_find_unless_fetching(Eigen3 NO_CMAKE_PACKAGE_REGISTRY)

    set(EIGEN_BUILD_CMAKE_PACKAGE ON)
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG ${SEQUANT_TRACKED_EIGEN_TAG}
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL
        SYSTEM
        # re:NO_CMAKE_PACKAGE_REGISTRY: eigen3 registers its *build* tree with the user package registry ...
        #                               to avoid issues with wiped build directory look for installed eigen
        # re: no version argument: the version files of Eigen 3.x and 5.x disagree on ranges, so no single
        #                          version spec is accepted by both; the minimum is checked below instead
        FIND_PACKAGE_ARGS NO_CMAKE_PACKAGE_REGISTRY
    )

    FetchContent_MakeAvailable(Eigen3)

    if (Eigen3_FOUND AND Eigen3_VERSION VERSION_LESS SEQUANT_OLDEST_EIGEN_VERSION)
        message(FATAL_ERROR "Found Eigen3 ${Eigen3_VERSION} at ${Eigen3_DIR}, but SeQuant requires Eigen3 ${SEQUANT_OLDEST_EIGEN_VERSION} or later")
    endif()
endif()

# postcond check
if (NOT TARGET Eigen3::Eigen)
    message(FATAL_ERROR "FindOrFetchEigen could not make TARGET Eigen3::Eigen available")
endif()
