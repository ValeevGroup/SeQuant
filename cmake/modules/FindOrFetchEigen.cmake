include(FetchContent)

if (NOT TARGET Eigen3::Eigen)
    # use TA's Eigen, if available
    if (TARGET TiledArray_Eigen)
        add_library(Eigen3::Eigen ALIAS TiledArray_Eigen)
    else()
		FetchContent_Declare(
			Eigen3
			GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
			GIT_TAG ${SEQUANT_TRACKED_EIGEN_TAG}
			GIT_SHALLOW TRUE
			EXCLUDE_FROM_ALL
			SYSTEM
		    # re:NO_CMAKE_PACKAGE_REGISTRY: eigen3 registers its *build* tree with the user package registry ...
		    #                               to avoid issues with wiped build directory look for installed eigen
			FIND_PACKAGE_ARGS ${SEQUANT_OLDEST_EIGEN_VERSION} NO_MODULE NO_CMAKE_PACKAGE_REGISTRY
		)

		FetchContent_MakeAvailable(Eigen3)
	endif()
endif()

# postcond check
if (NOT TARGET Eigen3::Eigen)
	message(FATAL_ERROR "FindOrFetchEigen could not make TARGET Eigen3::Eigen available")
endif()
