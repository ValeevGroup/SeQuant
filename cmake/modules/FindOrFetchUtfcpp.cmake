if (NOT TARGET utf8cpp::utf8cpp)
    include(FetchContent)

    sequant_find_unless_fetching(utf8cpp ${SEQUANT_OLDEST_UTFCPP_VERSION})

    FetchContent_Declare(
        utfcpp
        GIT_REPOSITORY "https://github.com/nemtrif/utfcpp.git"
		GIT_TAG "${SEQUANT_TRACKED_UTFCPP_TAG}"
        GIT_SHALLOW
        SYSTEM
        EXCLUDE_FROM_ALL
        FIND_PACKAGE_ARGS ${SEQUANT_OLDEST_UTFCPP_VERSION} NAMES utf8cpp
    )

    FetchContent_MakeAvailable(utfcpp)

	if (NOT TARGET utf8cpp::utf8cpp)
		# https://github.com/nemtrif/utfcpp/pull/146
		add_library(utf8cpp::utf8cpp ALIAS utf8cpp)
	endif()
endif()

# postcond check
if (NOT TARGET utf8cpp::utf8cpp)
    message(FATAL_ERROR "FindOrFetchUtfcpp could not make TARGET utf8cpp::utf8cpp available")
endif()
