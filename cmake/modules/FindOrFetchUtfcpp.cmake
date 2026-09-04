if (NOT TARGET utf8cpp::utf8cpp)
    include(FetchContent)

    FetchContent_Declare(
        utfcpp
        GIT_REPOSITORY "https://github.com/nemtrif/utfcpp.git"
		GIT_TAG "${SEQUANT_TRACKED_UTFCPP_TAG}"
        GIT_SHALLOW
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
