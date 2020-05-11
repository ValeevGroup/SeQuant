find_package(Range-v3 1.0.0 QUIET)

if (NOT TARGET range-v3)

  include(FetchContent)
  FetchContent_Declare(
      RangeV3
      GIT_REPOSITORY      https://github.com/ericniebler/range-v3.git
      GIT_TAG             ${SEQUANT_TRACKED_RANGEV3_TAG}
  )
  FetchContent_MakeAvailable(RangeV3)
  FetchContent_GetProperties(RangeV3
      SOURCE_DIR RANGEV3_SOURCE_DIR
      BINARY_DIR RANGEV3_BINARY_DIR
      )

endif(NOT TARGET range-v3)
