# -*- mode: cmake -*-

# Boost can be discovered by every (sub)package but only the top package can build it ...
# if we are the top package need to include the list of Boost components to be built
if("${CMAKE_PROJECT_NAME}" STREQUAL "${PROJECT_NAME}")
  set(required_components
          headers           # SeQuant TA, BTAS
          algorithm         # TA
          container         # SeQuant, TA, BTAS
          container_hash    # SeQuant
          core              # SeQuant
          filesystem
          fusion            # SeQuant
          iterator          # TA, BTAS
          locale            # SeQuant
          multiprecision    # SeQuant
#          numeric_conversion# SeQuant
#          numeric_interval  # SeQuant
          accumulators      # including numeric_* as targets is broken, instead include accumulators to ensure that numeric_{conversion,interval,ublas} are included
          random            # TA, BTAS
          range             # SeQuant
          regex             # SeQuant
          spirit            # SeQuant
          tuple             # TA
          variant           # SeQuant
  )
  if (DEFINED Boost_REQUIRED_COMPONENTS)
    list(APPEND Boost_REQUIRED_COMPONENTS
            ${required_components})
    list(REMOVE_DUPLICATES Boost_REQUIRED_COMPONENTS)
  else()
    set(Boost_REQUIRED_COMPONENTS "${required_components}" CACHE STRING "Components of Boost to discovered or built")
  endif()
  set(optional_components
          serialization     # BTAS
  )
  if (DEFINED Boost_OPTIONAL_COMPONENTS)
    list(APPEND Boost_OPTIONAL_COMPONENTS
            ${optional_components}
    )
    list(REMOVE_DUPLICATES Boost_OPTIONAL_COMPONENTS)
  else()
    set(Boost_OPTIONAL_COMPONENTS "${optional_components}" CACHE STRING "Optional components of Boost to discovered or built")
  endif()
endif()

if (NOT DEFINED Boost_FETCH_IF_MISSING)
  set(Boost_FETCH_IF_MISSING 1)
endif()

include(${vg_cmake_kit_SOURCE_DIR}/modules/FindOrFetchBoost.cmake)

# Boost.Move is broken in 1.77 and 1.78 unless using c++20
# fixed in 1.79 via https://github.com/boostorg/move/commit/78f26da1f3a5a3831e9e70efe83f9c56eef94e8c
if (CMAKE_CXX_STANDARD LESS 20)
  if (Boost_VERSION_MACRO GREATER_EQUAL 107700 AND Boost_VERSION_MACRO LESS 107900)
    message(FATAL_ERROR "Found Boost 1.77 <= version < 1.79, but its Boost.Move is broken with pre-C++20: use a version older than 1.77 or newer than 1.78")
  endif()
endif()
