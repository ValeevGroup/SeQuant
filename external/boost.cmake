# -*- mode: cmake -*-

# Boost can be discovered by every (sub)package but only the top package can *build* it ...
# in either case must declare the components used by SeQuant
set(required_components
        headers
        container
        container_hash
        core
        fusion
        locale
        multiprecision
        #          numeric_conversion
        #          numeric_interval
        accumulators      # including numeric_* as targets is broken, instead include accumulators to ensure that numeric_{conversion,interval,ublas} are included
        range
        regex
        spirit
        variant
)
if (DEFINED Boost_REQUIRED_COMPONENTS)
    list(APPEND Boost_REQUIRED_COMPONENTS
            ${required_components})
    list(REMOVE_DUPLICATES Boost_REQUIRED_COMPONENTS)
else ()
    set(Boost_REQUIRED_COMPONENTS "${required_components}" CACHE STRING "Components of Boost to discovered or built")
endif ()
set(optional_components)
if (DEFINED Boost_OPTIONAL_COMPONENTS)
    list(APPEND Boost_OPTIONAL_COMPONENTS
            ${optional_components}
    )
    list(REMOVE_DUPLICATES Boost_OPTIONAL_COMPONENTS)
else ()
    set(Boost_OPTIONAL_COMPONENTS "${optional_components}" CACHE STRING "Optional components of Boost to discovered or built")
endif ()

if (NOT DEFINED Boost_FETCH_IF_MISSING)
    set(Boost_FETCH_IF_MISSING 1)
endif ()

include(${vg_cmake_kit_SOURCE_DIR}/modules/FindOrFetchBoost.cmake)

if (Boost_BUILT_FROM_SOURCE)
    if ("${PROJECT_NAME}" STREQUAL "${CMAKE_PROJECT_NAME}")
        add_custom_target(build-boost-in-SeQuant)
        foreach(_target IN LISTS Boost_FOUND_TARGETS Boost_MODULAR_TARGETS_NOT_BUILT_BY_INSTALL)
            add_dependencies(build-boost-in-SeQuant Boost::${_target})
        endforeach()
    endif()
endif()

# Boost.Move is broken in 1.77 and 1.78 unless using c++20
# fixed in 1.79 via https://github.com/boostorg/move/commit/78f26da1f3a5a3831e9e70efe83f9c56eef94e8c
if (CMAKE_CXX_STANDARD LESS 20)
    if (Boost_VERSION_MACRO GREATER_EQUAL 107700 AND Boost_VERSION_MACRO LESS 107900)
        message(FATAL_ERROR "Found Boost 1.77 <= version < 1.79, but its Boost.Move is broken with pre-C++20: use a version older than 1.77 or newer than 1.78")
    endif ()
endif ()
