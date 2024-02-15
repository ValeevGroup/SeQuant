# -*- mode: cmake -*-

# update the Boost version that we can tolerate
if (NOT DEFINED Boost_OLDEST_BOOST_VERSION)
    set(Boost_OLDEST_BOOST_VERSION ${SEQUANT_OLDEST_BOOST_VERSION})
else()
    if (${Boost_OLDEST_BOOST_VERSION} VERSION_LESS ${SEQUANT_OLDEST_BOOST_VERSION})
        if (DEFINED CACHE{Boost_OLDEST_BOOST_VERSION})
            set(Boost_OLDEST_BOOST_VERSION "${SEQUANT_OLDEST_BOOST_VERSION}" CACHE STRING "Oldest Boost version to use" FORCE)
        else()
            set(Boost_OLDEST_BOOST_VERSION ${SEQUANT_OLDEST_BOOST_VERSION})
        endif()
    endif()
endif()

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
