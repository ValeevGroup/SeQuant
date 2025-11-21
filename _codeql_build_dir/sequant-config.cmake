# - CMAKE Config file for the SeQuant package
# This will define the following CMake cache variables
#
#    SEQUANT_FOUND        - true if SeQuant library were found
#    SEQUANT_VERSION      - the x.y.z SeQuant version
#    SEQUANT_EXT_VERSION  - the x.y.z(-buildid) SeQuant version, where the (optional) buildid is smth like beta.3
#
# and the following imported targets
#
#     SeQuant::SeQuant    - the SeQuant library
#

# Set package version
set(SEQUANT_VERSION "2.0.0")
set(SEQUANT_EXT_VERSION "2.0.0-alpha.1")


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was sequant-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################


# converts modular component _comp to list of targets defined by the component
# target name x means TARGET Boost::x is defined .. this is provided just in case on needs to map components to targets
macro(component_to_targets _comp _targets)
    if (${${_comp}} STREQUAL test)
        set(${_targets} unit_test_framework)
    else()
        set(${_targets} ${${_comp}})
    endif()
endmacro()

#########################################
# import boost components, if any missing
#########################################
set(Boost_IS_MODULARIZED ON)
if (Boost_IS_MODULARIZED)
  set(Boost_FOUND_COMPONENTS headers;container;container_hash;core;dynamic_bitset;fusion;hana;locale;multiprecision;accumulators;range;regex;spirit;unordered;variant)
  set(Boost_FOUND_TARGETS headers;container;container_hash;core;dynamic_bitset;fusion;hana;locale;multiprecision;accumulators;range;regex;spirit;unordered;variant)
else(Boost_IS_MODULARIZED)
  set(Boost_FOUND_COMPONENTS headers;locale;regex)
  set(Boost_FOUND_TARGETS headers;locale;regex)
endif(Boost_IS_MODULARIZED)

set(Boost_DEPS_LIBRARIES_NOT_FOUND_CHECK "NOT;TARGET;Boost::headers")
foreach(_tgt ${Boost_FOUND_TARGETS})
  list(APPEND Boost_DEPS_LIBRARIES_NOT_FOUND_CHECK "OR;NOT;TARGET;Boost::${_tgt}")
endforeach(_tgt)

if(${Boost_DEPS_LIBRARIES_NOT_FOUND_CHECK})
  include( CMakeFindDependencyMacro )
  set(Boost_BUILT_FROM_SOURCE ON)
  if (NOT Boost_BUILT_FROM_SOURCE)
    set(Boost_USE_CONFIG )
    # OPTIONAL_COMPONENTS in FindBoost available since 3.11
    cmake_minimum_required(VERSION 3.11.0)
    if (Boost_USE_CONFIG)
      set(Boost_CONFIG )
      if (NOT Boost_CONFIG OR NOT EXISTS ${Boost_CONFIG})
        message(FATAL_ERROR "Expected Boost config file at ${Boost_CONFIG}; directory moved since BTAS configuration?")
      endif()
      get_filename_component(Boost_DIR ${Boost_CONFIG} DIRECTORY)
      find_dependency(Boost REQUIRED OPTIONAL_COMPONENTS ${Boost_FOUND_COMPONENTS} PATHS ${Boost_DIR} NO_DEFAULT_PATH)
    else (Boost_USE_CONFIG)
      set(BOOST_INCLUDEDIR )
      set(BOOST_LIBRARYDIR )
      if (NOT BOOST_LIBRARYDIR OR NOT EXISTS ${BOOST_LIBRARYDIR})
        set(BOOST_LIBRARYDIR )
      endif()
      set(Boost_NO_SYSTEM_PATHS OFF)
      if (BOOST_LIBRARYDIR AND BOOST_INCLUDEDIR)
        if (EXISTS ${BOOST_LIBRARYDIR} AND EXISTS ${BOOST_INCLUDEDIR})
          set(Boost_NO_SYSTEM_PATHS ON)
        endif()
      endif()
      find_dependency(Boost REQUIRED OPTIONAL_COMPONENTS ${Boost_FOUND_COMPONENTS})
    endif (Boost_USE_CONFIG)
  else(NOT Boost_BUILT_FROM_SOURCE)
    foreach(_tgt IN LISTS Boost_FOUND_TARGETS)
      if (NOT TARGET Boost::${_tgt})
        find_dependency(boost_${_tgt} CONFIG REQUIRED)
      endif()
    endforeach(_tgt)
  endif(NOT Boost_BUILT_FROM_SOURCE)
endif(${Boost_DEPS_LIBRARIES_NOT_FOUND_CHECK})



set(SEQUANT_HAS_TILEDARRAY )
if(SEQUANT_HAS_TILEDARRAY AND NOT TARGET tiledarray)
  set(TiledArray_CONFIG )
  if (NOT TiledArray_CONFIG OR NOT EXISTS ${TiledArray_CONFIG})
    message(FATAL_ERROR "Expected TiledArray config file at ${TiledArray_CONFIG}; directory moved since SeQuant configuration?")
  endif()
  get_filename_component(TiledArray_DIR ${TiledArray_CONFIG} DIRECTORY)
  find_dependency(TiledArray CONFIG QUIET REQUIRED COMPONENTS tiledarray PATHS ${TiledArray_DIR} NO_DEFAULT_PATH)
endif()

set(SEQUANT_HAS_EIGEN )
if (NOT TARGET Eigen3::Eigen AND SEQUANT_HAS_EIGEN)
  if (TARGET TiledArray_Eigen)
     add_library(Eigen3::Eigen ALIAS TiledArray_Eigen)
  else()
    get_filename_component(Eigen3_DIR "" DIRECTORY)
     # re:NO_CMAKE_PACKAGE_REGISTRY: eigen3 registers its *build* tree with the user package registry ...
     #                               to avoid issues with wiped build directory look for installed eigen
     find_package(Eigen3 3.0 CONFIG REQUIRED NO_CMAKE_PACKAGE_REGISTRY HINTS "${Eigen3_DIR}")
  endif()
endif()

if (NOT TARGET range-v3::range-v3)
  get_filename_component(range-v3_DIR "/usr/local/lib/cmake/range-v3/range-v3-config.cmake" DIRECTORY)
  find_dependency(range-v3 QUIET REQUIRED HINTS "${range-v3_DIR}")
endif(NOT TARGET range-v3::range-v3)

if (NOT TARGET libperm::libperm)
  get_filename_component(libperm_DIR "" DIRECTORY)
  find_dependency(libperm QUIET REQUIRED HINTS "${libperm_DIR}")
endif(NOT TARGET libperm::libperm)

if (NOT TARGET polymorphic_variant::polymorphic_variant)
  get_filename_component(polymorphic_variant_DIR "" DIRECTORY)
  find_dependency(polymorphic_variant QUIET REQUIRED HINTS "${polymorphic_variant_DIR}")
endif(NOT TARGET polymorphic_variant::polymorphic_variant)

if (NOT TARGET Threads::Threads)
    find_dependency(Threads REQUIRED)
endif (NOT TARGET Threads::Threads)

# Include library IMPORT targets
if(NOT TARGET SeQuant::SeQuant)
  include("${CMAKE_CURRENT_LIST_DIR}/sequant-targets.cmake")
  if(NOT TARGET SeQuant::SeQuant)
    message(FATAL_ERROR "expected SeQuant::SeQuant among imported SeQuant targets")
  endif()
endif()

set(SEQUANT_FOUND TRUE)
