
cmake_minimum_required (VERSION 3.15.0) # need list(PREPEND for toolchains

###############################################################################
# Bring ValeevGroup cmake toolkit
###############################################################################
include(FetchContent)
FetchContent_Populate(
    vg_cmake_kit
    QUIET
    GIT_REPOSITORY      https://github.com/ValeevGroup/kit-cmake.git
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/cmake/vg
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/cmake/vg-build
    SUBBUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/cmake/vg-subbuild
)

# Set SeQuant version
set(SEQUANT_MAJOR_VERSION 2)
set(SEQUANT_MINOR_VERSION 0)
set(SEQUANT_MICRO_VERSION 0)
set(SEQUANT_PRERELEASE_ID alpha.1)
set(SEQUANT_VERSION "${SEQUANT_MAJOR_VERSION}.${SEQUANT_MINOR_VERSION}.${SEQUANT_MICRO_VERSION}")
if (SEQUANT_PRERELEASE_ID)
  set(SEQUANT_EXT_VERSION "${SEQUANT_VERSION}-${SEQUANT_PRERELEASE_ID}")
else(SEQUANT_PRERELEASE_ID)
  set(SEQUANT_EXT_VERSION "${SEQUANT_VERSION}")
endif(SEQUANT_PRERELEASE_ID)

# make SeQuant project
project(SeQuant LANGUAGES CXX VERSION "${SEQUANT_VERSION}" DESCRIPTION "SEcond QUANTization toolkit")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable ccache if not already enabled by symlink masquerading
if(NOT CMAKE_CXX_COMPILER MATCHES ".*/ccache$")
    find_program(CCACHE_EXECUTABLE ccache)
    if (CCACHE_EXECUTABLE)
        mark_as_advanced(CCACHE_EXECUTABLE)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ${CCACHE_EXECUTABLE})
        set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ${CCACHE_EXECUTABLE})
        message(STATUS "Found ccache: ${CCACHE_EXECUTABLE}")
    else()
        set(CCACHE_EXECUTABLE "")
    endif()
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules/")

set(FETCHCONTENT_UPDATES_DISCONNECTED OFF CACHE BOOL "Enables UPDATE_DISCONNECTED behavior for all content population")
include(FetchContent)
include(AddCustomTargetSubproject)
include(FeatureSummary)
include(CMakePackageConfigHelpers)
include(external/versions.cmake)

##########################
# Standard build variables
##########################
set(SEQUANT_INSTALL_BINDIR "bin"
    CACHE PATH "SeQuant BIN install directory")
set(SEQUANT_INSTALL_INCLUDEDIR "include"
    CACHE PATH "SeQuant INCLUDE install directory")
set(SEQUANT_INSTALL_LIBDIR "lib"
    CACHE PATH "SeQuant LIB install directory")
set(SEQUANT_INSTALL_SHAREDIR "share/sequant/${SEQUANT_MAJOR_VERSION}.${SEQUANT_MINOR_VERSION}.${SEQUANT_MICRO_VERSION}"
    CACHE PATH "SeQuant SHARE install directory")
set(SEQUANT_INSTALL_DATADIR "${SEQUANT_INSTALL_SHAREDIR}/data"
    CACHE PATH "SeQuant DATA install directory")
set(SEQUANT_INSTALL_DOCDIR "${SEQUANT_INSTALL_SHAREDIR}/doc"
    CACHE PATH "SeQuant DOC install directory")
set(SEQUANT_INSTALL_CMAKEDIR "lib/cmake/sequant"
    CACHE PATH "SeQuant CMAKE install directory")

# std::thread requires (on some platforms?) thread support which is not
# provided automatically
find_package(Threads REQUIRED)

# PSTL (used by g++-9 and clang++-8 in c++17 mode) needs TBB
if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" AND
    ${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 9) OR
    (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" AND
    ${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 8))
  set(PSTL_NEEDS_TBB TRUE)
  find_package(TBB REQUIRED)
endif()

# Ranges-V3
include(FindOrFetchRangeV3)

# need Boost.ContainerHash
# NB Boost.Container is broken in v1.70
find_package(Boost ${SEQUANT_TRACKED_BOOST_VERSION} CONFIG)
if (NOT TARGET Boost::boost)
  find_package(Boost ${SEQUANT_TRACKED_BOOST_VERSION} REQUIRED)
  set(Boost_USE_CONFIG FALSE)
  set(Boost_USE_CONFIG TRUE)
endif()
if (Boost_VERSION EQUAL 107000)
    message(FATAL_ERROR "Found Boost 1.70, but its Boost.Container is broken: "
        "use any version ${SEQUANT_TRACKED_BOOST_VERSION} "
        "or higher, other than 1.70")
endif()

# embedded bliss-0.73
add_library(SeQuant-bliss
        SeQuant/external/bliss/defs.cc
        SeQuant/external/bliss/defs.hh
        SeQuant/external/bliss/graph.cc
        SeQuant/external/bliss/graph.hh
        SeQuant/external/bliss/partition.cc
        SeQuant/external/bliss/partition.hh
        SeQuant/external/bliss/orbit.cc
        SeQuant/external/bliss/orbit.hh
        SeQuant/external/bliss/uintseqhash.cc
        SeQuant/external/bliss/uintseqhash.hh
        SeQuant/external/bliss/heap.cc
        SeQuant/external/bliss/heap.hh
        SeQuant/external/bliss/timer.cc
        SeQuant/external/bliss/timer.hh
        SeQuant/external/bliss/utils.cc
        SeQuant/external/bliss/utils.hh
        SeQuant/external/bliss/bliss_C.cc
        SeQuant/external/bliss/bliss_C.h
        )

add_library(SeQuant
        SeQuant/core/sequant.cpp
        SeQuant/core/sequant.hpp
        SeQuant/core/attr.hpp
        SeQuant/core/index.cpp
        SeQuant/core/index.hpp
        SeQuant/core/space.cpp
        SeQuant/core/space.hpp
        SeQuant/core/op.cpp
        SeQuant/core/op.hpp
        SeQuant/core/abstract_tensor.hpp
        SeQuant/core/abstract_tensor.cpp
        SeQuant/core/tensor.cpp
        SeQuant/core/tensor.hpp
        SeQuant/core/expr.cpp
        SeQuant/core/expr.hpp
        SeQuant/core/wick.hpp
        SeQuant/core/ranges.hpp
        SeQuant/core/container.hpp
        SeQuant/core/meta.hpp
        SeQuant/core/latex.hpp
        SeQuant/core/wolfram.hpp
        SeQuant/core/expr_algorithm.hpp
        SeQuant/core/wick.impl.hpp
        SeQuant/core/expr_operator.hpp
        SeQuant/core/hash.hpp
        SeQuant/core/tag.hpp
        SeQuant/core/any.hpp
        SeQuant/core/hugenholtz.hpp
        SeQuant/core/algorithm.hpp
        SeQuant/core/tensor_network.cpp
        SeQuant/core/tensor_network.hpp
        SeQuant/core/runtime.hpp
        SeQuant/core/utility.hpp
        SeQuant/core/bliss.hpp
        SeQuant/core/timer.hpp
        SeQuant/domain/evaluate/eval_fwd.hpp
        SeQuant/domain/evaluate/eval_tree.hpp
        SeQuant/domain/evaluate/eval_tree.cpp
        SeQuant/domain/evaluate/eval_tree_node.hpp
        SeQuant/domain/evaluate/eval_tree_node.cpp
        SeQuant/domain/factorize/factorizer.hpp
        SeQuant/domain/factorize/factorizer.cpp
        SeQuant/domain/evaluate/path_tree.hpp
        SeQuant/domain/evaluate/path_tree.cpp
        SeQuant/domain/mbpt/convention.cpp
        SeQuant/domain/mbpt/convention.hpp
        SeQuant/domain/mbpt/op.cpp
        SeQuant/domain/mbpt/op.hpp
        SeQuant/domain/mbpt/spin.hpp
        SeQuant/domain/mbpt/sr/sr.cpp
        SeQuant/domain/mbpt/sr/sr.hpp
        )
target_link_libraries(SeQuant range-v3::range-v3 Boost::boost SeQuant-bliss Threads::Threads)
if (PSTL_NEEDS_TBB)
    target_link_libraries(SeQuant ${TBB_LIBRARIES})
    target_include_directories(SeQuant PUBLIC ${TBB_INCLUDE_DIRS})
endif()
target_include_directories(SeQuant PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
target_compile_features(SeQuant INTERFACE "cxx_std_17")
install(TARGETS SeQuant-bliss EXPORT sequant COMPONENT sequant LIBRARY DESTINATION ${SEQUANT_INSTALL_LIBDIR})
install(TARGETS SeQuant EXPORT sequant COMPONENT sequant LIBRARY DESTINATION ${SEQUANT_INSTALL_LIBDIR})
install(DIRECTORY SeQuant
    COMPONENT sequant
    DESTINATION "${SEQUANT_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.hpp"
                   PATTERN "*.hh"
    )
add_library(SeQuant::SeQuant ALIAS SeQuant)  # to be able to use as subproject

### optional prereqs
if (NOT TARGET Libint2::libint2_cxx)
  find_package(Libint2 2.7.0 CONFIG)
endif(NOT TARGET Libint2::libint2_cxx)
if (NOT TARGET tiledarray)
  find_package(TiledArray CONFIG COMPONENTS tiledarray)
endif(NOT TARGET tiledarray)

enable_testing(true)

set(utests_src
         tests/unit/test_space.cpp
         tests/unit/test_index.cpp
         tests/unit/test_op.cpp
         tests/unit/test_wick.cpp
         tests/unit/test_tensor.cpp
         tests/unit/test_bliss.cpp
         tests/unit/test_expr.cpp
         tests/unit/test_iterator.cpp
         tests/unit/test_mbpt.cpp
         tests/unit/test_tensor_network.cpp
         tests/unit/test_tn_factorization.cpp
        )
set(utests_deps SeQuant)

if (TARGET tiledarray)
  list(APPEND utests_src
      ## unit tests for evaluation
      tests/unit/test_eval_tree.cpp
      tests/unit/test_path_tree.cpp
      )
  list(APPEND utests_deps tiledarray)
endif(TARGET tiledarray)

set(unit_test_executable unit_tests)
add_executable(${unit_test_executable} EXCLUDE_FROM_ALL
        tests/unit/test_main.cpp
        tests/unit/catch.hpp
        ${utests_src})
target_link_libraries(${unit_test_executable} ${utests_deps})
if (TARGET tiledarray)
  target_compile_definitions(${unit_test_executable} PRIVATE SEQUANT_HAS_TILEDARRAY)
endif(TARGET tiledarray)
add_test(unit/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR}
    --target ${unit_test_executable})

add_test(NAME unit/run
        COMMAND ${unit_test_executable})
set_tests_properties(unit/run
        PROPERTIES DEPENDS unit/build
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/unit)

####### Tests ########

# Single-Reference Coupled-Cluster equation generation (spin-orbital)
set(example0 srcc)
add_executable(${example0} EXCLUDE_FROM_ALL
        examples/${example0}/${example0}.cpp)
target_link_libraries(${example0} SeQuant ${TBB_LIBRARIES})

if (TARGET tiledarray AND TARGET Libint2::libint2_cxx)
  # SRCC spin-orbital equation evaluation
  set(example1 expr_evaluation)
  add_executable(${example1} EXCLUDE_FROM_ALL
      examples/evaluation/${example1}.cpp
      examples/contract/interpret/contract.hpp
      examples/contract/scf/hartree-fock.cc)
  target_link_libraries(${example1} SeQuant ${TBB_LIBRARIES}
      Libint2::libint2_cxx tiledarray)
  target_link_directories(${example1} PUBLIC ${TBB_LIBRARIES})
  target_compile_definitions(${example1} PRIVATE SEQUANT_HAS_TILEDARRAY)

  # SRCC spin-traced equation evaluation
  set(example2 rcc_evaluation)
  add_executable(${example2} EXCLUDE_FROM_ALL
      examples/evaluation/${example2}.cpp
      examples/contract/interpret/contract.hpp
      examples/contract/scf/hartree-fock.cc)
  target_link_libraries(${example2} SeQuant
      ${TBB_LIBRARIES} Libint2::libint2_cxx tiledarray)
  target_link_directories(${example2} PUBLIC ${TBB_LIBRARIES})
  target_compile_definitions(${example2} PRIVATE SEQUANT_HAS_TILEDARRAY)
endif (TARGET tiledarray AND TARGET Libint2::libint2_cxx)

set(example3 uccf12)
add_executable(${example3} EXCLUDE_FROM_ALL
    examples/${example3}/${example3}.cpp)
target_link_libraries(${example3} SeQuant ${TBB_LIBRARIES})

# add tests for running examples
set (lastexample 4)
foreach(i RANGE ${lastexample})
    if (TARGET ${example${i}})
        add_test(example/${example${i}}/build "${CMAKE_COMMAND}"
            --build ${PROJECT_BINARY_DIR} --target ${example${i}})
        add_test(NAME example/${example${i}}/run
                COMMAND ${example${i}}
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples)
        set_tests_properties(example/${example${i}}/run
                PROPERTIES DEPENDS example/${example${i}}/build)
    endif()
endforeach()

# target for running ALL tests
add_custom_target_subproject(sequant check
    USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} -V)

####### DOCS ########

add_subdirectory(doc)


##########################
# export SeQuant
##########################

if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
  find_package(Git REQUIRED)
  execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse -q HEAD
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE SEQUANT_REVISION )
  string(REGEX MATCH "[0-9a-f]*"
      SEQUANT_REVISION "${SEQUANT_REVISION}")
else()
  set(SEQUANT_REVISION "unknown")
endif()

configure_file(
    ${PROJECT_SOURCE_DIR}/SeQuant/version.hpp.in
    ${PROJECT_BINARY_DIR}/SeQuant/version.hpp
)
install(FILES ${PROJECT_BINARY_DIR}/SeQuant/version.hpp
    DESTINATION "${SEQUANT_INSTALL_INCLUDEDIR}/SeQuant")

# Create the version file
write_basic_package_version_file(sequant-config-version.cmake
    VERSION ${SEQUANT_VERSION} COMPATIBILITY AnyNewerVersion)

# Create the targets file
export(EXPORT sequant
    NAMESPACE SeQuant::
    FILE "${PROJECT_BINARY_DIR}/sequant-targets.cmake")

## Create the configure file
configure_package_config_file(cmake/sequant-config.cmake.in
    "${PROJECT_BINARY_DIR}/sequant-config.cmake"
    INSTALL_DESTINATION "${SEQUANT_INSTALL_CMAKEDIR}"
    PATH_VARS CMAKE_INSTALL_PREFIX SEQUANT_INSTALL_BINDIR
    SEQUANT_INSTALL_INCLUDEDIR SEQUANT_INSTALL_LIBDIR
    SEQUANT_INSTALL_DOCDIR SEQUANT_INSTALL_CMAKEDIR)

## Install config, version, and target files
install(EXPORT sequant
    FILE "sequant-targets.cmake"
    DESTINATION "${SEQUANT_INSTALL_CMAKEDIR}"
    NAMESPACE SeQuant::
    COMPONENT sequant-config)
install(FILES
    "${PROJECT_BINARY_DIR}/sequant-config.cmake"
    "${PROJECT_BINARY_DIR}/sequant-config-version.cmake"
    DESTINATION "${SEQUANT_INSTALL_CMAKEDIR}"
    COMPONENT sequant-config)
add_custom_target_subproject(sequant install-config
    COMMAND ${CMAKE_COMMAND} -DCOMPONENT=sequant-config -P ${PROJECT_BINARY_DIR}/cmake_install.cmake
    COMMENT "Installing SeQuant config components")

feature_summary(WHAT ALL
    DESCRIPTION "=== SeQuant Package/Feature Info ===")

###############################################################################
# appendix: misc details
###############################################################################
SET(CMAKE_COLOR_MAKEFILE ON)
