
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

project(SeQuant LANGUAGES CXX)
set(SEQUANT_VERSION 2.0.0-alpha.1)

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
include(external/versions.cmake)

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
if (NOT RANGEV3_DIR)
    message(FATAL_ERROR "Set cache variable RANGEV3_DIR "
        "to the top of the Range-V3 source directory")
endif(NOT RANGEV3_DIR)
if (NOT EXISTS "${RANGEV3_DIR}/include")
    message(FATAL_ERROR "Given RANGEV3_DIR=\"${RANGEV3_DIR}\" but"
        " did not find include subdirectory;"
        " make sure that RANGEV3_DIR points"
        " to the top of valid Range-V3 source directory")
endif()
# if Range-V3 is a git repo check the revision
if(EXISTS "${RANGEV3_DIR}/.git")
  find_package(Git REQUIRED)
  execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse -q HEAD
      WORKING_DIRECTORY "${RANGEV3_DIR}"
      OUTPUT_VARIABLE RANGEV3_REVISION )
  string(REGEX MATCH "[0-9a-f]*"
      RANGEV3_REVISION "${RANGEV3_REVISION}")
  if (NOT (${RANGEV3_REVISION} STREQUAL ${SEQUANT_TRACKED_RANGESV3_TAG}))
    message(FATAL_ERROR "RANGEV3_DIR=\"${RANGEV3_DIR}\" points to Range-V3"
        + " library with revision ${RANGEV3_REVISION}, "
        + " but need ${SEQUANT_TRACKED_RANGESV3_TAG}")
  endif ()
endif()
add_library(Range-V3 INTERFACE IMPORTED)
set_property(TARGET Range-V3 PROPERTY
             INTERFACE_INCLUDE_DIRECTORIES "${RANGEV3_DIR}/include")

# need Boost.ContainerHash
# NB Boost.Container is broken in v1.70
find_package(Boost ${SEQUANT_TRACKED_BOOST_VERSION} REQUIRED COMPONENTS serialization)
if (Boost_VERSION EQUAL 107000)
    message(FATAL_ERROR "Found Boost 1.70, but its Boost.Container is broken: "
        "use any version ${SEQUANT_TRACKED_BOOST_VERSION} "
        "or higher, other than 1.70")
endif()

# embedded bliss-0.73
add_library(bliss
        src/SeQuant/external/bliss/defs.cc
        src/SeQuant/external/bliss/defs.hh
        src/SeQuant/external/bliss/graph.cc
        src/SeQuant/external/bliss/graph.hh
        src/SeQuant/external/bliss/partition.cc
        src/SeQuant/external/bliss/partition.hh
        src/SeQuant/external/bliss/orbit.cc
        src/SeQuant/external/bliss/orbit.hh
        src/SeQuant/external/bliss/uintseqhash.cc
        src/SeQuant/external/bliss/uintseqhash.hh
        src/SeQuant/external/bliss/heap.cc
        src/SeQuant/external/bliss/heap.hh
        src/SeQuant/external/bliss/timer.cc
        src/SeQuant/external/bliss/timer.hh
        src/SeQuant/external/bliss/utils.cc
        src/SeQuant/external/bliss/utils.hh
        src/SeQuant/external/bliss/bliss_C.cc
        src/SeQuant/external/bliss/bliss_C.h
        )

add_library(SeQuant
        src/SeQuant/core/sequant.cpp
        src/SeQuant/core/sequant.hpp
        src/SeQuant/core/attr.hpp
        src/SeQuant/core/index.cpp
        src/SeQuant/core/index.hpp
        src/SeQuant/core/space.cpp
        src/SeQuant/core/space.hpp
        src/SeQuant/core/op.cpp
        src/SeQuant/core/op.hpp
        src/SeQuant/core/abstract_tensor.hpp
        src/SeQuant/core/abstract_tensor.cpp
        src/SeQuant/core/tensor.cpp
        src/SeQuant/core/tensor.hpp
        src/SeQuant/core/expr.cpp
        src/SeQuant/core/expr.hpp
        src/SeQuant/core/wick.hpp
        src/SeQuant/core/ranges.hpp
        src/SeQuant/core/container.hpp
        src/SeQuant/core/meta.hpp
        src/SeQuant/core/latex.hpp
        src/SeQuant/core/wolfram.hpp
        src/SeQuant/core/expr_algorithm.hpp
        src/SeQuant/core/wick.impl.hpp
        src/SeQuant/core/expr_operator.hpp
        src/SeQuant/core/hash.hpp
        src/SeQuant/core/tag.hpp
        src/SeQuant/core/any.hpp
        src/SeQuant/core/hugenholtz.hpp
        src/SeQuant/core/algorithm.hpp
        src/SeQuant/core/tensor_network.cpp
        src/SeQuant/core/tensor_network.hpp
        src/SeQuant/core/runtime.hpp
        src/SeQuant/core/utility.hpp
        src/SeQuant/core/bliss.hpp
        src/SeQuant/core/timer.hpp
        src/SeQuant/domain/mbpt/convention.cpp
        src/SeQuant/domain/mbpt/convention.hpp
        src/SeQuant/domain/mbpt/op.cpp
        src/SeQuant/domain/mbpt/op.hpp
        src/SeQuant/domain/mbpt/spin.hpp
        src/SeQuant/domain/mbpt/sr/sr.cpp
        src/SeQuant/domain/mbpt/sr/sr.hpp
        )
target_link_libraries(SeQuant Range-V3 Boost::boost bliss Threads::Threads)
if (PSTL_NEEDS_TBB)
    target_link_libraries(SeQuant ${TBB_LIBRARIES})
    target_include_directories(SeQuant PUBLIC ${TBB_INCLUDE_DIRS})
endif()
target_include_directories(SeQuant PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)

find_package(Libint2 2.6.0 CONFIG)

find_package(TiledArray CONFIG COMPONENTS tiledarray)

add_library(evaluator EXCLUDE_FROM_ALL
    src/SeQuant/domain/evaluate/eval_fwd.hpp
    src/SeQuant/domain/evaluate/eval_tree.hpp
    src/SeQuant/domain/evaluate/eval_tree.cpp
    src/SeQuant/domain/evaluate/eval_tree_node.hpp
    src/SeQuant/domain/evaluate/eval_tree_node.cpp
    src/SeQuant/domain/evaluate/factorizer.hpp
    src/SeQuant/domain/evaluate/factorizer.cpp
    # src/SeQuant/domain/evaluate/path_tree.hpp
    # src/SeQuant/domain/evaluate/path_tree.cpp
    )
target_include_directories(evaluator PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
target_link_libraries(evaluator SeQuant ${TBB_LIBRARIES})

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
        tests/unit/test_spin.cpp
        tests/unit/test_canonicalize.cpp
        # tests/unit/test_tensor_network.cpp
        tests/unit/test_tn_factorization.cpp
        )
set(utests_deps SeQuant evaluator)

if (TARGET tiledarray AND TARGET Libint2::cxx)
  list(APPEND utests_src
      ## unit tests for evaluation
      tests/unit/test_eval_tree.cpp
      # tests/unit/test_path_tree.cpp
      )
  list(APPEND utests_deps tiledarray)
endif(TARGET tiledarray AND TARGET Libint2::cxx)

set(unit_test_executable unit_tests)
add_executable(${unit_test_executable} EXCLUDE_FROM_ALL
        tests/unit/test_main.cpp
        tests/unit/catch.hpp
        ${utests_src})
target_link_libraries(${unit_test_executable} ${utests_deps})
add_test(unit/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR}
    --target ${unit_test_executable})

add_test(NAME unit/run
        COMMAND ${unit_test_executable})
set_tests_properties(unit/run
        PROPERTIES DEPENDS unit/build
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/unit)

####### Tests ########

# Single-Reference Coupled-Cluster
set(example0 srcc)
add_executable(${example0} EXCLUDE_FROM_ALL
        examples/${example0}/${example0}.cpp)
target_link_libraries(${example0} SeQuant ${TBB_LIBRARIES})

if (TARGET tiledarray AND TARGET Libint2::cxx)
  # sequant expression evaluation example
  set(example1 expr_evaluation)
  add_executable(${example1} EXCLUDE_FROM_ALL
      examples/evaluation/${example1}.cpp
      examples/contract/interpret/contract.hpp
      examples/contract/scf/hartree-fock.cc)
  target_link_libraries(${example1} SeQuant evaluator ${TBB_LIBRARIES}
      Libint2::cxx tiledarray)
  target_link_directories(${example1} PUBLIC ${TBB_LIBRARIES})
  target_compile_definitions(${example1} PRIVATE SEQUANT_HAS_TILEDARRAY)

  set(example2 rcc_evaluation)
  add_executable(${example2} EXCLUDE_FROM_ALL
      examples/evaluation/${example2}.cpp
      examples/contract/interpret/contract.hpp
      examples/contract/scf/hartree-fock.cc)
  target_link_libraries(${example2} SeQuant evaluator
      ${TBB_LIBRARIES} Libint2::cxx tiledarray)
  target_link_directories(${example2} PUBLIC ${TBB_LIBRARIES})
  target_compile_definitions(${example2} PRIVATE SEQUANT_HAS_TILEDARRAY)
endif (TARGET tiledarray AND TARGET Libint2::cxx)

set(example3 uccf12)
add_executable(${example3} EXCLUDE_FROM_ALL
    examples/${example3}/${example3}.cpp)
target_link_libraries(${example3} SeQuant ${TBB_LIBRARIES})

# add tests for running examples
set (lastexample 4)
foreach(i RANGE ${lastexample})
    if (TARGET ${example${i}})
        add_test(example/${example${i}}/build "${CMAKE_COMMAND}"
            --build ${CMAKE_BINARY_DIR} --target ${example${i}})
        add_test(NAME example/${example${i}}/run
                COMMAND ${example${i}}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples)
        set_tests_properties(example/${example${i}}/run
                PROPERTIES DEPENDS example/${example${i}}/build)
    endif()
endforeach()

# target for running ALL tests
add_custom_target_subproject(sequant check
    USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} -V)

####### DOCS ########

add_subdirectory(doc)
