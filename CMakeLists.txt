cmake_minimum_required(VERSION 3.9)
project(SeQuant LANGUAGES CXX)

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

# std::thread requires (on some platforms?) thread support which is not provided automatically
find_package(Threads REQUIRED)

# PSTL (used by g++-9 and clang++-8 in c++17 mode) needs TBB
if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 9) OR
    (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 8))
  set(PSTL_NEEDS_TBB TRUE)
  find_package(TBB REQUIRED)
endif()

if (NOT RANGEV3_DIR)
    message(FATAL_ERROR "Set cache variable RANGEV3_DIR to the top of the Range-V3 source directory")
endif(NOT RANGEV3_DIR)
add_library(Range-V3 INTERFACE IMPORTED)
set_property(TARGET Range-V3 PROPERTY
             INTERFACE_INCLUDE_DIRECTORIES ${RANGEV3_DIR}/include)

# need Boost.ContainerHash
# NB Boost.Container is broken in v1.70
find_package(Boost 1.67 REQUIRED)
if (Boost_VERSION EQUAL 107000)
    message(FATAL_ERROR "Found Boost 1.70, but its Boost.Container is broken: use versions 1.67, 1.68, or 1.69")
endif()
add_library(Boost INTERFACE IMPORTED)
set_property(TARGET Boost PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})

# look for BTAS, if needed
if (BTAS_INSTALL_DIR)
    find_package(BTAS)
endif(BTAS_INSTALL_DIR)

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
target_link_libraries(SeQuant Range-V3 Boost bliss Threads::Threads)
if (PSTL_NEEDS_TBB)
    target_link_libraries(SeQuant ${TBB_LIBRARIES})
    target_include_directories(SeQuant PUBLIC ${TBB_INCLUDE_DIRS})
endif()
target_include_directories(SeQuant PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )

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
        tests/unit/test_tensor_network.cpp
        tests/unit/test_mbpt.cpp)

#set(utests_src
#        tests/unit/test_mbpt.cpp)

set(unit_test_executable unit_tests)
add_executable(${unit_test_executable} EXCLUDE_FROM_ALL
        tests/unit/test_main.cpp
        tests/unit/catch.hpp
        ${utests_src})
target_link_libraries(${unit_test_executable} SeQuant)
add_test(unit/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${unit_test_executable})

add_test(NAME unit/run
        COMMAND ${unit_test_executable})
set_tests_properties(unit/run
        PROPERTIES DEPENDS unit/build
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/unit)

####### Tests ########

# Single-Reference Coupled-Cluster
set(example0 srcc)
add_executable(${example0} EXCLUDE_FROM_ALL
        examples/${example0}.cpp)
target_link_libraries(${example0} SeQuant ${TBB_LIBRARIES})
add_test(example/${example0}/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${example0})
add_test(NAME example/${example0}/run
        COMMAND ${example0})
set_tests_properties(example/${example0}/run
        PROPERTIES DEPENDS example/${example0}/build)

# contract method testing -- Bimal
if (BTAS_FOUND)
  find_package(Libint2 2.6.0 CONFIG REQUIRED)
  set(example1 test_contract)
  add_executable(${example1} EXCLUDE_FROM_ALL
          examples/${example1}.cpp examples/contractable_tensor.hpp examples/contract.hpp examples/hartree-fock.cc)
  target_link_libraries(${example1} SeQuant ${TBB_LIBRARIES} Libint2::cxx BTAS)
endif(BTAS_FOUND)

# # sequant::contractable::Tensor tests -- Bimal
# set(example2 test_ctensor)
# add_executable(${example2} EXCLUDE_FROM_ALL
#     examples/${example2}.cpp examples/contractable_tensor.hpp examples/contract.hpp)
# target_link_libraries(${example2} SeQuant)

# target for running ALL tests
add_custom_target(check USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} -V)

####### DOCS ########

find_package(Doxygen)
if (DOXYGEN_FOUND)
    set(DOXYGEN_GENERATE_HTML YES)
    doxygen_add_docs(
            html
            ${PROJECT_SOURCE_DIR}/src
            DOXYGEN_EXCLUDE_PATTERNS "catch.hpp"
            COMMENT "Generate html dox"
    )
endif(DOXYGEN_FOUND)

