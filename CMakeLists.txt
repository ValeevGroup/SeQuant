cmake_minimum_required(VERSION 3.9)
project(SeQuant)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(Range-V3 INTERFACE IMPORTED)
set_property(TARGET Range-V3 PROPERTY
             INTERFACE_INCLUDE_DIRECTORIES ${RANGEV3_DIR}/include)

# need Boost.ContainerHash
find_package(Boost 1.67 REQUIRED)
add_library(Boost INTERFACE IMPORTED)
set_property(TARGET Boost PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})

# embedded bliss-0.73
add_library(bliss
        external/bliss/defs.cc
        external/bliss/defs.hh
        external/bliss/graph.cc
        external/bliss/graph.hh
        external/bliss/partition.cc
        external/bliss/partition.hh
        external/bliss/orbit.cc
        external/bliss/orbit.hh
        external/bliss/uintseqhash.cc
        external/bliss/uintseqhash.hh
        external/bliss/heap.cc
        external/bliss/heap.hh
        external/bliss/timer.cc
        external/bliss/timer.hh
        external/bliss/utils.cc
        external/bliss/utils.hh
        external/bliss/bliss_C.cc
        external/bliss/bliss_C.h
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
        src/SeQuant/domain/mbpt/convention.cpp
        src/SeQuant/domain/mbpt/convention.hpp
        src/SeQuant/domain/mbpt/op.cpp
        src/SeQuant/domain/mbpt/op.hpp
        src/SeQuant/domain/mbpt/spin.hpp
        src/SeQuant/domain/mbpt/sr/sr.cpp
        src/SeQuant/domain/mbpt/sr/sr.hpp
        src/SeQuant/core/timer.hpp)
target_link_libraries(SeQuant Range-V3 Boost bliss)
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
target_link_libraries(${example0} SeQuant)
add_test(example/${example0}/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${example0})
add_test(NAME example/${example0}/run
        COMMAND ${example0})
set_tests_properties(example/${example0}/run
        PROPERTIES DEPENDS example/${example0}/build)

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

