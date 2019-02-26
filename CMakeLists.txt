cmake_minimum_required(VERSION 3.9)
project(SeQuant2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(Range-V3 INTERFACE IMPORTED)
set_property(TARGET Range-V3 PROPERTY
             INTERFACE_INCLUDE_DIRECTORIES ${RANGEV3_DIR}/include)

# need Boost.CallableTraits
find_package(Boost 1.68 REQUIRED)
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

add_library(SeQuant2
        src/SeQuant2/sequant.cpp
        src/SeQuant2/sequant.hpp
        src/SeQuant2/attr.hpp
        src/SeQuant2/index.hpp
        src/SeQuant2/space.cpp
        src/SeQuant2/space.hpp
        src/SeQuant2/op.cpp
        src/SeQuant2/op.hpp
        src/SeQuant2/tensor.cpp
        src/SeQuant2/tensor.hpp
        src/SeQuant2/expr.cpp
        src/SeQuant2/expr.hpp
        src/SeQuant2/wick.hpp
        src/SeQuant2/ranges.hpp
        src/SeQuant2/container.hpp
        src/SeQuant2/meta.hpp
        src/SeQuant2/latex.hpp
        src/SeQuant2/wolfram.hpp
        src/SeQuant2/expr_algorithm.hpp
        src/SeQuant2/wick.impl.hpp
        src/SeQuant2/expr_operator.hpp
        src/SeQuant2/hash.hpp
        src/SeQuant2/tag.hpp
        src/SeQuant2/hugenholtz.hpp
        src/SeQuant2/algorithm.hpp
        src/SeQuant2/mbpt/sr/sr.cpp
        src/SeQuant2/mbpt/sr/sr.hpp
        src/SeQuant2/tensor_network.cpp)
target_link_libraries(SeQuant2 Range-V3 Boost bliss)

enable_testing(true)
add_custom_target(check_unit USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} -V)

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
        tests/unit/test_mbpt.cpp
        tests/unit/timer.hpp)

#set(utests_src
#        tests/unit/test_wick.cpp
#        tests/unit/timer.hpp)

set(unit_test_executable unit_tests)
add_executable(${unit_test_executable} EXCLUDE_FROM_ALL
        tests/unit/test_main.cpp
        tests/unit/catch.hpp
        ${utests_src})
target_link_libraries(${unit_test_executable} SeQuant2)
add_test(unit/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${unit_test_executable})

add_test(NAME unit/run
        COMMAND ${unit_test_executable})
set_tests_properties(unit/run
        PROPERTIES DEPENDS unit/build
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/unit)

####### Tests ########

# PNO-CCSD
set(example1 pno-ccsd)
add_executable(${example1} EXCLUDE_FROM_ALL
        examples/${example1}.cpp)
target_link_libraries(${example1} SeQuant2)

# Traces -- disabled, will remove if bliss is sufficient
#add_library(Traces INTERFACE IMPORTED)
#set_property(TARGET Traces PROPERTY
#        INTERFACE_INCLUDE_DIRECTORIES /Users/evaleev/code/build/nauty27rc1)
#set_property(TARGET Traces PROPERTY
#        PROPERTY INTERFACE_LINK_LIBRARIES /Users/evaleev/code/build/nauty27rc1/nauty.a)
#
#set(example2 traces-test)
#add_executable(${example2} EXCLUDE_FROM_ALL
#        examples/${example2}.cpp)
#target_link_libraries(${example2} Traces)

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
