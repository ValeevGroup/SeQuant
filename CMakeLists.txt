cmake_minimum_required(VERSION 3.9)
project(SeQuant2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(Range-V3 INTERFACE IMPORTED)
set_property(TARGET Range-V3 PROPERTY
             INTERFACE_INCLUDE_DIRECTORIES ${RANGEV3_DIR}/include)

add_library(SeQuant2
        src/SeQuant2/sequant.cpp
        src/SeQuant2/sequant.hpp
        src/SeQuant2/index.hpp
        src/SeQuant2/space.cpp
        src/SeQuant2/space.hpp
        src/SeQuant2/op.hpp
        src/SeQuant2/vacuum.hpp
        src/SeQuant2/tensor.hpp
        src/SeQuant2/expr.hpp
        src/SeQuant2/wick.hpp
        src/SeQuant2/ranges.hpp
        src/SeQuant2/vector.hpp
        src/SeQuant2/meta.hpp
        src/SeQuant2/latex.hpp
        src/SeQuant2/expr_algorithm.hpp
        src/SeQuant2/wick.impl.hpp
        src/SeQuant2/expr_operator.hpp)
target_link_libraries(SeQuant2 Range-V3)

enable_testing(true)
add_custom_target(check_unit COMMAND ${CMAKE_CTEST_COMMAND} -V)

set(utests_src
        tests/unit/test_space.cpp
        tests/unit/test_index.cpp
        tests/unit/test_op.cpp
        tests/unit/test_wick.cpp
        tests/unit/test_tensor.cpp
        tests/unit/test_expr.cpp
        tests/unit/test_iterator.cpp
        tests/unit/timer.hpp)

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
