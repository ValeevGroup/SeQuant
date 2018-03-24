cmake_minimum_required(VERSION 3.9)
project(SeQuant2)

set(CMAKE_CXX_STANDARD 17)

add_library(SeQuant2
        src/SeQuant2/sequant.cpp
        src/SeQuant2/sequant.hpp
        src/SeQuant2/index.hpp
        src/SeQuant2/space.hpp
        src/SeQuant2/op.hpp
        src/SeQuant2/vacuum.hpp
        src/SeQuant2/tensor.hpp
        src/SeQuant2/expr.hpp
        src/SeQuant2/wick.hpp
        )

enable_testing(true)
add_custom_target(check_unit COMMAND ${CMAKE_CTEST_COMMAND} -V)

set(utests_src
        tests/unit/test_space.cpp
        tests/unit/test_index.cpp
        tests/unit/test_op.cpp
        tests/unit/test_wick.cpp
        tests/unit/test_tensor.cpp tests/unit/test_expr.cpp)

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
