project(python-sequant)
cmake_minimum_required(VERSION 3.11)

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/ValeevGroup/pybind11.git
    GIT_TAG dc7bb5db9a8688374bc16651fb81c964f9d62bf1 # ValeevGroup/v2.4.3
)
FetchContent_MakeAvailable(pybind11)

if (NOT(CMAKE_CXX_STANDARD GREATER_EQUAL 17))
  message(FATAL_ERROR "C++ 2017 ISO Standard or higher is required to compile PySeQuant")
endif()

# Python2.7 has register in Python.h
add_compile_options(-Wno-register)

pybind11_add_module(python-sequant MODULE src/sequant/_sequant.cc)

target_link_libraries(python-sequant PRIVATE SeQuant)

set_target_properties(
  python-sequant
  PROPERTIES
  #PREFIX ""
  OUTPUT_NAME _sequant
  )

if (BUILD_TESTING)
  add_test(sequant/unit/python/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target python-sequant)
  set_tests_properties(sequant/unit/python/build PROPERTIES FIXTURES_SETUP PYSEQUANT_MODULE)

  set(SEQUANT_TEST_PYTHON_ENVVARS "PYTHONPATH=${PROJECT_BINARY_DIR}:$ENV{PYTHONPATH}")
  if (DEFINED SEQUANT_TEST_EXTRA_ENVVARS)
    list(APPEND SEQUANT_TEST_PYTHON_ENVVARS ${SEQUANT_TEST_EXTRA_ENVVARS})
  endif()

  # basic test
  set(test_name sequant/unit/python/basic/run)
  add_test(NAME ${test_name}
           COMMAND ${PYTHON_EXECUTABLE} test_sequant.py
  )
  set_tests_properties(${test_name} PROPERTIES
           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
           FIXTURES_REQUIRED PYSEQUANT_MODULE
           ENVIRONMENT "${SEQUANT_TEST_PYTHON_ENVVARS}")

  unset(MPQC_TEST_PYTHON_ENVVARS)

  add_dependencies(check-sequant python-sequant)

endif (BUILD_TESTING)
