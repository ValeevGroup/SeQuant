project(python-sequant)
cmake_minimum_required(VERSION 3.11)

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/ValeevGroup/pybind11.git
    GIT_TAG dc7bb5db9a8688374bc16651fb81c964f9d62bf1 # ValeevGroup/v2.4.3
)
FetchContent_MakeAvailable(pybind11)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall)
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

enable_testing()

# add the executable
add_test(
  NAME test_sequant
  # need to use pytest to find sequant module properly
  COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test_sequant.py
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
