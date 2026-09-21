project(python-sequant)

if (NOT TARGET Python::Module)
  find_package(Python 3 COMPONENTS Interpreter Development REQUIRED)
endif()

include(FindOrFetchPybind11)

if (CMAKE_CXX_STANDARD VERSION_LESS 17)
  message(FATAL_ERROR "C++ 2017 ISO Standard or higher is required to compile PySeQuant")
endif()

pybind11_add_module(python-sequant MODULE src/sequant/_sequant.cc)

target_link_libraries(python-sequant PRIVATE SeQuant)
target_include_directories(python-sequant PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
# pybind11_add_module applies LTO to the module by itself (pybind11::lto), which is why the
# automatic default leaves it alone - but pybind11's envelope is narrower than ours (it skips
# RelWithDebInfo, and skips LTO entirely when CMAKE_INTERPROCEDURAL_OPTIMIZATION is defined).
# An explicit SEQUANT_LTO=ON must reach the module in any case: without "fat" objects the
# SeQuant archives are then bitcode-only and can only be linked with LTO enabled on the
# consumer too. Where both apply, the duplicate -flto is harmless.
if (SEQUANT_LTO)
  target_set_optimization_flags(python-sequant)
endif()

set_target_properties(
  python-sequant
  PROPERTIES
  #PREFIX ""
  OUTPUT_NAME _sequant
  )

if (BUILD_TESTING)
  add_test(sequant/unit/python/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target python-sequant)
  set_tests_properties(sequant/unit/python/build PROPERTIES FIXTURES_SETUP PYSEQUANT_MODULE)

  if (WIN32)
    set(SEQUANT_PYTHONPATH_SEP "$<SEMICOLON>")
  else()
    set(SEQUANT_PYTHONPATH_SEP ":")
  endif()
  set(SEQUANT_TEST_PYTHON_ENVVARS "PYTHONPATH=${PROJECT_BINARY_DIR}${SEQUANT_PYTHONPATH_SEP}$ENV{PYTHONPATH}")
  unset(SEQUANT_PYTHONPATH_SEP)
  if (DEFINED SEQUANT_TEST_EXTRA_ENVVARS)
    list(APPEND SEQUANT_TEST_PYTHON_ENVVARS ${SEQUANT_TEST_EXTRA_ENVVARS})
  endif()

  # basic test
  set(test_name sequant/unit/python/basic/run)
  add_test(NAME ${test_name}
           COMMAND ${Python_EXECUTABLE} test_sequant.py
  )
  set_tests_properties(${test_name} PROPERTIES
           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
           FIXTURES_REQUIRED PYSEQUANT_MODULE
           ENVIRONMENT "${SEQUANT_TEST_PYTHON_ENVVARS}")

  unset(SEQUANT_TEST_PYTHON_ENVVARS)

  add_dependencies(check-sequant python-sequant)

endif (BUILD_TESTING)
