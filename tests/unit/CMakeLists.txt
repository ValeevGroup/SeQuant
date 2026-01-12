include(FindOrFetchCatch2)
include(Catch)

add_executable(unit_tests-sequant ${BUILD_BY_DEFAULT}
    "test_abstract_tensor.cpp"
    "test_asy_cost.cpp"
    "test_binary_node.cpp"
    "test_biorthogonalization.cpp"
    "test_bliss.cpp"
    "test_canonicalize.cpp"
    "test_eval_expr.cpp"
    "test_eval_node.cpp"
    "test_export.cpp"
    "test_export.hpp"
    "test_export_python.cpp"
    "test_expr.cpp"
    "test_expr.cpp"
    "test_fusion.cpp"
    "test_index.cpp"
    "test_iterator.cpp"
    "test_latex.cpp"
    "test_main.cpp"
    "test_macros.cpp"
    "test_math.cpp"
    "test_mbpt.cpp"
    "test_mbpt_cc.cpp"
    "test_meta.cpp"
    "test_op.cpp"
    "test_optimize.cpp"
    "test_parse.cpp"
    "test_runtime.cpp"
    "test_space.cpp"
    "test_spin.cpp"
    "test_string.cpp"
    "test_tensor.cpp"
    "test_tensor_network.cpp"
    "test_utilities.cpp"
    "test_wick.cpp"
)

target_compile_definitions(unit_tests-sequant PRIVATE SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(unit_tests-sequant PRIVATE Boost::boost Eigen3::Eigen)

# Optional: Find Python for PythonEinsumGenerator validation tests
find_package(Python 3 COMPONENTS Interpreter QUIET)
if(Python_Interpreter_FOUND)
    message(STATUS "Found Python: ${Python_EXECUTABLE} (version ${Python_VERSION}) for export validation tests")

    # Check for NumPy
    execute_process(
        COMMAND "${Python_EXECUTABLE}" "-c" "import numpy; print(numpy.__version__)"
        RESULT_VARIABLE _numpy_status
        OUTPUT_VARIABLE _numpy_version
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(_numpy_status EQUAL 0)
        message(STATUS "NumPy ${_numpy_version} found for export validation tests")
        target_compile_definitions(unit_tests-sequant PRIVATE
            SEQUANT_HAS_NUMPY_FOR_VALIDATION=1)
    else()
        message(STATUS "NumPy not found - PythonEinsumGenerator NumPy validation tests will be disabled")
    endif()

    # Check for PyTorch
    execute_process(
        COMMAND "${Python_EXECUTABLE}" "-c" "import torch; print(torch.__version__)"
        RESULT_VARIABLE _torch_status
        OUTPUT_VARIABLE _torch_version
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(_torch_status EQUAL 0)
        message(STATUS "PyTorch ${_torch_version} found for export validation tests")
        target_compile_definitions(unit_tests-sequant PRIVATE SEQUANT_HAS_TORCH_FOR_VALIDATION=1)
    else()
        message(STATUS "PyTorch not found - PythonEinsumGenerator PyTorch validation tests will be disabled")
    endif()

    if (_numpy_status EQUAL 0 OR _torch_status EQUAL 0)
        target_compile_definitions(unit_tests-sequant PRIVATE
                SEQUANT_UNITTESTS_PYTHON_EXECUTABLE="${Python_EXECUTABLE}")
    endif()
else()
    message(STATUS "Python not found - PythonEinsumGenerator validation tests will be disabled")
endif()

if (SEQUANT_SKIP_LONG_TESTS)
    target_compile_definitions(unit_tests-sequant PRIVATE SEQUANT_SKIP_LONG_TESTS=1)
endif()

if (TARGET tiledarray)
    target_sources(unit_tests-sequant
        PRIVATE
            "test_cache_manager.cpp"
            "test_eval_btas.cpp"
            "test_eval_ta.cpp"
    )
    set_source_files_properties(
        "test_eval_btas.cpp"
        "test_eval_ta.cpp"
        "test_main.cpp"
        PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
    )
    target_link_libraries(unit_tests-sequant PRIVATE tiledarray)
    target_compile_definitions(unit_tests-sequant PRIVATE SEQUANT_HAS_TILEDARRAY)
endif (TARGET tiledarray)

target_link_libraries(unit_tests-sequant PRIVATE SeQuant::SeQuant Catch2::Catch2 dtl::dtl)

target_set_warning_flags(unit_tests-sequant)

if (SEQUANT_TESTS)
    catch_discover_tests(
        unit_tests-sequant
        TEST_PREFIX "sequant/unit/"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        PROPERTIES FIXTURES_REQUIRED "SEQUANT_BUILD_ALL"
    )
else()
    add_test(
        NAME "sequant/unit"
        COMMAND unit_tests-sequant
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )

    set(TEST_NAMES "sequant/unit")
    build_test_as_needed(unit_tests-sequant "sequant/unit" TEST_NAMES)
endif()
