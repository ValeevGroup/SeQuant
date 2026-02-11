include(FindOrFetchCatch2)
include(Catch)

##########################
# OBJECT library: Symbolic layer tests (SeQuant-symb)
##########################
set(symb_test_sources
    "test_abstract_tensor.cpp"
    "test_asy_cost.cpp"
    "test_binary_node.cpp"
    "test_bliss.cpp"
    "test_canonicalize.cpp"
    "test_expr.cpp"
    "test_index.cpp"
    "test_iterator.cpp"
    "test_latex.cpp"
    "test_macros.cpp"
    "test_math.cpp"
    "test_meta.cpp"
    "test_op.cpp"
    "test_parse.cpp"
    "test_runtime.cpp"
    "test_space.cpp"
    "test_tensor.cpp"
    "test_tensor_network.cpp"
    "test_utilities.cpp"
    "test_wick.cpp"
)

add_library(unit_tests-sequant-symb-obj OBJECT ${symb_test_sources})
target_link_libraries(unit_tests-sequant-symb-obj
    PUBLIC SeQuant::symb SeQuant::bliss
    PRIVATE Catch2::Catch2 dtl::dtl)
target_compile_definitions(unit_tests-sequant-symb-obj PRIVATE
    SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

##########################
# OBJECT library: Eval framework tests (SeQuant-eval)
##########################
set(eval_test_sources
    "test_eval_expr.cpp"
    "test_eval_node.cpp"
)
add_library(unit_tests-sequant-eval-obj OBJECT ${eval_test_sources})
target_link_libraries(unit_tests-sequant-eval-obj
    PUBLIC SeQuant::eval
    PRIVATE Catch2::Catch2 dtl::dtl)
target_compile_definitions(unit_tests-sequant-eval-obj PRIVATE
    SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

##########################
# OBJECT library: Optimize tests (SeQuant-optimize)
##########################
set(optimize_test_sources
    "test_fusion.cpp"
    "test_optimize.cpp"
)
add_library(unit_tests-sequant-optimize-obj OBJECT ${optimize_test_sources})
target_link_libraries(unit_tests-sequant-optimize-obj
    PUBLIC SeQuant::optimize SeQuant::bliss
    PRIVATE Catch2::Catch2 dtl::dtl)
target_compile_definitions(unit_tests-sequant-optimize-obj PRIVATE
    SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

##########################
# OBJECT library: Export tests (SeQuant-export)
##########################
set(export_test_sources
    "test_export.cpp"
    "test_export.hpp"
    "test_export_python.cpp"
)
add_library(unit_tests-sequant-export-obj OBJECT ${export_test_sources})
target_link_libraries(unit_tests-sequant-export-obj
    PUBLIC SeQuant::export SeQuant::optimize
    PRIVATE Catch2::Catch2 dtl::dtl Eigen3::Eigen)
target_compile_definitions(unit_tests-sequant-export-obj PRIVATE
    SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

##########################
# OBJECT library: MBPT tests
##########################
set(mbpt_test_sources
    "test_biorthogonalization.cpp"
    "test_mbpt.cpp"
    "test_mbpt_cc.cpp"
    "test_spin.cpp"
)

add_library(unit_tests-sequant-mbpt-obj OBJECT ${mbpt_test_sources})
target_link_libraries(unit_tests-sequant-mbpt-obj
    PUBLIC SeQuant::mbpt
    PRIVATE Catch2::Catch2 dtl::dtl)
target_compile_definitions(unit_tests-sequant-mbpt-obj PRIVATE
    SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

##########################
# OBJECT library: TiledArray backend tests (conditional)
##########################
if (SEQUANT_HAS_TILEDARRAY)
    set(eval_ta_test_sources
        "test_eval_ta.cpp"
        "test_cache_manager.cpp"
    )
    add_library(unit_tests-sequant-eval-ta-obj OBJECT ${eval_ta_test_sources})
    target_link_libraries(unit_tests-sequant-eval-ta-obj
        PUBLIC SeQuant::eval::ta SeQuant::mbpt
        PRIVATE Catch2::Catch2 dtl::dtl)
    target_compile_definitions(unit_tests-sequant-eval-ta-obj PRIVATE
        SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    set_source_files_properties(${eval_ta_test_sources}
        PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
endif()

##########################
# OBJECT library: BTAS backend tests (conditional)
##########################
if (SEQUANT_HAS_BTAS)
    set(eval_btas_test_sources "test_eval_btas.cpp")
    # Include cache_manager test if TA is not available (avoid duplicate)
    if (NOT SEQUANT_HAS_TILEDARRAY)
        list(APPEND eval_btas_test_sources "test_cache_manager.cpp")
    endif()
    add_library(unit_tests-sequant-eval-btas-obj OBJECT ${eval_btas_test_sources})
    target_link_libraries(unit_tests-sequant-eval-btas-obj
        PUBLIC SeQuant::eval::btas SeQuant::mbpt
        PRIVATE Catch2::Catch2 dtl::dtl)
    target_compile_definitions(unit_tests-sequant-eval-btas-obj PRIVATE
        SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    set_source_files_properties(${eval_btas_test_sources}
        PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
endif()

##########################
# OBJECT library: TAPP backend tests (conditional)
##########################
if (SEQUANT_HAS_TAPP)
    set(eval_tapp_test_sources "test_eval_tapp.cpp")
    add_library(unit_tests-sequant-eval-tapp-obj OBJECT ${eval_tapp_test_sources})
    target_link_libraries(unit_tests-sequant-eval-tapp-obj
        PUBLIC SeQuant::eval::tapp SeQuant::mbpt
        PRIVATE Catch2::Catch2 dtl::dtl)
    target_compile_definitions(unit_tests-sequant-eval-tapp-obj PRIVATE
        SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    set_source_files_properties(${eval_tapp_test_sources}
        PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
endif()

##########################
# Single executable combining all OBJECT libraries
# NB linking OBJECT libraries (since CMake 3.21) both includes their object
# files AND propagates their PUBLIC/INTERFACE usage requirements, so the
# SeQuant module link deps flow through automatically.
##########################
add_executable(unit_tests-sequant ${BUILD_BY_DEFAULT}
    "test_main.cpp"
)
target_link_libraries(unit_tests-sequant PRIVATE
    unit_tests-sequant-symb-obj
    unit_tests-sequant-eval-obj
    unit_tests-sequant-optimize-obj
    unit_tests-sequant-export-obj
    unit_tests-sequant-mbpt-obj
    Catch2::Catch2
    dtl::dtl)
if (SEQUANT_HAS_TILEDARRAY)
    target_link_libraries(unit_tests-sequant PRIVATE unit_tests-sequant-eval-ta-obj)
endif()
if (SEQUANT_HAS_BTAS)
    target_link_libraries(unit_tests-sequant PRIVATE unit_tests-sequant-eval-btas-obj)
endif()
if (SEQUANT_HAS_TAPP)
    target_link_libraries(unit_tests-sequant PRIVATE unit_tests-sequant-eval-tapp-obj)
endif()

target_compile_definitions(unit_tests-sequant PRIVATE
    SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

# Skip unity build for test_main.cpp
set_source_files_properties("test_main.cpp"
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)

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
        target_compile_definitions(unit_tests-sequant-export-obj PRIVATE
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
        target_compile_definitions(unit_tests-sequant-export-obj PRIVATE
            SEQUANT_HAS_TORCH_FOR_VALIDATION=1)
    else()
        message(STATUS "PyTorch not found - PythonEinsumGenerator PyTorch validation tests will be disabled")
    endif()

    if (_numpy_status EQUAL 0 OR _torch_status EQUAL 0)
        target_compile_definitions(unit_tests-sequant-export-obj PRIVATE
            SEQUANT_UNITTESTS_PYTHON_EXECUTABLE="${Python_EXECUTABLE}")
    endif()
else()
    message(STATUS "Python not found - PythonEinsumGenerator validation tests will be disabled")
endif()

if (SEQUANT_INTERNAL_SKIP_LONG_TESTS)
    target_compile_definitions(unit_tests-sequant-symb-obj PRIVATE SEQUANT_SKIP_LONG_TESTS=1)
    target_compile_definitions(unit_tests-sequant-mbpt-obj PRIVATE SEQUANT_SKIP_LONG_TESTS=1)
endif()

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
