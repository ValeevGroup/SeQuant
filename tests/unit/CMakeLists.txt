include(FindOrFetchCatch2)
include(Catch)

add_executable(unit_tests-sequant ${BUILD_BY_DEFAULT}
    "test_asy_cost.cpp"
    "test_binary_node.cpp"
	"test_biorthogonalization.cpp"
    "test_bliss.cpp"
    "test_canonicalize.cpp"
    "test_eval_expr.cpp"
    "test_eval_node.cpp"
    "test_export.cpp"
    "test_expr.cpp"
    "test_expr.cpp"
    "test_fusion.cpp"
    "test_index.cpp"
    "test_iterator.cpp"
    "test_latex.cpp"
    "test_main.cpp"
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

target_link_libraries(unit_tests-sequant PRIVATE SeQuant Catch2 dtl::dtl)

if (SEQUANT_TESTS)
    catch_discover_tests(
        unit_tests-sequant
        TEST_PREFIX "sequant/unit/"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
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
