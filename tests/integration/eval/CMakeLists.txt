option(SEQUANT_EVAL_TESTS "Enable building of evaluation tests (if true, will look for and/or build TiledArray)" OFF)
add_feature_info(EVAL_TESTS SEQUANT_EVAL_TESTS "Build evaluation tests (if true, will look for and/or build TiledArray)")

# uccf12 example moved to MPQC

add_library(eval_shared STATIC EXCLUDE_FROM_ALL
        "eval_utils.hpp"
        "options.hpp"
        "options.cpp"
        "data_info.hpp"
        "data_info.cpp"
        "calc_info.hpp"
        "calc_info.cpp"
        "scf.hpp"
        "scf.cpp"
)
target_link_libraries(eval_shared PUBLIC SeQuant)
target_include_directories(eval_shared PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

if (TARGET tiledarray)
    set(example5 eval_ta)
    add_executable(eval_ta ${BUILD_BY_DEFAULT}
            "ta/data_world_ta.hpp"
            "ta/scf_ta.hpp"
            "ta/main.cpp"
    )
    target_link_libraries(eval_ta PRIVATE eval_shared tiledarray)
    target_compile_definitions(eval_ta PRIVATE SEQUANT_HAS_TILEDARRAY)

    set(test_name "sequant/integration/eval_ta")
    add_test(
            NAME "${test_name}"
            COMMAND eval_ta
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )

    build_test_as_needed(eval_ta "${test_name}" test_name)

    add_executable(eval_btas ${BUILD_BY_DEFAULT}
            "btas/data_world_btas.hpp"
            "btas/scf_btas.hpp"
            "btas/main.cpp"
    )
    target_include_directories(eval_btas PRIVATE ${BTAS_SOURCE_DIR})
    target_link_libraries(eval_btas PRIVATE eval_shared)

    set(test_name "sequant/integration/eval_btas")
    add_test(
            NAME "${test_name}"
            COMMAND eval_btas
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )

    build_test_as_needed(eval_btas "${test_name}" test_name)
endif (TARGET tiledarray)
