add_executable(
  scratch EXCLUDE_FROM_ALL
  ../examples/eval/read_tensor.hpp
  ../examples/eval/read_tensor.cpp

  ../examples/eval/read_tensor_btas.hpp
  ../examples/eval/read_tensor_btas.cpp
  ../examples/eval/eval_btas.cpp

  # ../examples/eval/read_tensor_ta.cpp
  # ../examples/eval/read_tensor_ta.hpp
  # ../examples/eval/eval_ta.cpp
)

target_link_libraries(scratch SeQuant ${TBB_LIBRARIES})

if (TARGET tiledarray)
    target_compile_definitions(scratch PRIVATE SEQUANT_HAS_TILEDARRAY)
    target_link_libraries(scratch Libint2::libint2_cxx tiledarray)
endif(TARGET tiledarray)

if (BTAS_SOURCE_DIR AND NOT TARGET tiledarray)
    target_include_directories(scratch PUBLIC ${BTAS_SOURCE_DIR} ${Boost_INCLUDE_DIRS})
endif(BTAS_SOURCE_DIR AND NOT TARGET tiledarray)
