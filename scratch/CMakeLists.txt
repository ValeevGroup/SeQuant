add_executable(
  scratch EXCLUDE_FROM_ALL

  # ../examples/eval/btas_eval.cpp

  ../examples/eval/ta_eval.cpp
  ../examples/hf/hartree-fock.h
  ../examples/hf/hartree-fock.cc
)

target_link_libraries(scratch SeQuant ${TBB_LIBRARIES})

if (TARGET tiledarray)
    target_compile_definitions(scratch PRIVATE SEQUANT_HAS_TILEDARRAY)
    target_link_libraries(scratch Libint2::libint2_cxx tiledarray)
endif(TARGET tiledarray)

if (BTAS_SOURCE_DIR AND NOT TARGET tiledarray)
    target_include_directories(scratch PUBLIC ${BTAS_SOURCE_DIR} ${Boost_INCLUDE_DIRS})
endif(BTAS_SOURCE_DIR AND NOT TARGET tiledarray)
