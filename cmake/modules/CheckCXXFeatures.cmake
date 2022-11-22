include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

# if TBB_FOUND is true will check for usable <execution> without TBB, then with TBB
macro(check_cxx_execution_header _prefix)

  ##############################################
  # compilation checks
  ##############################################
  set(_prereq_list "_STANDALONE")
  if (TARGET TBB::tbb)
    list(APPEND _prereq_list _WITH_TBB)
  endif ()

  foreach (_prereq ${_prereq_list})
    cmake_push_check_state()

    if (_prereq STREQUAL _WITH_TBB)
      list(APPEND CMAKE_REQUIRED_LIBRARIES TBB::tbb)
    endif ()

    CHECK_CXX_SOURCE_COMPILES(
        "
  #include <algorithm>
  #include <vector>
  #include <execution>
  int main(int argc, char** argv) {
    std::vector<int> v{0,1,2};
    std::for_each(std::execution::par_unseq, begin(v), end(v),
                  [](auto&& i) {i *= 2;});
    return 0;
  }
  " ${_prefix}_HAS_EXECUTION_HEADER${_prereq})

    cmake_pop_check_state()
    if (${_prefix}_HAS_EXECUTION_HEADER${_prereq})
      break()
    endif ()

  endforeach (_prereq)

endmacro(check_cxx_execution_header)
