
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

if (NOT Boost_FOUND)
  find_package(Boost REQUIRED)
endif()

if (BTAS_INCLUDE_DIRS)
  set (BTAS_FOUND 1)

else (BTAS_INCLUDE_DIRS)
  set (BTAS_FOUND 0)

  if(BTAS_INSTALL_DIR)
    set(_BTAS_INSTALL_DIR ${BTAS_INSTALL_DIR})
  endif(BTAS_INSTALL_DIR)

  set(BTAS_INC_SEARCH_DIR ${_BTAS_INSTALL_DIR})

  find_path(BTAS_INCLUDE_DIR NAMES btas/btas.h
    HINTS
    ${BTAS_INC_SEARCH_DIR})

  mark_as_advanced(BTAS_INCLUDE_DIR)

  # CODA
  if(BTAS_INCLUDE_DIR)
  
    # validate version, etc. by compiling tests
    cmake_push_check_state()
  
    list(APPEND CMAKE_REQUIRED_INCLUDES ${Boost_INCLUDE_DIRS} ${BTAS_INCLUDE_DIR})
    list(APPEND CMAKE_REQUIRED_DEFINITIONS ${CMAKE_CXX14_STANDARD_COMPILE_OPTION})
    
    # sanity check
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <btas/btas.h>
      #include <array>
      #include <complex>
      int main(int argc, char** argv) {
        btas::Tensor<std::array<std::complex<double>,3>> T1(2,3,4);
        T1.fill({{{1.0,2.0}, {2.0,1.0}, {2.0,3.0} }});
        return 0;
      }
      "  BTAS_COMPILES)
    if (NOT BTAS_COMPILES)
      message(FATAL_ERROR "Could not compile BTAS test program.\nSee CMakeFiles/CMakeError.log for details")
    endif()
    
    cmake_pop_check_state()
  
    set(BTAS_FOUND 1)
    set(BTAS_INCLUDE_DIRS ${BTAS_INCLUDE_DIR})
    mark_as_advanced(BTAS_INCLUDE_DIRS)
    message(STATUS "Found BTAS library: include ${BTAS_INCLUDE_DIRS}")
  endif(BTAS_INCLUDE_DIR)


endif(BTAS_INCLUDE_DIRS)

# create BTAS interface library
if (BTAS_FOUND)
  # look for CBLAS + LAPACKE
  find_package(CBLAS)
  find_package(LAPACKE)
  add_library(BTAS INTERFACE IMPORTED)
  if (${CBLAS_FOUND} AND ${LAPACKE_FOUND})
    set(BTAS_DEFINITIONS "_CBLAS_HEADER=\"${CBLAS_INCLUDE_FILE}\";_LAPACKE_HEADER=\"${LAPACKE_INCLUDE_FILE}\";BTAS_HAS_CBLAS=1")
    if (MKL_FOUND)
      set(BTAS_DEFINITIONS "_HAS_INTEL_MKL=1;${BTAS_DEFINITIONS}")
    endif(MKL_FOUND)
    message(STATUS "BTAS_DEFINITIONS=${BTAS_DEFINITIONS}")
    set_target_properties(BTAS PROPERTIES
            INTERFACE_LINK_LIBRARIES "${LAPACKE_LIBRARIES};${CBLAS_LIBRARIES}"
            INTERFACE_COMPILE_DEFINITIONS "${BTAS_DEFINITIONS}"
            INTERFACE_INCLUDE_DIRECTORIES "${BTAS_INCLUDE_DIRS};${LAPACKE_INCLUDE_DIR};${CBLAS_INCLUDE_DIR}"
            )
  else()
    set_target_properties(BTAS PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES ${BTAS_INCLUDE_DIRS}
            )
  endif()
endif(BTAS_FOUND)
