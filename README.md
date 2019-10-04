SeQuant: second quantization toolkit in C++
===========================================

prerequisites:
  * CMake 3.9
  * a C++17 compiler
  * a recent (1.67 or later) Boost library (N.B. Boost.Container is broken in 1.70)
  * the HEAD version of Range-V3 library (see https://github.com/ericniebler/range-v3.git)
  * optional (if `test_contract` is needed):
    * BTAS library + its prerequisites:
      * CBLAS + LAPACKE interfaces to BLAS + LAPACK libraries, as provided by e.g. MKL
    * Libint library + its prerequisites:
      * Eigen

for the impatient:
  * `cmake . -DRANGEV3_DIR=...`
    * `RANGEV3_DIR` = path to the top of Range-V3 source tree
  * `make check`

CMake variables:
  * mandatory
    * `RANGEV3_DIR` = path to the top of Range-V3 source tree
  * optional (if `test_contract` is needed):
    * `BTAS_INSTALL_DIR` = path to the top of the BTAS source directory
    * `BOOST_ROOT` = path to the top of the Boost install tree
    * `CMAKE_PREFIX_PATH` = add Libint install directory to this

