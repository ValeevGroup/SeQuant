SeQuant: installation guide
===========================

prerequisites:
  * mandatory:
    * CMake 3.15 or later
    * a C++17 compiler
    * [Boost](https://www.boost.org/), version 1.67 or higher (N.B. critical bugs make the following versions unusable: 1.70, 1.77, 1.78)
    * [Range-V3](https://github.com/ericniebler/range-v3.git), tag d800a032132512a54c291ce55a2a43e0460591c7, *if not found, SeQuant will download and build Range-V3*
  * optional:
    * unless CMake variable `BUILD_TESTING` is set to `FALSE`
      * [Eigen](http://eigen.tuxfamily.org/), version 3
    * for building coupled-cluster evaluation tests:
      * [TiledArray](https://github.com/ValeevGroup/tiledarray.git), tag 36e2ad205c21c339434dd0ef8f4f1467e7e26037

for the impatient (from the top of the SeQuant source directory):
  * `cmake -B build -S .`
  * `cmake --build build --target check-sequant`

useful CMake variables:
  * `BUILD_TESTING` --- enables unit tests targets, e.g. `check-sequant` [default=ON]
  * `CMAKE_CXX_COMPILER` --- specifies the C++ compiler to use
  * `CMAKE_PREFIX_PATH` --- this semicolon-separated list specifies search paths for dependencies (Boost, Range-V3, etc.)
