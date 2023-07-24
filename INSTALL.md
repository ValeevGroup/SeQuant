SeQuant: installation guide
===========================

prerequisites:
  * mandatory:
    * CMake 3.15 or later
    * a C++17 compiler
    * [Boost](https://www.boost.org/), version 1.67 or higher (N.B. critical bugs make the following versions unusable: 1.70, 1.77, 1.78); the following non-header-only Boost libraries are required:
      - [Boost.Regex](https://www.boost.org/doc/libs/master/libs/regex/doc/html/index.html)
      - [Boost.Locale](https://www.boost.org/doc/libs/master/libs/locale/doc/html/index.html)
    * [Range-V3](https://github.com/ericniebler/range-v3.git), tag d800a032132512a54c291ce55a2a43e0460591c7, *if not found, SeQuant will download and build Range-V3*
  * optional:
    * for building coupled-cluster evaluation tests:
      * [TiledArray](https://github.com/ValeevGroup/tiledarray.git), tag e52c733e9f1400e70425d5be171a8e7cd5b56876
    * for building `stcc*` example programs
        * [Eigen](http://eigen.tuxfamily.org/), version 3

for the impatient (from the top of the SeQuant source directory):
  * `cmake -B build -S .`
  * `cmake --build build --target check-sequant`

useful CMake variables:
  * `BUILD_TESTING` --- enables unit tests targets, e.g. `check-sequant` [default=ON]
  * `CMAKE_CXX_COMPILER` --- specifies the C++ compiler to use
  * `CMAKE_PREFIX_PATH` --- this semicolon-separated list specifies search paths for dependencies (Boost, Range-V3, etc.)
  * `SEQUANT_MIMALLOC` --- use [mimalloc](https://github.com/microsoft/mimalloc) for fast memory allocation
  * `SEQUANT_EVAL_TRACE` --- enables tracing of expression interpretation; especially useful in combination with TiledArray's memory tracing mechanism (configure TiledArray with `TA_TENSOR_MEM_PROFILE=ON` to enable that)
