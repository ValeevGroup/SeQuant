# Installation Guide


### TL;DR
While in the SeQuant source directory:
* `cmake -B build -S .`
* `cmake --build build --target check-sequant`

### Prerequisites
  * mandatory:
    * CMake 3.15 or later
    * a C++17 compiler
    * [Boost](https://www.boost.org/), version 1.81 or higher. *SeQuant can download and build Boost if configured with `Boost_FETCH_IF_MISSING=ON`, but this is not recommended.* The following non-header-only Boost libraries are required, hence Boost must be configured/built:
      - [Boost.Regex](https://www.boost.org/doc/libs/master/libs/regex/doc/html/index.html)
      - [Boost.Locale](https://www.boost.org/doc/libs/master/libs/locale/doc/html/index.html)
    * [Range-V3](https://github.com/ericniebler/range-v3.git), tag 0.12.0, *if not found, SeQuant will download and build Range-V3*
  * optional:
    * for building coupled-cluster evaluation tests:
      * [TiledArray](https://github.com/ValeevGroup/tiledarray.git), tag 3f6629db047417e814b75ad5069b7f4ce26428e7
    * for building `stcc*` example programs
        * [Eigen](http://eigen.tuxfamily.org/), version 3

### Configure

From the SeQuant source directory run the following command to configure the build harness:
* `cmake -B build -S .`

#### Useful CMake variables
  * [`BUILD_TESTING`](https://cmake.org/cmake/help/latest/module/CTest.html) --- enables unit tests targets, e.g. `check-sequant` [default=ON]
  * [`CMAKE_CXX_COMPILER`](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER.html#variable:CMAKE_%3CLANG%3E_COMPILER) --- specifies the C++ compiler to use
  * [`CMAKE_PREFIX_PATH`](https://cmake.org/cmake/help/latest/variable/CMAKE_PREFIX_PATH.html) --- this semicolon-separated list specifies search paths for dependencies (Boost, Range-V3, etc.)
  * [`CMAKE_INSTALL_PREFIX`](https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html) --- the installation path for SeQuant
  * `SEQUANT_MIMALLOC` --- use [mimalloc](https://github.com/microsoft/mimalloc) for fast memory allocation
  * `SEQUANT_EVAL_TRACE` --- enables tracing of expression interpretation; especially useful in combination with TiledArray's memory tracing mechanism (configure TiledArray with `TA_TENSOR_MEM_PROFILE=ON` to enable that)
  * `SEQUANT_PYTHON` --- enables building of Pythin bindings
  * `Boost_FETCH_IF_MISSING` --- if set to `ON`, SeQuant will download and build Boost if it is not found by `find_package(Boost ...)`; this is not recommended. [default=OFF]

### Build

To build, test, and install SeQuant, run the following commands:
* `cmake --build build`
* (optional) `cmake --build build --target check-sequant`
* `cmake --build build --target install`
