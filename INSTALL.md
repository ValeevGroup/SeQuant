# Installation Guide


### TL;DR
While in the SeQuant source directory:
* `cmake -B build -S .`
* `cmake --build build --target check-sequant`

### Prerequisites
  * mandatory:
    * CMake 3.27 or later
    * a C++20 compiler
    * [Boost](https://www.boost.org/), version 1.81 or higher. *SeQuant can download and build Boost if configured with `Boost_FETCH_IF_MISSING=ON`, but the use of Boost provided by the system package manager is recommended.* The following non-header-only Boost libraries are required, hence Boost must be configured/built:
      - [Boost.Regex](https://www.boost.org/doc/libs/master/libs/regex/doc/html/index.html)
      - [Boost.Locale](https://www.boost.org/doc/libs/master/libs/locale/doc/html/index.html)
    * [Range-V3](https://github.com/ericniebler/range-v3.git), tag 0.12.0, *if not found, SeQuant will download and build Range-V3*
    * [libPerm](https://github.com/Krzmbrzl/libPerm), tag cada3e185549896203cf4d0c7f26ea22c7de428f, *if not found, SeQuant will download and build libPerm*
  * optional:
    * for building coupled-cluster evaluation tests:
      * [TiledArray](https://github.com/ValeevGroup/tiledarray.git), tag 5477174cd8e857f76d29cf6c0dcc80277393a9eb
    * for building `stcc*` test programs
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
  * `SEQUANT_PYTHON` --- enables building of Python bindings
  * `Boost_FETCH_IF_MISSING` --- if set to `ON`, SeQuant will download and build Boost if it is not found by `find_package(Boost ...)`; this is not recommended. [default=OFF]
  * `SEQUANT_CONTEXT_MANIPULATION_THREADSAFE` --- set to `OFF` if you will initialize the default context and change its state from a single thread only (most users will want to do this). [default=ON]

### Build

To build, test, and install SeQuant, run the following commands:
* `cmake --build build`
* (optional) `cmake --build build --target check-sequant`
* `cmake --build build --target install`
