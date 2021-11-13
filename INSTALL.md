SeQuant: installation guide
===========================

prerequisites:
  * mandatory:
    * CMake 3.15 or later
    * a C++17 compiler
    * a recent (1.67 or later) Boost library (N.B. Boost.Container is broken in 1.70)
    * a recent Range-V3 library (see https://github.com/ericniebler/range-v3.git), *if not found, SeQuant will download and build Range-V3*
    * Eigen, version 3.
  * optional:
    * for building coupled-cluster evaluation tests:
      * TiledArray library + its prerequisites

for the impatient:
  * `cmake -B build -S .`
  * `cmake --build build --target check`

CMake variables:
  * optional:
    * `CMAKE_PREFIX_PATH` = add Boost, Range-V3, TiledArray, Libint install prefix paths to this (semicolon-separated)
