SeQuant: installation guide
===========================

prerequisites:
  * mandatory:
    * CMake 3.15 or later
    * a C++17 compiler
    * Boost, version 1.67 or higher (N.B. critical bugs make the following versions unusable: 1.70, 1.77, 1.78)
    * [Range-V3](https://github.com/ericniebler/range-v3.git), tag 2e0591c57fce2aca6073ad6e4fdc50d841827864, *if not found, SeQuant will download and build Range-V3*
    * Eigen, version 3.
  * optional:
    * for building coupled-cluster evaluation tests:
      * [TiledArray](https://github.com/ValeevGroup/tiledarray.git), tag 5c768a7b121886dfe406c6dd6a38acaa8782ae6e

for the impatient:
  * `cmake -B build -S .`
  * `cmake --build build --target check`

CMake variables:
  * optional:
    * `CMAKE_PREFIX_PATH` = add Boost, Range-V3, TiledArray, Libint install prefix paths to this (semicolon-separated)
