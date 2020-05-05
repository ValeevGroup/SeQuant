SeQuant: installation guide
===========================

prerequisites:
  * CMake 3.15
  * a C++17 compiler
  * a recent (1.67 or later) Boost library (N.B. Boost.Container is broken in 1.70)
  * the HEAD version of Range-V3 library (see https://github.com/ericniebler/range-v3.git)
  * optional:
    * TiledArray library + its prerequisites
    * Libint library + its prerequisites:

for the impatient:
  * `cmake . -DRANGEV3_DIR=...`
    * `RANGEV3_DIR` = path to the top of Range-V3 source tree
  * `make check`

CMake variables:
  * mandatory:
    * `RANGEV3_DIR` = path to the top of Range-V3 source tree
  * optional:
    * `CMAKE_PREFIX_PATH` = add TiledArray + Libint install directory to this

See `config.sh` for an example CMake config generation (example is for a Linux system).
