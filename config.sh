#!/usr/bin/env bash


rm -rf ./build

# INTEL_DIR=/opt/intel
# source /opt/intel/tbb/bin/tbbvars.sh "intel64"

cmake -B build                            \
      -D RANGEV3_DIR=/opt/range-v3        \
      -D BTAS_INSTALL_DIR=/opt/BTAS       \
      -D CMAKE_PREFIX_PATH=/opt/libint    \
      -D MKL_THREADING=TBB                \
      -D CMAKE_EXPORT_COMPILE_COMMANDS=ON \
      -D CMAKE_BUILD_TYPE=Debug           \
      --verbose                           \
      -S .
      # -G"Ninja"                          \
      # -D TiledArray_DIR=/opt/tiledarray  \
      # -D BOOST_ROOT=/usr                 \
