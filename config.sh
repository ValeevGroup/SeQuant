#!/usr/bin/env bash

# set up MKL variables

rm -rf ./build

INTEL_DIR=/opt/intel
${INTEL_DIR}/bin/compilervars.sh -arch intel64 -platform linux
cmake -Bbuild \
      -DRANGEV3_DIR=/opt/range-v3 \
      -DBTAS_INSTALL_DIR=/opt/BTAS \
      -DBOOST_ROOT=/usr \
      -DCMAKE_PREFIX_PATH=/opt/libint \
      -DMKL_THREADING=TBB \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
      -G"Ninja"
