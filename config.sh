#!/usr/bin/env bash

# set up MKL variables

rm -rf ./build

INTEL_DIR=/opt/intel
source ${INTEL_DIR}/compilers_and_libraries/linux/mkl/bin/mklvars.sh "intel64"
source /opt/intel/compilers_and_libraries/linux/tbb/bin/tbbvars.sh "intel64"
cmake -Bbuild \
      -DRANGEV3_DIR=/opt/range-v3 \
      -DBTAS_INSTALL_DIR=/opt/BTAS \
      -DBOOST_ROOT=/usr \
      -DCMAKE_PREFIX_PATH=/opt/libint \
      -DMKL_THREADING=TBB

