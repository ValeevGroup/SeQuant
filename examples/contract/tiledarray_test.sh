#!/usr/bin/env bash

INTEL_DIR=/opt/intel

source ${INTEL_DIR}/compilers_and_libraries/linux/mkl/bin/mklvars.sh "intel64"

source ${INTEL_DIR}/tbb/bin/tbbvars.sh "intel64"

MKL_NUM_THREADS=7

cmake --build ../../build \
      --target cc_tiledarray \
      --verbose \
      -j 7 \
      && ../../build/cc_tiledarray ./h2o.xyz
