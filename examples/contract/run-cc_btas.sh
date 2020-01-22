#!/usr/bin/env bash

INTEL_DIR=/opt/intel

source ${INTEL_DIR}/mkl/bin/mklvars.sh "intel64"

source ${INTEL_DIR}/tbb/bin/tbbvars.sh "intel64"

export MKL_NUM_THREADS=7

build_dir=$(realpath $0)
build_dir=${build_dir%/examples*}
build_dir=${build_dir}/cmake-build-debug

[[ -d ${build_dir} ]] ||
    build_dir=${build_dir%/cmake-build-debug}/build

cmake_source_dir=${build_dir%/*}

target=$(basename $0)
target=${target%.*}
target=${target#run-}

executable=${build_dir}/${target}

geometry_input=${build_dir%/*}/examples/h2o.xyz

function run_exec(){
cmake --build ${build_dir} \
      --target ${target} \
      --verbose \
      -j 7 \
      && ${executable} ${geometry_input}
}


if [[ -d ${build_dir} ]] && [[ -x ${executable} ]]; then
    echo "The build_dir is ${build_dir}"
    echo "CMAKE source dir is ${cmake_source_dir}"
    echo "The target is ${target}"
    echo "The exec.  is ${executable}"
    echo "The geom.  is ${geometry_input}"
    run_exec
else
    echo
    echo "Executable ${executable} not found! Target ${target} needs to be built before."
    exit 1
fi
