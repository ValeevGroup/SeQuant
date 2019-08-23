#!/usr/bin/env bash
source /opt/intel/tbb/bin/tbbvars.sh "intel64"
cmake --build ../build --target srcc -v -j7 && ../build/srcc
