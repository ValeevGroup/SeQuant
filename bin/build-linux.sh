#! /bin/bash

# wipe out cache if commit message contains [travis delete-cache]
# see https://github.com/travis-ci/travis-ci/issues/2245
ls ${INSTALL_PREFIX}
if [[ "${TRAVIS_COMMIT_MESSAGE}" = *"[travis delete-cache]"* ]]; then
  echo "Clearing out cache due to [travis delete-cache] found in the commit message"
  rm -rf ${INSTALL_PREFIX}/*
  ccache -C
  ccache -c
fi

${TRAVIS_BUILD_DIR}/bin/build-rangev3-$TRAVIS_OS_NAME.sh
${TRAVIS_BUILD_DIR}/bin/build-btas-$TRAVIS_OS_NAME.sh
${TRAVIS_BUILD_DIR}/bin/build-eigen3-$TRAVIS_OS_NAME.sh
${TRAVIS_BUILD_DIR}/bin/build-libint-$TRAVIS_OS_NAME.sh

# Exit on error
set -ev

#
# Environment variables
#
export CXX_FLAGS="-mno-avx -ftemplate-depth=1024"
if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
    # ggc-min params to try to reduce peak memory consumption ... typical ICE under Travis is due to insufficient memory
    export EXTRAFLAGS="--param ggc-min-expand=20 --param ggc-min-heapsize=2048000"
else
    export CC=/usr/bin/clang-$CLANG_VERSION
    export CXX=/usr/bin/clang++-$CLANG_VERSION
    export EXTRAFLAGS="-Wno-unused-command-line-argument"
fi

echo $($CC --version)
echo $($CXX --version)

# where to install
export INSTALL_DIR=${INSTALL_PREFIX}/SeQuant

# make build dir
cd ${BUILD_PREFIX}
mkdir -p SeQuant
cd SeQuant

# configure CodeCov only for g++-8 debug build
if [ "$BUILD_TYPE" = "Debug" ] && [ "$GCC_VERSION" = 8 ] && [ "$CXX" = "g++" ]; then
    export CODECOVCXXFLAGS="--coverage -O0"
fi

cmake ${TRAVIS_BUILD_DIR} \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DRANGEV3_DIR="${INSTALL_PREFIX}/range-v3" \
    -DBTAS_INSTALL_DIR="${INSTALL_PREFIX}/BTAS" \
    -DCMAKE_PREFIX_PATH="${INSTALL_PREFIX}/libint;${INSTALL_PREFIX}/eigen3" \
    -DCMAKE_CXX_FLAGS="${CXX_FLAGS} ${EXTRAFLAGS} ${CODECOVCXXFLAGS}"

### test within build tree
make -j2 check VERBOSE=1

# print ccache stats
ccache -s

