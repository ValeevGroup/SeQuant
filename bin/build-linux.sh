#! /bin/sh

${TRAVIS_BUILD_DIR}/bin/build-rangev3-$TRAVIS_OS_NAME.sh

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
export INSTALL_DIR=${INSTALL_PREFIX}/SeQuant2

# make build dir
cd ${BUILD_PREFIX}
mkdir -p SeQuant2
cd SeQuant2

# configure CodeCov only for gcc-7 debug build
if [ "$BUILD_TYPE" = "Debug" ] && [ "$GCC_VERSION" = 7 ]; then
    export CODECOVCXXFLAGS="--coverage -O0"
fi

cmake ${TRAVIS_BUILD_DIR} \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DRANGEV3_DIR="${INSTALL_DIR}/range-v3" \
    -DCMAKE_CXX_FLAGS="${CXX_FLAGS} ${EXTRAFLAGS} ${CODECOVCXXFLAGS}"

### test within build tree
make -j2 check VERBOSE=1
make srcc VERBOSE=1
examples/srcc 3

# print ccache stats
ccache -s

