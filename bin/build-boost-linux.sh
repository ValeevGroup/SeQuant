#! /bin/sh

export BOOST_VERSION=1_73_0

# Exit on error
set -ev

if [ "$CXX" = "g++" ]; then
    export CXX=/usr/bin/g++-$GCC_VERSION
    export CXXFLAGS="-mno-avx"
    export BOOST_TOOLSET=gcc
else
    export CXX=/usr/bin/clang++-$CLANG_VERSION
    export CXXFLAGS="-mno-avx -stdlib=libc++"
    export BOOST_TOOLSET=clang
fi

if [ "X$BUILD_TYPE" = "XDebug" ]; then
    export BOOST_VARIANT="debug"
else
    export BOOST_VARIANT="release"
fi

# download+unpack (but not build!) Boost unless previous install is cached ... must manually wipe cache on version bump or toolchain update
export INSTALL_DIR=${INSTALL_PREFIX}/boost
if [ ! -d "${INSTALL_DIR}" ]; then
    wget -q https://dl.bintray.com/boostorg/release/1.73.0/source/boost_${BOOST_VERSION}.tar.bz2
    tar -xjf boost_${BOOST_VERSION}.tar.bz2
    cd boost_${BOOST_VERSION}
    cat > user-config.jam << END
using ${BOOST_TOOLSET} : : ${CXX} :
      <cxxflags>"${CXXFLAGS}"
      <linkflags>"${CXXFLAGS}" ;
END
    ./bootstrap.sh --prefix=${INSTALL_DIR} --with-libraries=serialization
    ./b2 --user-config=`pwd`/user-config.jam toolset=${BOOST_TOOLSET} link=static variant=${BOOST_VARIANT}
    ./b2 install
else
    echo "Boost already installed ..."
fi
