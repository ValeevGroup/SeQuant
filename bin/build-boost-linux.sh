#! /bin/sh

export BOOST_VERSION=1_73_0

# Exit on error
set -ev

# download+unpack (but not build!) Boost unless previous install is cached ... must manually wipe cache on version bump or toolchain update
export INSTALL_DIR=${INSTALL_PREFIX}/boost
if [ ! -d "${INSTALL_DIR}" ]; then
    wget -q https://dl.bintray.com/boostorg/release/1.73.0/source/boost_${BOOST_VERSION}.tar.bz2
    tar -xjf boost_${BOOST_VERSION}.tar.bz2
    cd boost_${BOOST_VERSION}
    ./bootstrap.sh --prefix=${INSTALL_DIR}
    ./b2 --with-serialization
    ./b2 --with-serialization install
else
    echo "Boost already installed ..."
fi
