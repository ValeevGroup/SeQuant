#! /bin/sh

export BOOST_VERSION=1_72_0

# Exit on error
set -ev

# download+unpack (but not build!) Boost unless previous install is cached ... must manually wipe cache on version bump or toolchain update
export INSTALL_DIR=${INSTALL_PREFIX}/boost_${BOOST_VERSION}
if [ ! -d "${INSTALL_DIR}" ]; then
    cd ${INSTALL_PREFIX}
    wget -q https://dl.bintray.com/boostorg/release/1.72.0/source/boost_${BOOST_VERSION}.tar.bz2
    tar -xjf boost_${BOOST_VERSION}.tar.bz2
    # make shortcut
    rm -f boost
    ln -s boost_${BOOST_VERSION} boost
else
    echo "Boost already installed ..."
fi
