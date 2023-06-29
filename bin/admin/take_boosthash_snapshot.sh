#!/bin/sh

set -e

VERSION=1.81.0
VERSION_="$(echo ${VERSION} | sed 's/\./_/g')"
TARBALL=boost_${VERSION_}.tar.bz2

#### download
if ! test -f ${TARBALL}; then
  BOOST_URL=https://boostorg.jfrog.io/artifactory/main/release/${VERSION}/source/${TARBALL}
  curl -o ${TARBALL} -L ${BOOST_URL} 
fi

#### unpack
if ! test -d boost_${VERSION_}; then
  tar -xjf ${TARBALL}
fi

#### build bcp
cd boost_${VERSION_}
./bootstrap.sh && ./b2 tools/bcp

#### extract boost/preprocessor.hpp and dependencies
mkdir tmp
bin.v2/tools/bcp/*darwin*/release/link-static/bcp --namespace=sequant_boost --unix-lines boost/container_hash/hash.hpp tmp

#### specialize to SeQuant by
# 1. prepending all includes with SeQuant/external/boost
find tmp/boost -type f -print0 | xargs -0 sed -i "" 's/ <boost\// <SeQuant\/external\/boost\//g'
# 2. renaming all header guards
find tmp/boost -type f -print0 | xargs -0 sed -i "" 's/ BOOST_\(.*\)_INCLUDED/ SEQUANT_BOOST_\1_INCLUDED/g'
find tmp/boost -type f -print0 | xargs -0 sed -i "" 's/ BOOST_\(.*\)_HPP/ SEQUANT_BOOST_\1_HPP/g'
# 3. renaming other macros
find tmp/boost -type f -print0 | xargs -0 sed -i "" 's/BOOST_HAS_TRIVIAL_/SEQUANT_BOOST_HAS_TRIVIAL_/g'
find tmp/boost -type f -print0 | xargs -0 sed -i "" 's/BOOST_HAS_NOTHROW_/SEQUANT_BOOST_HAS_NOTHROW_/g'
find tmp/boost/version.hpp -type f -print0 | xargs -0 sed -i "" 's/ BOOST_VERSION/ SEQUANT_BOOST_VERSION/g'
find tmp/boost/version.hpp -type f -print0 | xargs -0 sed -i "" 's/ BOOST_LIB_VERSION/ SEQUANT_BOOST_LIB_VERSION/g'
cd tmp
tar -cvzf boost.tar.gz boost/
mv boost.tar.gz ../../../../SeQuant/external
cd ../..

#### cleanup
rm -rf boost_${VERSION_} boost_${VERSION_}.tar.bz2
