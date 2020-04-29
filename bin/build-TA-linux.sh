#! /bin/sh

# Exit on error
set -ev

# ninja for some reason fails to compile MADNESS although changes are seemingly mundane
# https://travis-ci.com/ValeevGroup/mpqc4/jobs/279131477
export CMAKE_GENERATOR_TA="Unix Makefiles"

# Environment variables
export EXTRA_CMAKE_FLAGS=""
if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
    export EXTRACXXFLAGS="-mno-avx"
else
    export CC=/usr/bin/clang-$CLANG_VERSION
    export CXX=/usr/bin/clang++-$CLANG_VERSION
    export EXTRACXXFLAGS="-mno-avx -stdlib=libc++"
fi
# need Fortran for Elemental and to help BTAS find lapacke based on netlib
export F77=gfortran-$GCC_VERSION

echo $($CC --version)
echo $($CXX --version)

# cmake
export PATH=${INSTALL_PREFIX}/cmake/bin:${PATH}
cmake --version

export MPI_HOME=${INSTALL_PREFIX}/mpich
export MPICC=$MPI_HOME/bin/mpicc
export MPICXX=$MPI_HOME/bin/mpicxx
export LD_LIBRARY_PATH=/usr/lib/lapack:/usr/lib/libblas:$LD_LIBRARY_PATH
export INSTALL_DIR=${INSTALL_PREFIX}/TA

# extract the tracked tag of TA
export TA_TAG=`grep 'set(MPQC_TRACKED_TILEDARRAY_TAG ' ${TRAVIS_BUILD_DIR}/external/versions.cmake | awk '{print $2}' | sed s/\)//g`
echo "required TA revision = ${TA_TAG}"

# make sure installed TiledArray tag matches the required tag, if not, remove INSTALL_DIR (will cause reinstall)
if [ -f "${INSTALL_DIR}/include/TiledArray/version.h" ]; then
  export INSTALLED_TA_TAG=`grep 'define TILEDARRAY_REVISION' ${INSTALL_DIR}/include/TiledArray/version.h | awk '{print $3}' | sed s/\"//g`
  echo "installed TA revision = ${INSTALLED_TA_TAG}"
  if [ "${TA_TAG}" != "${INSTALLED_TA_TAG}" ]; then
    rm -rf "${INSTALL_DIR}"
  fi
fi

# Install TA unless previous install with required revision is cached ... still, must manually wipe cache on toolchain update!
if [ ! -d "${INSTALL_DIR}" ]; then

  # Configure TiledArray
  cd ${BUILD_PREFIX}
  mkdir -p TA
  cd TA

  # check out the tracked version of TiledArray
  git clone https://github.com/ValeevGroup/tiledarray.git ta_src
  cd ta_src && git checkout ${TA_TAG} && cd ..

  # always compile Elemental in Release mode to avoid non-reentrancy problems
  cmake ta_src \
      -G "${CMAKE_GENERATOR_TA}" \
      -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_Fortran_COMPILER=$F77 \
      -DMPI_CXX_COMPILER=$MPICXX \
      -DMPI_C_COMPILER=$MPICC \
      -DBUILD_SHARED_LIBS=${BUILD_SHARED} \
      -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
      -DCMAKE_CXX_FLAGS="${EXTRACXXFLAGS}" \
      -DCMAKE_PREFIX_PATH=${INSTALL_PREFIX} \
      -DENABLE_ELEMENTAL=ON \
      -DMADNESS_CMAKE_EXTRA_ARGS="-DELEMENTAL_CMAKE_BUILD_TYPE=Release;-DELEMENTAL_MATH_LIBS='-L/usr/lib/libblas -L/usr/lib/lapack -lblas -llapack';-DELEMENTAL_CMAKE_EXTRA_ARGS=-DCMAKE_Fortran_COMPILER=$F77" \
      -DCMAKE_TOOLCHAIN_FILE=cmake/vg/toolchains/travis.cmake \
      $EXTRA_CMAKE_FLAGS

  # Build all libraries, examples, and applications
  cmake --build . -j 2
  cmake --build . --target install
fi
