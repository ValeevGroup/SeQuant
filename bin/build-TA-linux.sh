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

# Install TA unless previous install with required revision is cached ... still, must manually wipe cache on toolchain update!
if [ ! -d "${INSTALL_DIR}" ]; then

  # Configure TiledArray
  cd ${BUILD_PREFIX}
  mkdir -p TA
  cd TA

  # check out the master of TiledArray
  git clone https://github.com/ValeevGroup/tiledarray.git ta_src

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
      -DCMAKE_PREFIX_PATH="${INSTALL_PREFIX}/boost;${INSTALL_PREFIX}/eigen3" \
      -DCMAKE_TOOLCHAIN_FILE=cmake/vg/toolchains/travis.cmake \
      $EXTRA_CMAKE_FLAGS

  # Build all libraries, examples, and applications
  cmake --build . -j 2
  cmake --build . --target install
fi
