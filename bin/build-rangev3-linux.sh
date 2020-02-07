#! /bin/sh

# Exit on error
set -ev

export INSTALL_DIR=${INSTALL_PREFIX}/range-v3

# extract the tracked tag of Range-V3
export RANGEV3_TAG=`grep 'set(SEQUANT_TRACKED_RANGESV3_TAG ' ${TRAVIS_BUILD_DIR}/external/versions.cmake | awk '{print $2}' | sed s/\)//g`
echo "required Range-V3 revision = ${RANGEV3_TAG}"

# make sure installed Range-V3 tag matches the required tag, if not, remove INSTALL_DIR (will cause reinstall)
# if can't determine the tag, remove and reinstall just in case
if [ -f "${INSTALL_DIR}/revision.cmake" ]; then
  export INSTALLED_RANGEV3_TAG=`grep 'set(SEQUANT_INSTALLED_RANGESV3_TAG ' ${INSTALL_DIR}/revision.cmake | awk '{print $2}' | sed s/\)//g`
  echo "installed Range-V3 revision = ${INSTALLED_RANGEV3_TAG}"
  if [ "${RANGEV3_TAG}" != "${INSTALLED_RANGEV3_TAG}" ]; then
    rm -rf "${INSTALL_DIR}"
  fi
else
  rm -rf "${INSTALL_DIR}"
fi

# Clone Range-V3 unless previous install is cached ... must manually wipe cache on version bump or toolchain update
if [ ! -d "${INSTALL_DIR}" ]; then
    cd ${INSTALL_PREFIX}
    git clone https://github.com/ericniebler/range-v3.git
    cd range-v3 && git checkout ${RANGEV3_TAG} && cd ..
    echo "set(SEQUANT_INSTALLED_RANGESV3_TAG ${RANGEV3_TAG})" > revision.cmake
else
    echo "Range-V3 already installed ..."
fi
