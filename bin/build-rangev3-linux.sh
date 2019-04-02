#! /bin/sh

# Exit on error
set -ev

# Clone Range-V3 unless previous install is cached ... must manually wipe cache on version bump or toolchain update
export INSTALL_DIR=${INSTALL_PREFIX}/range-v3
if [ ! -d "${INSTALL_DIR}" ]; then
    cd ${INSTALL_PREFIX}
    git clone https://github.com/ericniebler/range-v3.git
else
    echo "Range-V3 already installed ..."
fi
