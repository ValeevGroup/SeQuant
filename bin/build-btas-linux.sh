#! /bin/sh

# Exit on error
set -ev

# Clone BTAS unless previous install is cached ... must manually wipe cache on version bump or toolchain update
export INSTALL_DIR=${INSTALL_PREFIX}/BTAS
if [ ! -d "${INSTALL_DIR}" ]; then
    cd ${INSTALL_PREFIX}
    git clone https://github.com/BTAS/BTAS.git
else
    echo "BTAS already installed ..."
fi
