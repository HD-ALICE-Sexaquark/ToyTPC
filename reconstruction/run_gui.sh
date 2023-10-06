#!/bin/bash

BKG_PDG_CODE=${1}

CURRENT_DIR=${PWD}

# define build dir
BUILD_DIR=$(readlink -f ./G4Setup_build)

cd ${BUILD_DIR}

./main ${BKG_PDG_CODE}

cd ${CURRENT_DIR}
