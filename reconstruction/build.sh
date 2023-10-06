#!/bin/bash

if [[ ! $(geant4-config) ]]; then
    echo "ERROR: make sure to set GEANT4"
    exit 1
fi

CURRENT_DIR=${PWD}

# get geant4 lib dir
G4_LIB_DIR=$(readlink -f "$(geant4-config --prefix)/lib")

#define source dir
SOURCE_DIR=$(readlink -f ./G4Setup)

# define and create build dir
BUILD_DIR=$(readlink -f ./G4Setup_build)
mkdir -p ${BUILD_DIR}

cd ${BUILD_DIR}

# if makefile doesn't exist
if [[ ! -e ${BUILD_DIR}/Makefile ]]; then
    cmake -DGeant4_DIR=${G4_LIB_DIR} ${SOURCE_DIR}
fi

# if executable already exists
if [[ -e ${BUILD_DIR}/exampleB2a ]]; then
    make clean
fi

make -j 8

cd ${CURRENT_DIR}
