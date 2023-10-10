#!/bin/bash

CURRENT_DIR=${PWD}

# define build dir
BUILD_DIR=$(readlink -f ./G4Setup_build)

cd ${BUILD_DIR}

./main

cd ${CURRENT_DIR}
