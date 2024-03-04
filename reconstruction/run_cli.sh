#!/bin/bash

CURRENT_DIR=${PWD}

# define build dir
BUILD_DIR=$(readlink -f ./G4Setup_build)

cd ${BUILD_DIR}

./main ../mc.csv ../traj.csv ../its.csv ../tpc.csv

cd ${CURRENT_DIR}
