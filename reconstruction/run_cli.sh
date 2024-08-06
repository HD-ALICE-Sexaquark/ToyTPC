#!/bin/bash

CURRENT_DIR=${PWD}

# define build dir
BUILD_DIR=$(readlink -f ./G4Setup_build)

cd ${BUILD_DIR}

./main ../mc.csv ../traj.csv ../its.csv ../tpc.csv 0.2 2>&1 | tee ../rec.log

cd ${CURRENT_DIR}
