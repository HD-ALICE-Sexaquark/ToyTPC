#!/bin/bash

CURRENT_DIR=${PWD}

# define build dir
BUILD_DIR=$(readlink -f ./G4Setup_build)

cd ${BUILD_DIR}

# GEANT4_NUM_THREADS=1 ./main ../mc.csv ../traj_noTC.csv ../its_noTC.csv ../tpc_noTC.csv 0 0.2 2>&1 | tee ../rec_noTC.log
GEANT4_NUM_THREADS=1 ./main ../mc.csv ../traj.csv ../its.csv ../tpc.csv 1 0.2 2>&1 | tee ../rec.log

cd ${CURRENT_DIR}
