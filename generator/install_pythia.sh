#!/bin/bash

# download and extract tar
wget https://pythia.org/download/pythia83/pythia8310.tgz
tar xvfz pythia8310.tgz

# set source dir
mv pythia8310 pythia8310_src
export PYTHIA_SRC=$(readlink -f pythia8310_src)

# set installation dir
mkdir pythia8310
export PYTHIA_ROOT=$(readlink -f pythia8310)

# set current dir
export CURRENT_DIR=$(readlink -f ${PWD})

# build
cd ${PYTHIA_SRC}
bash configure --prefix=${PYTHIA_ROOT} # --with-hepmc3=${HEPMC3_ROOT} # CREATED: bin/ Makefile.inc
make -j 8 # CREATED: examples tmp lib
make install

# come back to current dir
cd ${CURRENT_DIR}

# copy Makefiles
cp ${PYTHIA_ROOT}/share/Pythia8/examples/Makefile .
cp ${PYTHIA_ROOT}/share/Pythia8/examples/Makefile.inc .
