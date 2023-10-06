#!/bin/bash

# download and extract tar
wget https://pythia.org/download/pythia83/pythia8309.tgz
tar xvfz pythia8309.tgz

# set source dir
mv pythia8309 pythia8309_src
export PYTHIA_SRC=$(readlink -f pythia8309_src)

# set installation dir
mkdir pythia8309_root
export PYTHIA_ROOT=$(readlink -f pythia8309_root)

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
