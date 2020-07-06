#!/bin/sh
if [ ! $?LD_LIBRARY_PATH ]; then
  export LD_LIBRARY_PATH=/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/lib
fi
if [ $?LD_LIBRARY_PATH ]; then
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/lib
fi
export PYTHIA8DATA=${PYTHIA8_HOME}/xmldoc
