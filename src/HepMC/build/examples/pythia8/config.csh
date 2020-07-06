#!/bin/csh
if( ! $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/lib
else
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/lib
endif
setenv PYTHIA8DATA ${PYTHIA8_HOME}/xmldoc
