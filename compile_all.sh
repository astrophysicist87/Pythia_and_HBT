#! /usr/bin/env bash

#=====================================
# Load header info
#source env.sh
source scripts/env.sh

default_num_threads=1
export OMP_NUM_THREADS="${1:-$default_num_threads}"
echo 'OMP_NUM_THREADS =' $OMP_NUM_THREADS

success=0

#=====================================
# Compile Pythia
echo '#====================================='
echo '# Compiling Pythia'
echo '#====================================='
cd $PYTHIA_DIRECTORY
echo 'In directory='`pwd`':'
echo '#====================================='
rm main_BEeffects_OpenMP main_BEeffects_arbtryParticle_OpenMP
make main_BEeffects_OpenMP main_BEeffects_arbtryParticle_OpenMP
success=$[success+`echo $?`]

#=====================================
# Compile HBT_event_generator
echo '#====================================='
echo '# Compiling HBT_event_generator'
echo '#====================================='
cd $HBT_EVENT_GEN_DIRECTORY
echo 'In directory='`pwd`':'
echo '#====================================='
gmake distclean
gmake target
success=$[success+`echo $?`]

#=====================================
# Compile fit_correlation_function
echo '#====================================='
echo '# Compiling fit_correlation_function'
echo '#====================================='
cd $HBT_FITCF_DIRECTORY
echo 'In directory='`pwd`':'
echo '#====================================='
gmake distclean
gmake all
success=$[success+`echo $?`]
echo '#====================================='

#=====================================
# Compile SV.e
echo '#====================================='
echo '# Compiling SV.e'
echo '#====================================='
cd $HBT_SV_DIRECTORY
echo 'In directory='`pwd`':'
echo '#====================================='
gmake distclean
gmake all
success=$[success+`echo $?`]
echo '#====================================='

#=====================================
# Check success
cd $HOME_DIRECTORY
echo 'In directory='`pwd`':'
echo '#====================================='
if [ "$success" -eq "0" ]
then
	echo 'END RESULT: Everything compiled successfully!'
else
	echo 'END RESULT: There were problems compiling.'
	exit $success
fi

# End of file
