#! /usr/bin/env bash

source ../../omp_env.sh
export OMP_NUM_THREADS=$chosen_OMP_NUM_THREADS
echo 'OMP_NUM_THREADS =' $OMP_NUM_THREADS

# make sure results directory exists
DIRECTORY="$5"
if [ ! -d "$DIRECTORY" ]; then
	mkdir -p $DIRECTORY
fi

echo 'Running ./main_BEeffects_OpenMP' $@

# time and run
nohup time valgrind ./main_BEeffects_OpenMP $@\
			1> $DIRECTORY/main_BEeffects_OpenMP.out\
			2> $DIRECTORY/main_BEeffects_OpenMP.err
