#! /usr/bin/env bash

export OMP_NUM_THREADS=12

rm -rf results checks
mkdir results checks
DIRECTORY=`get_dirname checks/results "-" true`
mv results $DIRECTORY

# time and run
nohup time ./run_checks.e \
		1> HBT_event_generator.out \
		2> HBT_event_generator.err

cp ../parameters.dat .
cp ../parameters.dat $DIRECTORY

# End of file
