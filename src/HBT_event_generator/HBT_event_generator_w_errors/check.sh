#! /usr/bin/env bash

export OMP_NUM_THREADS=12

rm -rf results checks
mkdir results checks

i=1

# time and run
nohup time ./run_checks.e \
		bin_epsilon=0.000025 \
		n_KT_pts=2 \
		KTmin=0.0 \
		KTmax=10000000.0 \
		n_KL_pts=2 \
		KLmin=-10000000.0 \
		KLmax=10000000.0 \
		n_qo_pts=22 \
		n_qs_pts=22 \
		n_ql_pts=22 \
		n_Q_pts=4 \
		q_mode=1 \
		scalar_mode=1 \
		method_mode=2 \
		BE_mode=0 \
		RNG_Nev=100 \
		RNG_xDir=1 \
		RNG_yDir=1 \
		RNG_zDir=1 \
		RNG_mult=100000 \
		RNG_nLoops=10 \
		1> check_HBT_event_generator.out \
		2> check_HBT_event_generator.err

DIRECTORY=checks/results-${i}
mv results $DIRECTORY
mv check_HBT_event_generator.* $DIRECTORY
cp parameters.dat $DIRECTORY

# End of file
