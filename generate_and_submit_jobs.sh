#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=12
	DIRECTORY=RESULTS_pp_shiftMode_1_test_Nev100000_woLinearInterp_newCompensation

	./generate_jobs.sh $NTHREADS $DIRECTORY

	./submit_jobs.sh $DIRECTORY
) &> /dev/null &

# End of file
