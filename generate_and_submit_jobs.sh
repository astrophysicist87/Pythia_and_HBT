#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=12
	DIRECTORY=RESULTS_pp_Nev1000_wBEeffects_DUMMY

	./generate_jobs.sh $NTHREADS $DIRECTORY

	./submit_jobs.sh $DIRECTORY
) &> /dev/null &

# End of file
