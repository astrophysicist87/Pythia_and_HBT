#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=12
	DIRECTORY=RESULTS_pp_Nev10000000_wBA

	./generate_jobs.sh $NTHREADS $DIRECTORY

	./submit_jobs.sh $DIRECTORY
) &> /dev/null &

# End of file
