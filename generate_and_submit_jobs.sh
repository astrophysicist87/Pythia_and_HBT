#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=4
	DIRECTORY=RESULTS_pp_Nev10000000_MB_wBA_noE1E2

	./generate_jobs.sh $NTHREADS $DIRECTORY

	./submit_jobs.sh $DIRECTORY
) &> /dev/null &

# End of file
