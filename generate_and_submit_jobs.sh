#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=1
	DIRECTORY=RESULTS_pp_13TeV_v8244/

	./generate_jobs.sh $NTHREADS $DIRECTORY

	./submit_jobs.sh $DIRECTORY
) &> /dev/null &

# End of file
