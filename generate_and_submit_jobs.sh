#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=1
	DIRECTORY=RESULTS_test_new_format

	./generate_jobs.sh $NTHREADS $DIRECTORY

	./submit_jobs.sh $DIRECTORY
) &> /dev/null &

# End of file
