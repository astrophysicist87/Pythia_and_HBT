#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=40
	DIRECTORY=RESULTS

	./generate_pbs_jobs.sh $NTHREADS $DIRECTORY

	./submit_pbs_jobs.sh $DIRECTORY
) &> /dev/null &

# End of file
