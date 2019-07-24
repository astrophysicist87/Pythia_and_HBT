#! /usr/bin/env bash

(
	./generate_jobs.sh 12 RESULTS

	./submit_jobs.sh RESULTS
) &> /dev/null &

# End of file
