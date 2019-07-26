#! /usr/bin/env bash

(
	./generate_jobs.sh 12 RESULTS_pPb

	./submit_jobs.sh RESULTS_pPb
) &> /dev/null &

# End of file
