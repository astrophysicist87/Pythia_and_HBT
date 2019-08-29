#! /usr/bin/env bash

(
	./generate_jobs.sh 12 RESULTS_PbPb

	./submit_jobs.sh RESULTS_PbPb
) &> /dev/null &

# End of file
