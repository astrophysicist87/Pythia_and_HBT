#! /usr/bin/env bash

(
	./generate_jobs.sh 16 RESULTS_PbPb

	./submit_jobs.sh RESULTS_PbPb
) &> /dev/null &

# End of file
