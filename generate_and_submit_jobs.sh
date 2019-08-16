#! /usr/bin/env bash

(
	./generate_jobs.sh 12 RESULTS_PbPb_N50_100

	./submit_jobs.sh RESULTS_PbPb_N50_100
) &> /dev/null &

# End of file
