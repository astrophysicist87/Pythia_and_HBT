#! /usr/bin/env bash

(
	./generate_jobs.sh 12 RESULTS_pp_N50_100_wCollectivity

	./submit_jobs.sh RESULTS_pp_N50_100_wCollectivity
) &> /dev/null &

# End of file
