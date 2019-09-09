#! /usr/bin/env bash

(
	./generate_jobs.sh 1 RESULTS_pp_SV_noDecay

	./submit_jobs.sh RESULTS_pp_SV_noDecay
) &> /dev/null &

# End of file
