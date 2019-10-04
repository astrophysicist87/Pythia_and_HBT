#! /usr/bin/env bash

(
	./generate_jobs.sh 1 RESULTS_AuAu_200GeV

	./submit_jobs.sh RESULTS_AuAu_200GeV
) &> /dev/null &

# End of file
