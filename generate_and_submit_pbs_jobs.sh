#! /usr/bin/env bash

(
	./generate_pbs_jobs.sh 40 RESULTS

	./submit_pbs_jobs.sh RESULTS
) &> /dev/null &

# End of file
