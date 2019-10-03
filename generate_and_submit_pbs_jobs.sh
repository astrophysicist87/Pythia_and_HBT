#! /usr/bin/env bash

(
	./generate_pbs_jobs.sh 40 RESULTS_PbPb_5020GeV_pions

	./submit_pbs_jobs.sh RESULTS_PbPb_5020GeV_pions
) &> /dev/null &

# End of file
