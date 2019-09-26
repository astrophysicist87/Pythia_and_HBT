#! /usr/bin/env bash

(
	./generate_jobs.sh 12 RESULTS_PbPb_5020GeV

	./submit_jobs.sh RESULTS_PbPb_5020GeV
) &> /dev/null &

# End of file
