#! /usr/bin/env bash

(
	./generate_jobs.sh 12 RESULTS_PbPb_5020GeV_kaons

	./submit_jobs.sh RESULTS_PbPb_5020GeV_kaons
) &> /dev/null &

# End of file
