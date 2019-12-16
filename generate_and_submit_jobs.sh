#! /usr/bin/env bash

(
	# info for this run
	#NTHREADS=12
	#DIRECTORY=RESULTS_pp_13TeV_v8243_QRef_RMS_SpatialSep_inRF/
	NTHREADS=1
	DIRECTORY=RESULTS_pp_shiftMode_2_test

	./generate_jobs.sh $NTHREADS $DIRECTORY

	./submit_jobs.sh $DIRECTORY
) &> /dev/null &

# End of file
