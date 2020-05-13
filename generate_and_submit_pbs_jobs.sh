#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=1
	NDATASETS=2
	DIRECTORY=RESULTS_test_NDATASETS
	PYTHIA_WALLTIME='01:00:00'
	HBT_WALLTIME='01:00:00'

	# Set job specifications here
	declare -a specs=(
		'projectile="p" target="p" beamEnergy="13000" chosenHBTparticle="211" Nevents=10000 runHBTEG="false" runFitCF="false" runSV="false"'
	)

	#-----------------------------------------------------
	# Shouldn't have to change anything below this point
	#-----------------------------------------------------
	# Make sure all scripts are executable!
	find . -name "*.sh" | xargs chmod 755

	# Compile source code
	./compile_all.sh $NTHREADS &> compile_all.out

	# Export job specifications
	source scripts/env.sh
	source scripts/export_specs.sh
	SPECS_FULL_PATH=`readlink -f scripts/specs.sh`
	export_specs_array "${specs[@]}" > $SPECS_FULL_PATH

	# Generate jobs for these specifications
	DIRECTORY_FULL_PATH=`readlink -f $DIRECTORY`
	./scripts/generate_pbs_jobs.sh $NTHREADS $NDATASETS $DIRECTORY_FULL_PATH $PYTHIA_WALLTIME $HBT_WALLTIME

	# Submit generated jobs
	./scripts/submit_pbs_jobs.sh $DIRECTORY_FULL_PATH

) &>/dev/null &

# End of file
