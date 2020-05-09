#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=12
	DIRECTORY=RESULTS

	# Set job specifications here
	declare -a specs=(
		'projectile="Pb" target="Pb" beamEnergy="2760" chosenHBTparticle="211" eventClassSelectionMode="centrality" Nevents=25000 storeBjorkenCoordinates="false" BEeffects="off" runSV="false"'
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
	./scripts/generate_jobs.sh $NTHREADS $DIRECTORY_FULL_PATH

	# Submit generated jobs
	./scripts/submit_jobs.sh $DIRECTORY_FULL_PATH

) &> /dev/null &

# End of file
