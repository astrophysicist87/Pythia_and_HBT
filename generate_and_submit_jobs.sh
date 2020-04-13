#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=12
	DIRECTORY=RESULTS

	# Set job specifications here
	declare -a specs=(
		'useArbitraryParticle=true projectile="Pb" target="Pb" beamEnergy="2760.0" chosenHBTparticle="211" Nevents=50 storeBjorkenCoordinates="false" BEeffects="off"'
	)

	#-----------------------------------------------------
	# Shouldn't have to change anything below this point
	#-----------------------------------------------------
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
