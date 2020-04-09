#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=1
	DIRECTORY=RESULTS_test_new_format

	# Set job specifications here
	declare -a specs=(
		'useArbitraryParticle=true projectile="p" target="p" beamEnergy="7000.0" chosenHBTparticle="211" Nevents=10 storeBjorkenCoordinates="false" BEeffects="off"'
	)

	#-----------------------------------------------------
	# Shouldn't have to change anything below this point
	#-----------------------------------------------------
	source scripts/export_specs.sh
	export_specs_array "${array[@]}" > scripts/specs.sh

	# Generate jobs for these specifications
	./generate_jobs.sh $NTHREADS $DIRECTORY

	# Submit generated jobs
	./submit_jobs.sh `readlink -f $DIRECTORY`
) &> /dev/null &

# End of file
