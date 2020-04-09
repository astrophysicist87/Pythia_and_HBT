#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=40
	DIRECTORY=RESULTS2
	PYTHIA_WALLTIME='12:00:00'
	HBT_WALLTIME='48:00:00'

	# Set job specifications here
	declare -a specs=(
		'useArbitraryParticle=true projectile="p" target="p" beamEnergy="7000.0" chosenHBTparticle="211" Nevents=10000000 storeBjorkenCoordinates="false" BEeffects="off"'
		'useArbitraryParticle=true projectile="p" target="Pb" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=1000000 storeBjorkenCoordinates="false" BEeffects="off"'
		'useArbitraryParticle=true projectile="Pb" target="Pb" beamEnergy="2760.0" chosenHBTparticle="211" Nevents=100000 storeBjorkenCoordinates="false" BEeffects="off"'
	)

	#-----------------------------------------------------
	# Shouldn't have to change anything below this point
	#-----------------------------------------------------
	source scripts/export_specs.sh
	export_specs_array "${array[@]}" > scripts/specs.sh

	# Generate jobs for these specifications
	./scripts/generate_pbs_jobs.sh $NTHREADS $DIRECTORY $PYTHIA_WALLTIME $HBT_WALLTIME

	# Submit generated jobs
	./scripts/submit_pbs_jobs.sh `readlink -f $DIRECTORY`
) &> /dev/null &

# End of file
