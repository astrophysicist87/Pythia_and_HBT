#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=1
	DIRECTORY=RESULTS_pp_7TeV
	PYTHIA_WALLTIME='01:00:00'
	HBT_WALLTIME='01:00:00'

	# Set job specifications here
	declare -a specs=(
		'useArbitraryParticle=true projectile="p" target="p" beamEnergy="7000.0" chosenHBTparticle="211" eventClassSelectionMode="multiplicity" Nevents=250 storeBjorkenCoordinates="false" BEeffects="off" runSV="false"'
		#'useArbitraryParticle=true projectile="p" target="Pb" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=1000000 storeBjorkenCoordinates="false" BEeffects="off" runSV="false'
		#'useArbitraryParticle=true projectile="Pb" target="Pb" beamEnergy="2760.0" chosenHBTparticle="211" Nevents=10000 storeBjorkenCoordinates="false" BEeffects="off"'
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
	./scripts/generate_pbs_jobs.sh $NTHREADS $DIRECTORY_FULL_PATH $PYTHIA_WALLTIME $HBT_WALLTIME

	# Submit generated jobs
	./scripts/submit_pbs_jobs.sh $DIRECTORY_FULL_PATH

) &> /dev/null &

# End of file
