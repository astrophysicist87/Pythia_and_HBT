#! /usr/bin/env bash

(
	# info for this run
	CLUSTER='mesabi'
	NTHREADS=1
	NDATASETS=1
	DIRECTORY=RESULTS_test_pp_shifting_modes
	PYTHIA_WALLTIME='01:00:00'
	HBT_WALLTIME='01:00:00'

	# Set job specifications here
	declare -a specs=(
		'projectile="p" target="p" beamEnergy="13000" Nevents=1000'
		#'projectile="p" target="p" beamEnergy="13000" chosenHBTparticle="211" Nevents=1000000 BEeffects="on" shiftingSet="1" compensationSet="0" compensationMode="1" eventClassSelectionMode="multiplicity" runSV="false"'
		#'projectile="p" target="p" beamEnergy="13000" chosenHBTparticle="211" Nevents=1000000 BEeffects="on" shiftingSet="1" compensationSet="2" compensationMode="1" eventClassSelectionMode="multiplicity" runSV="false"'
		#'projectile="p" target="p" beamEnergy="13000" chosenHBTparticle="211" Nevents=1000000 BEeffects="on" shiftingSet="0" compensationSet="1" compensationMode="0" eventClassSelectionMode="multiplicity" runSV="false"'
		#'projectile="p" target="p" beamEnergy="13000" chosenHBTparticle="211" Nevents=1000000 BEeffects="on" shiftingSet="1" compensationSet="1" compensationMode="0" eventClassSelectionMode="multiplicity" runSV="false"'
	)

	#-----------------------------------------------------
	# Shouldn't have to change anything below this point
	#-----------------------------------------------------
	# Make sure all scripts are executable!
	find . -name "*.sh" | xargs chmod 755

	# Load any relevant defaults for this machine
	source scripts/pbs_env.sh
	setup_env $CLUSTER
	chosen_NTHREADS="${NTHREADS:-$def_nthreads}"
	chosen_PYTHIA_WALLTIME="${PYTHIA_WALLTIME:-$def_pythia_walltime}"
	chosen_HBT_WALLTIME="${HBT_WALLTIME:-$def_HBT_walltime}"
	echo '-------------------------' &> compile_all.out
	echo 'Working on' $CLUSTER &>> compile_all.out
	echo $CLUSTER': default number of threads =' $def_nthreads &>> compile_all.out
	echo $CLUSTER': chosen number of threads =' $chosen_NTHREADS &>> compile_all.out
	echo $CLUSTER': default Pythia walltime =' $def_pythia_walltime &>> compile_all.out
	echo $CLUSTER': chosen Pythia walltime =' $chosen_PYTHIA_WALLTIME &>> compile_all.out
	echo $CLUSTER': default HBT walltime =' $def_HBT_walltime &>> compile_all.out
	echo $CLUSTER': chosen HBT walltime =' $chosen_HBT_WALLTIME &>> compile_all.out
	echo '-------------------------' &>> compile_all.out

	# Compile source code
	./compile_all.sh $chosen_NTHREADS &>> compile_all.out

	# Export job specifications
	source scripts/env.sh
	source scripts/export_specs.sh
	SPECS_FULL_PATH=`readlink -f scripts/specs.sh`
	export_specs_array "${specs[@]}" > $SPECS_FULL_PATH

	# Generate jobs for these specifications
	DIRECTORY_FULL_PATH=`readlink -f $DIRECTORY`
	./scripts/generate_pbs_jobs.sh $chosen_NTHREADS $NDATASETS $DIRECTORY_FULL_PATH $chosen_PYTHIA_WALLTIME $chosen_HBT_WALLTIME

	# Submit generated jobs
	./scripts/submit_pbs_jobs.sh $DIRECTORY_FULL_PATH

) &>/dev/null &

# End of file
