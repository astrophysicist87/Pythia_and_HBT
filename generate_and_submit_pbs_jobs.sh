#! /usr/bin/env bash

(
	# info for this run
	CLUSTER='mangi'
	NTHREADS=128
	NDATASETS=10
	DIRECTORY=RESULTS_test_pp_shifting_modes
	PYTHIA_WALLTIME='96:00:00'
	HBT_WALLTIME='72:00:00'

	# Set job specifications here
	declare -a class_ranges=("1-11" "12-16" "17-22" "23-28" "29-34" "35-41" "42-51" "52-151" "152-1000000")
	#declare -a class_ranges=("0-5%" "5-10%" "10-20%" "20-30%" "30-40%" "40-60%" "60-100%")
	declare -a specs=(
		'projectile="p" target="p" beamEnergy="13000" chosenHBTparticle="211" Nevents=1000000 BEeffects="on" shiftingSet="1" compensationSet="0" compensationMode="1" eventClassSelectionMode="multiplicity" runSV="false"'
		'projectile="p" target="p" beamEnergy="13000" chosenHBTparticle="211" Nevents=1000000 BEeffects="on" shiftingSet="1" compensationSet="2" compensationMode="1" eventClassSelectionMode="multiplicity" runSV="false"'
		'projectile="p" target="p" beamEnergy="13000" chosenHBTparticle="211" Nevents=1000000 BEeffects="on" shiftingSet="0" compensationSet="1" compensationMode="0" eventClassSelectionMode="multiplicity" runSV="false"'
		'projectile="p" target="p" beamEnergy="13000" chosenHBTparticle="211" Nevents=1000000 BEeffects="on" shiftingSet="1" compensationSet="1" compensationMode="0" eventClassSelectionMode="multiplicity" runSV="false"'
	)

	#-----------------------------------------------------
	# Shouldn't have to change anything below this point
	#-----------------------------------------------------
	# Make sure all scripts are executable!
	find . -name "*.sh" | xargs chmod 755

	# Load any relevant defaults for this machine
	source scripts/pbs_env.sh
	setup_env $CLUSTER
	chosen_QUEUENAME="${QUEUENAME:-$def_queuename}"
	chosen_NTHREADS="${NTHREADS:-$def_nthreads}"
	chosen_PYTHIA_WALLTIME="${PYTHIA_WALLTIME:-$def_pythia_walltime}"
	chosen_HBT_WALLTIME="${HBT_WALLTIME:-$def_HBT_walltime}"
	echo '-------------------------' &> compile_all.out
	echo 'Working on' $CLUSTER &>> compile_all.out
	echo $CLUSTER': default queue name =' $def_queuename &>> compile_all.out
	echo $CLUSTER': chosen queue name =' $chosen_QUEUENAME &>> compile_all.out
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
	export_class_range_array "${class_ranges[@]}" >> $SPECS_FULL_PATH

	# Generate jobs for these specifications
	DIRECTORY_FULL_PATH=`readlink -f $DIRECTORY`
	./scripts/generate_pbs_jobs.sh $chosen_NTHREADS $NDATASETS $DIRECTORY_FULL_PATH \
                                   $chosen_QUEUENAME $chosen_PYTHIA_WALLTIME $chosen_HBT_WALLTIME

	# Submit generated jobs
	./scripts/submit_pbs_jobs.sh $DIRECTORY_FULL_PATH

) &>/dev/null &

# End of file
