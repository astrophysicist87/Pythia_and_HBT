#! /usr/bin/env bash

(
	# info for this run
	CLUSTER='qgp'
	NTHREADS=12
	NDATASETS=1
	DIRECTORY=PbPb2760GeV
	PYTHIA_WALLTIME='12:00:00'
	HBT_WALLTIME='12:00:00'


	# Set job specifications here
	declare -a class_ranges=("0-100%")
	declare -a specs=(
		'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=100 BEeffects="on" BEEnhancementMode="0" useDistribution="on" useHepMCOutputFormat="false"'
	)

	#-----------------------------------------------------
	# Shouldn't have to change anything below this point
	#-----------------------------------------------------
	# Make sure all scripts are executable!
	find . -name "*.sh" | xargs chmod 755

	# Load any relevant defaults for this machine
	source scripts/slurm_env.sh
	setup_env $CLUSTER
	chosen_QUEUENAME="${QUEUENAME:-$def_queuename}"
	chosen_NTHREADS="${NTHREADS:-$def_nthreads}"
	chosen_PYTHIA_WALLTIME="${PYTHIA_WALLTIME:-$def_pythia_walltime}"
	chosen_HBT_WALLTIME="${HBT_WALLTIME:-$def_HBT_walltime}"
	echo '-------------------------'                                     &>  compile_all.out
	echo 'Working on' $CLUSTER                                           &>> compile_all.out
	echo $CLUSTER': default queue name        =' $def_queuename          &>> compile_all.out
	echo $CLUSTER': chosen queue name         =' $chosen_QUEUENAME       &>> compile_all.out
	echo '--------'                                                      &>> compile_all.out
	echo $CLUSTER': default number of threads =' $def_nthreads           &>> compile_all.out
	echo $CLUSTER': chosen number of threads  =' $chosen_NTHREADS        &>> compile_all.out
	echo '--------'                                                      &>> compile_all.out
	echo $CLUSTER': default Pythia walltime   =' $def_pythia_walltime    &>> compile_all.out
	echo $CLUSTER': chosen Pythia walltime    =' $chosen_PYTHIA_WALLTIME &>> compile_all.out
	echo '--------'                                                      &>> compile_all.out
	echo $CLUSTER': default HBT walltime      =' $def_HBT_walltime       &>> compile_all.out
	echo $CLUSTER': chosen HBT walltime       =' $chosen_HBT_WALLTIME    &>> compile_all.out
	echo '-------------------------'                                     &>> compile_all.out

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
	./scripts/generate_slurm_jobs.sh $chosen_NTHREADS $NDATASETS $DIRECTORY_FULL_PATH \
                                   $chosen_QUEUENAME $chosen_PYTHIA_WALLTIME $chosen_HBT_WALLTIME

	# Submit generated jobs
	./scripts/submit_slurm_jobs.sh $DIRECTORY_FULL_PATH

) &>/dev/null &

# End of file
