#! /usr/bin/env bash

source $SCRIPTS_DIRECTORY/pbs_env.sh

# Set PBS directives
export chosen_QUEUENAME=$4
export chosen_Pythia_walltime=$5
export chosen_HBT_walltime_per_event_class=$6

########################################
# Fix OpenMP settings and compile
chosen_OMP_NUM_THREADS=$1
echo 'export chosen_OMP_NUM_THREADS='$chosen_OMP_NUM_THREADS > $SCRIPTS_DIRECTORY/omp_env.sh

source $SCRIPTS_DIRECTORY/specs.sh

########################################
# total number of jobs
nJobs=${#specs[@]}

# make sure main results directory exists
HOME_RESULTS_DIRECTORY=$3
if [ ! -d "$HOME_RESULTS_DIRECTORY" ]
then
	mkdir $HOME_RESULTS_DIRECTORY
	#echo 'Created' $HOME_RESULTS_DIRECTORY
fi

# fix number of datasets to use
NDATASETS=$2

# loop over jobs
for ((i=0; i<$nJobs; i++))
do
	job=$[i+1]
	mkdir $HOME_RESULTS_DIRECTORY/job-${job}
	mkdir $HOME_RESULTS_DIRECTORY/job-${job}/scripts
	cp -r $SOURCE_DIRECTORY $HOME_RESULTS_DIRECTORY/job-${job}/
	cp $HOME_DIRECTORY/parameters.dat $HOME_RESULTS_DIRECTORY/job-${job}/src/HBT_event_generator
	cp $HOME_DIRECTORY/parameters.dat $HOME_RESULTS_DIRECTORY/job-${job}
	cp $HOME_DIRECTORY/parameters.dat $PYTHIA_DIRECTORY

	# Before submitting Pythia jobs, clear out any *catalogue.dat files
	find $HOME_RESULTS_DIRECTORY/job-${job} -name "*catalogue.dat" | xargs rm

	# Generate submission scripts
	executableString="./driver_pbs.sh NDATASETS=$NDATASETS chosen_QUEUENAME=$chosen_QUEUENAME chosen_Pythia_walltime=$chosen_Pythia_walltime chosen_HBT_walltime_per_event_class=$chosen_HBT_walltime_per_event_class ${specs[i]} &> driver_pbs.out"
	generate_sh $executableString > $HOME_RESULTS_DIRECTORY/job-${job}/submit.sh

	cp $SCRIPTS_DIRECTORY/driver_pbs.sh           $HOME_RESULTS_DIRECTORY/job-${job}
	#cp $SCRIPTS_DIRECTORY/run_Pythia.sh           $HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/run_Pythia.pbs          $HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/run_HBT_analysis.pbs    $HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/post_process_Pythia.pbs $HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/rerun_pbs.sh            $HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/defaults.sh \
       $SCRIPTS_DIRECTORY/specs.sh    \
       $SCRIPTS_DIRECTORY/env.sh      \
       $SCRIPTS_DIRECTORY/omp_env.sh  \
       $SCRIPTS_DIRECTORY/pbs_env.sh              $HOME_RESULTS_DIRECTORY/job-${job}/scripts
done

# End of file
