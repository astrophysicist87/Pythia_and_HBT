#! /usr/bin/env bash

source $SCRIPTS_DIRECTORY/pbs_env.sh

# Set PBS walltime directives
export chosen_Pythia_walltime=$3
export chosen_HBT_walltime_per_centrality=$4

########################################
# Fix OpenMP settings and compile
chosen_OMP_NUM_THREADS=$1
echo 'export chosen_OMP_NUM_THREADS='$chosen_OMP_NUM_THREADS > $SCRIPTS_DIRECTORY/omp_env.sh

source $SCRIPTS_DIRECTORY/specs.sh

########################################
# total number of jobs
nJobs=${#specs[@]}

# make sure main results directory exists
HOME_RESULTS_DIRECTORY=$2
if [ ! -d "$HOME_RESULTS_DIRECTORY" ]
then
	mkdir $HOME_RESULTS_DIRECTORY
	#echo 'Created' $HOME_RESULTS_DIRECTORY
fi

# loop over jobs
for ((i=0; i<$nJobs; i++))
do
	job=$[i+1]
	mkdir $HOME_RESULTS_DIRECTORY/job-${job}
	mkdir $HOME_RESULTS_DIRECTORY/job-${job}/scripts
	cp -r $SOURCE_DIRECTORY $HOME_RESULTS_DIRECTORY/job-${job}/
	
	executableString="./driver_pbs.sh ${specs[i]} &> driver_pbs.out"
	generate_pbs $chosen_Pythia_walltime $chosen_OMP_NUM_THREADS $executableString > $HOME_RESULTS_DIRECTORY/job-${job}/submit.pbs
	
	cp $SCRIPTS_DIRECTORY/driver_pbs.sh			$HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/run_Pythia.sh			$HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/run_HBT_analysis.pbs	$HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/rerun_pbs.sh			$HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/defaults.sh $SCRIPTS_DIRECTORY/specs.sh \
		$SCRIPTS_DIRECTORY/env.sh $SCRIPTS_DIRECTORY/omp_env.sh \
		$SCRIPTS_DIRECTORY/pbs_env.sh 			$HOME_RESULTS_DIRECTORY/job-${job}/scripts
done

# End of file
