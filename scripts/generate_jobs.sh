#! /usr/bin/env bash

########################################
# Fix OpenMP settings and compile
chosen_OMP_NUM_THREADS=$1
echo 'export chosen_OMP_NUM_THREADS='$chosen_OMP_NUM_THREADS > $SCRIPTS_DIRECTORY/omp_env.sh

########################################
# set up array of job specifications
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
	cp $HOME_DIRECTORY/parameters.dat $HOME_RESULTS_DIRECTORY/job-${job}/src/HBT_event_generator

	echo "./driver.sh ${specs[i]} &> driver.out" > $HOME_RESULTS_DIRECTORY/job-${job}/submit.sh
	chmod 755 $HOME_RESULTS_DIRECTORY/job-${job}/submit.sh	# set correct permissions!

	cp $SCRIPTS_DIRECTORY/driver.sh				$HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/run_Pythia.sh			$HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/run_HBT_analysis.sh	$HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/rerun.sh				$HOME_RESULTS_DIRECTORY/job-${job}
	cp $SCRIPTS_DIRECTORY/defaults.sh \
		$SCRIPTS_DIRECTORY/specs.sh $SCRIPTS_DIRECTORY/env.sh $SCRIPTS_DIRECTORY/omp_env.sh \
							$HOME_RESULTS_DIRECTORY/job-${job}/scripts
done

# End of file
