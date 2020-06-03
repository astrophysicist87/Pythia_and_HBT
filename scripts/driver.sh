#! /usr/bin/env bash

echo '--------------------------------------------------------------------------------'	# width is 80 spaces
echo '| '`basename "$0"`': Getting started!'

# Load header info
echo '| '`basename "$0"`': loading scripts/env.sh...'
source scripts/env.sh

# Load default variable values
echo '| '`basename "$0"`': loading scripts/defaults.sh...'
source scripts/defaults.sh

# Load OpenMP info
echo '| '`basename "$0"`': loading scripts/omp_env.sh...'
source scripts/omp_env.sh
export OMP_NUM_THREADS=$chosen_OMP_NUM_THREADS
echo '| '`basename "$0"`': OMP_NUM_THREADS =' $OMP_NUM_THREADS

# Load job specifications and event class ranges
echo '| '`basename "$0"`': loading scripts/specs.sh...'
source scripts/specs.sh

# Update any variables set from the command line
echo '| '`basename "$0"`': loading settings from command-line...'
for var in "$@"
do
    export "$var"
done

# Save the settings this job was run with (for future defaults)
echo '| '`basename "$0"`': saving all loaded settings to settings.sh.'
output_settings > settings.sh


echo '| '`basename "$0"`': storing all job specifications in this_job.txt.'
output_this_job > this_job.txt


echo '| '`basename "$0"`': everything being run from '$HOSTNAME'.'


#if [ -z ${chosen_OMP_NUM_THREADS+x} ]
#then
#	echo "chosen_OMP_NUM_THREADS is unset"
#else
#	echo "chosen_OMP_NUM_THREADS is set to '$chosen_OMP_NUM_THREADS'"
#fi


# truth value of $runPythia evaluated inside
# note command-line argument corresponding to seed=0
echo '| '`basename "$0"`': running Pythia!'
./run_Pythia.sh 0



# apply HBT analysis to each chosen event class
for eventClassCutString in "${class_ranges[@]}"
do
	echo '| '`basename "$0"`': running HBT analysis for event class = '$eventClassCutString'!'
	./run_HBT_analysis.sh $eventClassCutString	# do NOT submit Bash scripts in background
done	# all centralities finished


#PYTHIA_RESULTS_DIRECTORY=$MAIN_RESULTS_DIRECTORY/Pythia_results/dataset_0
#zipFilename=$PYTHIA_RESULTS_DIRECTORY".zip"

#zip -r $zipFilename $PYTHIA_RESULTS_DIRECTORY && rm -rf $PYTHIA_RESULTS_DIRECTORY


echo '| '`basename "$0"`': Finished everything!'
echo '--------------------------------------------------------------------------------'	# width is 80 spaces


# End of file
