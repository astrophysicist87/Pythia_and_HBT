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

# Update any variables set from the command line
echo '| '`basename "$0"`': loading settings from command-line...'
for var in "$@"
do
    export "$var"
done

# Save the settings this job was run with (for future defaults)
echo '| '`basename "$0"`': saving all loaded settings to settings.sh.'
output_settings > settings.sh


#if [ -z ${chosen_OMP_NUM_THREADS+x} ]
#then
#	echo "chosen_OMP_NUM_THREADS is unset"
#else
#	echo "chosen_OMP_NUM_THREADS is set to '$chosen_OMP_NUM_THREADS'"
#fi


# truth value of $runPythia evaluated inside
echo '| '`basename "$0"`': running Pythia!'
./run_Pythia.sh



# apply HBT analysis to each chosen event class
#for eventClassCutString in "0-100%" "0-0.1%" "0-1%" "0-5%" "5-10%" "0-10%" "10-20%" "20-30%" "30-40%" "40-50%" "50-60%" "60-70%" "70-80%" "80-90%" "90-100%"
for eventClassCutString in "1-11" "12-16" "17-22" "23-28" "29-34" "35-41" "42-51" "52-151"
do
	echo '| '`basename "$0"`': running HBT analysis for event class = '$eventClassCutString'!'
	./run_HBT_analysis.sh $eventClassCutString	# do NOT submit Bash scripts in background
done	# all centralities finished



#zipFilename=$CURRENT_RESULTS_DIRECTORY".zip"

#zip -r $zipFilename $CURRENT_RESULTS_DIRECTORY


echo '| '`basename "$0"`': Finished everything!'
echo '--------------------------------------------------------------------------------'	# width is 80 spaces


# End of file
