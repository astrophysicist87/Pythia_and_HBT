#! /usr/bin/env bash

# Load header info
source scripts/env.sh

# Load default variable values
source scripts/defaults.sh

# Load OpenMP info
source scripts/omp_env.sh
export OMP_NUM_THREADS=$chosen_OMP_NUM_THREADS
echo 'OMP_NUM_THREADS =' $OMP_NUM_THREADS

# Load PBS script defaults
source scripts/pbs_env.sh

# Update any variables set from the command line
for var in "$@"
do
    export "$var"
done


# Save the settings this job was run with (for future defaults)
output_settings > settings.sh


# truth value of $runPythia evaluated inside
./run_Pythia.sh $seed


# apply HBT analysis to each chosen event class
#for eventClassCutString in "1-11" "12-16" "17-22" "23-28" "29-34" "35-41" "42-51" "52-151" "152-1000000"
#for eventClassCutString in "0-20%" "20-40%" "40-60%" "60-90%"
for eventClassCutString in "0-100%"
do
	echo "Submitting qsub -l walltime=$chosen_HBT_walltime_per_event_class -l nodes=1:ppn=$OMP_NUM_THREADS run_HBT_analysis.pbs"
	qsub -l walltime=$chosen_HBT_walltime_per_event_class \
		-l nodes=1:ppn=$OMP_NUM_THREADS \
		-v "eventClassCutString=$eventClassCutString" \
		run_HBT_analysis.pbs
done	# all event classes finished


#zipFilename=$CURRENT_RESULTS_DIRECTORY".zip"

#zip -r $zipFilename $CURRENT_RESULTS_DIRECTORY


echo 'Finished everything!'


# End of file
