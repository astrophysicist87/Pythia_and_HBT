#! /usr/bin/env bash

# Load header info
source scripts/env.sh

# Load default variable values
source scripts/defaults.sh

# Load OpenMP info
source scripts/omp_env.sh
export OMP_NUM_THREADS=$chosen_OMP_NUM_THREADS
echo 'driver_pbs.sh: OMP_NUM_THREADS =' $OMP_NUM_THREADS

# Load PBS script defaults
source scripts/pbs_env.sh

# Load job specifications and event class ranges
source scripts/specs.sh

# Update any variables set from the command line
for var in "$@"
do
    export "$var"
done


# Save the settings this job was run with (for future defaults)
output_settings > settings.sh

# Total number of events = NDATASETS * Nevents
export NTotalEvents=$[NDATASETS*Nevents]

# N.B. - truth value of $runPythia evaluated inside run_Pythia.pbs

#------------------------
# Submit all Pythia datasets
echo "Submitting qsub -q $chosen_QUEUENAME -l walltime=$chosen_Pythia_walltime -l nodes=1:ppn=$OMP_NUM_THREADS -v datasetSeed=0 -V run_Pythia.pbs"
jobids=`qsub -q $chosen_QUEUENAME -l walltime=$chosen_Pythia_walltime -l nodes=1:ppn=$OMP_NUM_THREADS -v datasetSeed=0 -V run_Pythia.pbs`
echo '--------'
for iDataset in $(seq 1 $[NDATASETS-1])
do
	echo "Submitting qsub -q $chosen_QUEUENAME -l walltime=$chosen_Pythia_walltime -l nodes=1:ppn=$OMP_NUM_THREADS -v datasetSeed=$iDataset -V run_Pythia.pbs"
	jobid=`qsub -q $chosen_QUEUENAME -l walltime=$chosen_Pythia_walltime -l nodes=1:ppn=$OMP_NUM_THREADS -v datasetSeed=$iDataset -V run_Pythia.pbs`
	jobids+=":${jobid}"
	echo '--------'
done

echo
echo

#------------------------
# Wait until Pythia jobs finish, then post-process before running HBT analyses
echo "Submitting qsub -l walltime=00:15:00 -l nodes=1:ppn=$OMP_NUM_THREADS -W depend=afterok:${jobids} -V post_process_Pythia.pbs"
jobid=`qsub -q $chosen_QUEUENAME -l walltime=00:15:00 -l nodes=1:ppn=$OMP_NUM_THREADS -W depend=afterok:${jobids} -V post_process_Pythia.pbs`

echo 'jobid =' $jobid

echo
echo

usingFilemode1=false

if $usingFilemode
then

	#------------------------
	# apply HBT analysis to each chosen event class after 
	# post-processing of Pythia datasets is complete
	for eventClassCutString in "${class_ranges[@]}"
	do
		varString="eventClassCutString=${eventClassCutString},datasetSeed=0,HBTfilemode=1"
		echo "varString =" $varString
		echo "Submitting qsub -q $chosen_QUEUENAME -l walltime=$chosen_HBT_walltime_per_event_class -l nodes=1:ppn=$OMP_NUM_THREADS -v $varString -W depend=afterok:${jobid} run_HBT_analysis.pbs"
		qsub -q $chosen_QUEUENAME                            \
			-l walltime=$chosen_HBT_walltime_per_event_class \
			-l nodes=1:ppn=$OMP_NUM_THREADS                  \
			-v $varString                                    \
			-W depend=afterok:${jobid}                       \
			run_HBT_analysis.pbs
		echo '--------'
	done	# all event classes finished

else

	echo 'TESTING FILEMODE = 2'

	#------------------------
	# apply HBT analysis to each chosen event class after 
	# post-processing of Pythia datasets is complete
	for eventClassCutString in "${class_ranges[@]}"
	do
		varString="eventClassCutString=${eventClassCutString},datasetSeed=0,HBTfilemode=2"
		jobids=`qsub -q $chosen_QUEUENAME                        \
				-l walltime=$chosen_HBT_walltime_per_event_class \
				-l nodes=1:ppn=$OMP_NUM_THREADS                  \
				-v $varString                                    \
				-W depend=afterok:${jobid}                       \
				run_HBT_analysis.pbs`
		echo '--------'
		for iDataset in $(seq 1 $[NDATASETS-1])
		do
			varString="eventClassCutString=${eventClassCutString},datasetSeed=${iDataset},HBTfilemode=2"
			newjobid=`qsub -q $chosen_QUEUENAME                          \
						-l walltime=$chosen_HBT_walltime_per_event_class \
						-l nodes=1:ppn=$OMP_NUM_THREADS                  \
						-v $varString                                    \
						-W depend=afterok:${jobid}                       \
						run_HBT_analysis.pbs`
			jobids+=":${newjobid}"
			echo '--------'
		done
	done	# all event classes finished

fi


#zipFilename=$CURRENT_RESULTS_DIRECTORY".zip"

#zip -r $zipFilename $CURRENT_RESULTS_DIRECTORY

echo

echo 'driver_pbs.sh: Finished everything!'


# End of file
