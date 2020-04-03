#! /usr/bin/env bash

# Load header info
source env.sh

# Load default variable values
source defaults.sh

# Load OpenMP info
source omp_env.sh
export OMP_NUM_THREADS=$chosen_OMP_NUM_THREADS
echo 'OMP_NUM_THREADS =' $OMP_NUM_THREADS

# Update any variables set from the command line
for var in "$@"
do
    export "$var"
done



./run_Pythia.sh



for centralityCutString in "0-100%" "0-10%" "10-20%" "20-40%" "40-60%" "60-100%"
do
	qsub -l walltime=24:00:00 -l nodes=1 -l ppn=40 -v "centralityCutString=$centralityCutString" run_HBT_anlysis.pbs
done	# all centralities finished


#zipFilename=$CURRENT_RESULTS_DIRECTORY".zip"

#zip -r $zipFilename $CURRENT_RESULTS_DIRECTORY


echo 'Finished everything!'


# End of file