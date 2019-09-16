#! /usr/bin/env bash

source pbs_env.sh

########################################
# Fix OpenMP settings and compile
chosen_OMP_NUM_THREADS=$1
echo 'chosen_OMP_NUM_THREADS='$chosen_OMP_NUM_THREADS > omp_env.sh

./compile_all.sh			\
	$chosen_OMP_NUM_THREADS	\
	&> compile_all.out

########################################
# set up array of job specifications

declare -a specs=(	'useParallel=true projectile="Pb" target="Pb" beamEnergy="2760.0" Nevents=100000 ThermalOnly="true" SetPartonVertices="off"'
					'useParallel=true projectile="Pb" target="Pb" beamEnergy="2760.0" Nevents=100000 ThermalOnly="false" SetPartonVertices="off"'
					'useParallel=true projectile="Pb" target="Pb" beamEnergy="2760.0" Nevents=100000 ThermalOnly="true" SetPartonVertices="on"'
					'useParallel=true projectile="Pb" target="Pb" beamEnergy="2760.0" Nevents=100000 ThermalOnly="false" SetPartonVertices="on"'
				)

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
	cp -r src $HOME_RESULTS_DIRECTORY/job-${job}/
	
	executableString="./driver.sh ${specs[i]} &> driver.out"
	generate_pbs '18:00:00' $chosen_OMP_NUM_THREADS $executableString > $HOME_RESULTS_DIRECTORY/job-${job}/submit.pbs
	
	cp scripts/driver.sh $HOME_RESULTS_DIRECTORY/job-${job}
	cp defaults.sh env.sh omp_env.sh $HOME_RESULTS_DIRECTORY/job-${job}
done

# End of file