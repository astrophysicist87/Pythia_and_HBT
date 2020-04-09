#! /usr/bin/env bash

source pbs_env.sh

# Set PBS walltime directives
export chosen_Pythia_walltime=$3
export chosen_HBT_walltime_per_centrality=$4

########################################
# Fix OpenMP settings and compile
chosen_OMP_NUM_THREADS=$1
echo 'export chosen_OMP_NUM_THREADS='$chosen_OMP_NUM_THREADS > omp_env.sh

../compile_all.sh	\
	$chosen_OMP_NUM_THREADS	\
	&> ../compile_all.out

########################################
# set up array of job specifications

#declare -a specs=(
#    'useArbitraryParticle=true projectile="p" target="p" beamEnergy="7000.0" chosenHBTparticle="211" Nevents=10000000 storeBjorkenCoordinates="false" BEeffects="off"'
#    'useArbitraryParticle=true projectile="p" target="Pb" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=1000000 storeBjorkenCoordinates="false" BEeffects="off"'
#	'useArbitraryParticle=true projectile="Pb" target="Pb" beamEnergy="2760.0" chosenHBTparticle="211" Nevents=100000 storeBjorkenCoordinates="false" BEeffects="off"'
#	)
source specs.sh

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
	
	executableString="./driver_pbs.sh ${specs[i]} &> driver_pbs.out"
	generate_pbs $chosen_Pythia_walltime $chosen_OMP_NUM_THREADS $executableString > $HOME_RESULTS_DIRECTORY/job-${job}/submit.pbs
	
	cp driver_pbs.sh		$HOME_RESULTS_DIRECTORY/job-${job}
	cp run_Pythia.sh		$HOME_RESULTS_DIRECTORY/job-${job}
	cp run_HBT_analysis.pbs	$HOME_RESULTS_DIRECTORY/job-${job}
	cp rerun_pbs.sh			$HOME_RESULTS_DIRECTORY/job-${job}
	cp defaults.sh specs.sh \
		env.sh omp_env.sh pbs_env.sh \
							$HOME_RESULTS_DIRECTORY/job-${job}
done

# End of file
