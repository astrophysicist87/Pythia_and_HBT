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

declare -a specs=(
	'useArbitraryParticle=true projectile="C" target="C" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=0.0 bMax=1.0'
	'useArbitraryParticle=true projectile="O" target="O" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=0.0 bMax=1.0'
	'useArbitraryParticle=true projectile="Cu" target="Cu" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=0.0 bMax=1.0'
	'useArbitraryParticle=true projectile="Xe" target="Xe" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=0.0 bMax=1.0'
	'useArbitraryParticle=true projectile="Au" target="Au" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=0.0 bMax=1.0'
	'useArbitraryParticle=true projectile="Pb" target="Pb" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=0.0 bMax=1.0'
	'useArbitraryParticle=true projectile="U" target="U" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=0.0 bMax=1.0'

	'useArbitraryParticle=true projectile="C" target="C" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=1.0 bMax=5.0'
	'useArbitraryParticle=true projectile="O" target="O" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=1.0 bMax=5.0'
	'useArbitraryParticle=true projectile="Cu" target="Cu" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=1.0 bMax=5.0'
	'useArbitraryParticle=true projectile="Xe" target="Xe" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=1.0 bMax=5.0'
	'useArbitraryParticle=true projectile="Au" target="Au" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=1.0 bMax=5.0'
	'useArbitraryParticle=true projectile="Pb" target="Pb" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=1.0 bMax=5.0'
	'useArbitraryParticle=true projectile="U" target="U" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=1.0 bMax=5.0'

	'useArbitraryParticle=true projectile="C" target="C" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=5.0 bMax=20.0'
	'useArbitraryParticle=true projectile="O" target="O" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=5.0 bMax=20.0'
	'useArbitraryParticle=true projectile="Cu" target="Cu" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=5.0 bMax=20.0'
	'useArbitraryParticle=true projectile="Xe" target="Xe" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=5.0 bMax=20.0'
	'useArbitraryParticle=true projectile="Au" target="Au" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=5.0 bMax=20.0'
	'useArbitraryParticle=true projectile="Pb" target="Pb" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=5.0 bMax=20.0'
	'useArbitraryParticle=true projectile="U" target="U" beamEnergy="5020.0" chosenHBTparticle="211" Nevents=10 runHBTEG=false runFitCF=false bMin=5.0 bMax=20.0'
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
	generate_pbs '00:15:00' $chosen_OMP_NUM_THREADS $executableString > $HOME_RESULTS_DIRECTORY/job-${job}/submit.pbs
	
	cp scripts/driver.sh $HOME_RESULTS_DIRECTORY/job-${job}
	cp scripts/rerun.sh $HOME_RESULTS_DIRECTORY/job-${job}
	cp defaults.sh env.sh omp_env.sh $HOME_RESULTS_DIRECTORY/job-${job}
done

# End of file
