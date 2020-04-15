#! /usr/bin/env bash

#=====================================
# Header info
#=====================================
#export HOME_DIRECTORY=`pwd`
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export HOME_DIRECTORY=`readlink -f $DIR/..`
export SCRIPTS_DIRECTORY=$HOME_DIRECTORY/scripts
export SOURCE_DIRECTORY=$HOME_DIRECTORY/src
# Pythia
export versionNumber=8243
export PYTHIA_DIRECTORY=$SOURCE_DIRECTORY/pythia${versionNumber}/examples
# HBT event generator
export HBT_DIRECTORY=$SOURCE_DIRECTORY/HBT_event_generator
# HBT event generator
export HBT_EVENT_GEN_DIRECTORY=$HBT_DIRECTORY/HBT_event_generator_w_errors
# Fit correlation function
export HBT_FITCF_DIRECTORY=$HBT_DIRECTORY/fit_correlation_function
# Source variances/HBT radii
export HBT_SV_DIRECTORY=$HBT_DIRECTORY/source_variances

# Main results directory
export MAIN_RESULTS_DIRECTORY=$HOME_DIRECTORY/results

# convenient boolean conversion table
declare -A boolVal=( ["true"]="1" ["false"]="0")

#============================
# Some function definitions
#============================
check_success () {
	cdef="="
	pdef=''
	c="${3:-$cdef}"
	prefix="${4:-$pdef}"
	if [ "$2" -eq "0" ]
	then
		#echo '===================================================='
		#echo '== '$1': run completed successfully!'
		#echo '===================================================='
		echo $prefix`printf '%0.s'$c {1..60}`
		string=`printf '%0.s'$c {1..2} && echo ' '$1': run completed successfully!'`
		echo $prefix$string
		echo $prefix`printf '%0.s'$c {1..60}`
	else
		#echo '===================================================='
		#echo '== '$1': problems encountered!'
		#echo '===================================================='
		echo $prefix`printf '%0.s'$c {1..60}`
		string=`printf '%0.s'$c {1..2} && echo ' '$1': problems encountered!'`
		echo $prefix$string
		echo $prefix`printf '%0.s'$c {1..60}`
		exit $2
	fi
}

clean_directory () {
	rm -f $1/*.out $1/*.err $1/*.txt
	#rm $1/*catalogue.dat
	#rm $1/parameters.dat
	rm -rf $1/results
}


output_settings () {

echo '#! /usr/bin/env bash'
echo
echo 'export runPythia='$runPythia
#echo 'export useParallel='$useParallel
echo 'export eventClassSelectionInPythia='$eventClassSelectionInPythia
echo 'export runHBTEG='$runHBTEG
echo 'export runFitCF='$runFitCF
echo 'export runSV='$runSV
#echo 'export runBF='$runBF
echo
echo 'export projectile='$projectile
echo 'export target='$target
echo 'export beamEnergy='$beamEnergy
echo 'export Nevents='$Nevents
echo 'export chosenHBTparticle='$chosenHBTparticle
echo 'export eventClassSelectionMode='$eventClassSelectionMode
echo 'export eventClass='$eventClass
echo 'export bMin='$bMin
echo 'export bMax='$bMax
echo
echo 'export QRefValue='$QRefValue
echo 'export BEeffects='$BEeffects
echo 'export BEEnhancementMode='$BEEnhancementMode
echo 'export SetFragmentationVertices='$SetFragmentationVertices
echo 'export SetPartonVertices='$SetPartonVertices
echo 'export ThermalOnly='$ThermalOnly
echo
echo 'export useInvariantSourceSize='$useInvariantSourceSize
echo 'export useDistribution='$useDistribution
echo 'export useRelativeDistance='$useRelativeDistance
echo 'export useRestFrame='$useRestFrame
echo 'export includePhaseSpace='$includePhaseSpace
echo 'export linearInterpolateCDF='$linearInterpolateCDF
echo 'export usePositiveShiftsForCompensation='$usePositiveShiftsForCompensation
echo 'export computeBEEnhancementExactly='$computeBEEnhancementExactly
echo
echo 'export UseColorReconnection='$UseColorReconnection
echo 'export UseRopeHadronization='$UseRopeHadronization
echo 'export IncludeStringShoving='$IncludeStringShoving
echo 'export IncludeFlavourRopesMechanism='$IncludeFlavourRopesMechanism
echo 'export storeBjorkenCoordinates='$storeBjorkenCoordinates

}

export -f check_success
export -f clean_directory
export -f output_settings

# End of file
