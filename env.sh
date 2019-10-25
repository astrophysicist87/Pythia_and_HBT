#=====================================
# Header info
#=====================================
HOME_DIRECTORY=`pwd`
SOURCE_DIRECTORY=$HOME_DIRECTORY/src
# Pythia
versionNumber=8244
PYTHIA_DIRECTORY=$SOURCE_DIRECTORY/pythia${versionNumber}/examples
# HBT event generator
HBT_DIRECTORY=$SOURCE_DIRECTORY/HBT_event_generator
# HBT event generator
HBT_EVENT_GEN_DIRECTORY=$HBT_DIRECTORY/HBT_event_generator_w_errors
# Fit correlation function
HBT_FITCF_DIRECTORY=$HBT_DIRECTORY/fit_correlation_function
# Source variances/HBT radii
HBT_SV_DIRECTORY=$HBT_DIRECTORY/source_variances

# Main results directory
MAIN_RESULTS_DIRECTORY=$HOME_DIRECTORY/results
#CURRENT_RESULTS_DIRECTORY=$MAIN_RESULTS_DIRECTORY/results

# convenient boolean conversion table
declare -A boolVal=( ["true"]="1" ["false"]="0")

#============================
# Some function definitions
#============================
check_success () {
	if [ "$2" -eq "0" ]
	then
		echo '===================================================='
		echo '== '$1': run completed successfully!'
		echo '===================================================='
	else
		echo '===================================================='
		echo '== '$1': problems encountered!'
		echo '===================================================='
		exit $2
	fi
}

clean_directory () {
	rm $1/*.out $1/*.err $1/*.txt
	#rm $1/*catalogue.dat
	#rm $1/parameters.dat
	rm -rf $1/results
}


output_settings () {

echo 'runPythia='$runPythia
#echo 'useParallel='$useParallel
echo 'centralitySelectionInPythia='$centralitySelectionInPythia
echo 'runHBTEG='$runHBTEG
echo 'runFitCF='$runFitCF
echo 'runSV='$runSV
#echo 'runBF='$runBF
echo
echo 'projectile='$projectile
echo 'target='$target
echo 'beamEnergy='$beamEnergy
echo 'Nevents='$Nevents
echo 'chosenHBTparticle='$chosenHBTparticle
echo 'centralityClass='$centralityClass
echo 'bMin='$bMin
echo 'bMax='$bMax
echo
echo 'QRefValue='$QRefValue
echo 'BEeffects='$BEeffects
echo 'BEEnhancementMode='$BEEnhancementMode
echo 'SetFragmentationVertices='$SetFragmentationVertices
echo 'SetPartonVertices='$SetPartonVertices
echo 'ThermalOnly='$ThermalOnly
echo
echo 'UseColorReconnection='$UseColorReconnection
echo 'UseRopeHadronization='$UseRopeHadronization
echo 'IncludeStringShoving='$IncludeStringShoving
echo 'IncludeFlavourRopesMechanism='$IncludeFlavourRopesMechanism
echo 'StoreBjorkenCoordinates='$StoreBjorkenCoordinates

}


# End of file
