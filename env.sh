#=====================================
# Header info
#=====================================
HOME_DIRECTORY=`pwd`
SOURCE_DIRECTORY=$HOME_DIRECTORY/src
# Pythia
PYTHIA_DIRECTORY=$SOURCE_DIRECTORY/pythia8235/examples
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
	rm $1/*catalogue.dat
	rm $1/parameters.dat
	rm -rf $1/results
}


# End of file
