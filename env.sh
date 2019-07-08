#=====================================
# Header info
#=====================================
HOME_DIRECTORY=`pwd`
# Pythia
PYTHIA_DIRECTORY=$HOME_DIRECTORY/pythia8235/examples
# HBT event generator
HBT_EVENT_GEN_DIRECTORY=$HOME_DIRECTORY/HBT_event_generator/HBT_event_generator_w_errors
# Fit correlation function
HBT_FITCF_DIRECTORY=$HOME_DIRECTORY/HBT_event_generator/fit_correlation_function
# Source variances/HBT radii
HBT_SV_DIRECTORY=$HOME_DIRECTORY/HBT_event_generator/source_variances

# Main results directory
MAIN_RESULTS_DIRECTORY=$HOME_DIRECTORY/RESULTS
CURRENT_RESULTS_DIRECTORY=$MAIN_RESULTS_DIRECTORY/results

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
