#! /usr/bin/env bash
#-------------------

CWD=`pwd`
echo '| ------------------------------------------------------------'
echo '| - '`basename "$0"`': Executing this script in the following directory:'
echo '| | '$CWD

declare -A boolVal=( ["true"]="1" ["false"]="0")

# make sure main results directory exists
if [ ! -d "$MAIN_RESULTS_DIRECTORY" ]
then
	mkdir $MAIN_RESULTS_DIRECTORY
	echo '| - '`basename "$0"`': Created directory:' `realpath --relative-to="${PWD}" "$MAIN_RESULTS_DIRECTORY"`
fi


#=================================================
# Miscellaneous set-up
CURRENT_RESULTS_DIRECTORY=$MAIN_RESULTS_DIRECTORY

#=====================================
# Clean all working directories
#clean_directory $HBT_DIRECTORY
#clean_directory $HBT_EVENT_GEN_DIRECTORY
#clean_directory $HBT_FITCF_DIRECTORY
#clean_directory $HBT_SV_DIRECTORY
#=====================================

#=================================================
# Set event class from command-line argument
#=================================================
eventClassCutString=$1	# should be passed in as argument directly

#=================================================
# process event class information
success=0

declare -A selectionTokens=( ["centrality"]="C" ["multiplicity"]="N")
declare -A selectionModes=( ["centrality"]="0" ["multiplicity"]="1")
selectionToken="${selectionTokens[$eventClassSelectionMode]}"
selectionMode="${selectionModes[$eventClassSelectionMode]}"

eventClassCut=(`echo $eventClassCutString | sed 's/-/ /g' | sed 's/%//g'`)
thisEventClass=${selectionToken}${eventClassCut[0]}"_"${eventClassCut[1]}


#=================================================
# if Pythia was minimum bias (default), do event class selection in subsequent codes
# otherwise, just do whatever events have been produced
lowerLimit=${eventClassCut[0]}
upperLimit=${eventClassCut[1]}
if $eventClassSelectionInPythia
then
	lowerLimit=0
	upperLimit=100
fi


#=================================================
# set names of sub-directories
HBT_RESULTS_DIRECTORY=$CURRENT_RESULTS_DIRECTORY/HBT_results
HBT_CEN_RESULTS_DIRECTORY=$HBT_RESULTS_DIRECTORY/$thisEventClass

# make sure HBT results directory exists
if [ -d "$HBT_CEN_RESULTS_DIRECTORY" ]
then
	rm -rf $HBT_CEN_RESULTS_DIRECTORY
fi
mkdir -p $HBT_CEN_RESULTS_DIRECTORY
echo '| - '`basename "$0"`': Created directory:' `realpath --relative-to="${PWD}" "$HBT_CEN_RESULTS_DIRECTORY"`
#=================================================


collisionSystemStem=$projectile$target"_"`echo $beamEnergy`"GeV_Nev"$Nevents


#=====================================
#=====================================
#=====================================

#===================
# Main calculation
#===================

echo '| - '`basename "$0"`': Processing' \
		$Nevents $projectile'+'$target \
		'collisions at' $beamEnergy 'GeV' \
		'in event class' $eventClassCutString \
		'(selection mode = '$eventClassSelectionMode')'


# Run HBT_event_generator
if $runHBTEG
then
(

	echo '| - '`basename "$0"`': Now entering '`realpath --relative-to="${PWD}" "$HBT_EVENT_GEN_DIRECTORY"`
	cd $HBT_EVENT_GEN_DIRECTORY

	# using OpenMP (leave a couple cores free)
	export OMP_NUM_THREADS=$chosen_OMP_NUM_THREADS

	CF_RESULTS_DIRECTORY=$HBT_CEN_RESULTS_DIRECTORY/CF_results

	if [ ! -d "$CF_RESULTS_DIRECTORY" ]; then
		mkdir -p $CF_RESULTS_DIRECTORY
		mkdir -p $CF_RESULTS_DIRECTORY/state
		#echo '| - '`basename "$0"`': Created directory "results" in' `realpath --relative-to="${PWD}" "$HBT_EVENT_GEN_DIRECTORY"`
		echo '| - '`basename "$0"`': Created directory "CF_results" in' `realpath --relative-to="${PWD}" "$HBT_CEN_RESULTS_DIRECTORY"`
	fi

	#cp ../parameters.dat .
	cp ../parameters.dat ./*catalogue.dat $CF_RESULTS_DIRECTORY

	#default: assume BE effects are turned off in Pythia
	chosen_BE_mode=0
	if [ "$BEeffects" == "on" ]; then
			chosen_BE_mode=1
	fi

	# time and run
	nohup time nice -n 10 ./run_HBT_event_generator.e \
			$CF_RESULTS_DIRECTORY \
			$CF_RESULTS_DIRECTORY/parameters.dat \
			$CF_RESULTS_DIRECTORY/particle_catalogue.dat \
			$CF_RESULTS_DIRECTORY/catalogue.dat \
			$CF_RESULTS_DIRECTORY/ensemble_catalogue.dat \
			selection_mode=$selectionMode \
			event_class_minimum=$lowerLimit \
			event_class_maximum=$upperLimit \
			BE_mode=$chosen_BE_mode \
			chosen_MCID=$chosenHBTparticle \
			store_Bjorken_coordinates="${boolVal[$storeBjorkenCoordinates]}" \
			1> $CF_RESULTS_DIRECTORY/HBT_event_generator.out \
			2> $CF_RESULTS_DIRECTORY/HBT_event_generator.err
	# N.B. - centralities determined either here or in Pythia, but not both

	# check and report whether run was successful
	runSuccess=`echo $?`
	check_success 'HBT_event_generator' $runSuccess '-' '| | - '

	# copy results
	#cp HBT_event_generator.[oe]* ./results
	#mkdir $HBT_CEN_RESULTS_DIRECTORY/CF_results
	#cp -r ./results/* $HBT_CEN_RESULTS_DIRECTORY/CF_results

	readlink -f $CF_RESULTS_DIRECTORY/HBT_pipiCF.dat > $HBT_FITCF_DIRECTORY/catalogue.dat

	#exit $runSuccess
)
fi

#=====================================
#=====================================
#=====================================
# Run fit_correlation_function
if $runFitCF
then
(

	echo '| - '`basename "$0"`': Now entering '`realpath --relative-to="${PWD}" "$HBT_FITCF_DIRECTORY"`
	cd $HBT_FITCF_DIRECTORY

	FIT_RESULTS_DIRECTORY=$HBT_CEN_RESULTS_DIRECTORY/fit_results

	if [ ! -d "$FIT_RESULTS_DIRECTORY" ]
	then
		mkdir $FIT_RESULTS_DIRECTORY
		#echo '| - '`basename "$0"`': Created directory "results" in' `realpath --relative-to="${PWD}" "$HBT_FITCF_DIRECTORY"`
		echo '| - '`basename "$0"`': Created directory "fit_results" in' `realpath --relative-to="${PWD}" "$HBT_CEN_RESULTS_DIRECTORY"`
	fi

	#cp ../parameters.dat .
	cp ../parameters.dat ./*catalogue.dat $FIT_RESULTS_DIRECTORY

	# time and run
	nohup time ./run_fit_correlation_function.e \
			$FIT_RESULTS_DIRECTORY \
			$FIT_RESULTS_DIRECTORY/parameters.dat \
			$FIT_RESULTS_DIRECTORY/particle_catalogue.dat \
			$FIT_RESULTS_DIRECTORY/catalogue.dat \
			1> $FIT_RESULTS_DIRECTORY/fit_correlation_function.out \
			2> $FIT_RESULTS_DIRECTORY/fit_correlation_function.err

	# check and report whether run was successful
	runSuccess=`echo $?`
	check_success 'fit_correlation_function' $runSuccess '-' '| | - '

	# copy results
	#cp fit_correlation_function.[oe]* ./results
	#mkdir $HBT_CEN_RESULTS_DIRECTORY/fit_results
	#cp -r ./results/* $HBT_CEN_RESULTS_DIRECTORY/fit_results

	#exit $runSuccess
)
fi

#=====================================
#=====================================
#=====================================
# Run SV.e
if $runSV
then
(

	echo '| - '`basename "$0"`': Now entering '`realpath --relative-to="${PWD}" "$HBT_SV_DIRECTORY"`
	cd $HBT_SV_DIRECTORY

	SV_RESULTS_DIRECTORY=$HBT_CEN_RESULTS_DIRECTORY/SV_results

	if [ ! -d "$SV_RESULTS_DIRECTORY" ]
	then
		mkdir $SV_RESULTS_DIRECTORY
		#echo '| - '`basename "$0"`': Created directory "results" in' `realpath --relative-to="${PWD}" "$HBT_SV_DIRECTORY"`
		echo '| - '`basename "$0"`': Created directory "SV_results" in' `realpath --relative-to="${PWD}" "$HBT_CEN_RESULTS_DIRECTORY"`
	fi

	#cp ../parameters.dat .
	cp ../parameters.dat ./*catalogue.dat $SV_RESULTS_DIRECTORY

	# time and run
	nohup time nice -n 10 ./SV.e \
			$SV_RESULTS_DIRECTORY \
			$SV_RESULTS_DIRECTORY/parameters.dat \
			$SV_RESULTS_DIRECTORY/particle_catalogue.dat \
			$SV_RESULTS_DIRECTORY/catalogue.dat \
			$SV_RESULTS_DIRECTORY/ensemble_catalogue.dat \
			selection_mode=$selectionMode \
			event_class_minimum=$lowerLimit \
			event_class_maximum=$upperLimit \
			chosen_MCID=$chosenHBTparticle \
			store_Bjorken_coordinates="${boolVal[$storeBjorkenCoordinates]}" \
			1> $SV_RESULTS_DIRECTORY/SV_record.out \
			2> $SV_RESULTS_DIRECTORY/SV_record.err

	# check and report whether run was successful
	runSuccess=`echo $?`
	check_success 'source_variances' $runSuccess '-' '| | - '

	# copy results
	#cp SV_record.[oe]* ./results
	#mkdir $HBT_CEN_RESULTS_DIRECTORY/SV_results
	#cp -r ./results/* $HBT_CEN_RESULTS_DIRECTORY/SV_results

	#exit $runSuccess
)
fi

# Overkill, but copy it over anyway
cp $HBT_DIRECTORY/parameters.dat $HBT_CEN_RESULTS_DIRECTORY

echo '| - '`basename "$0"`': Finished everything for event class =' $eventClassCutString'!'
echo '| ------------------------------------------------------------'


# End of file
