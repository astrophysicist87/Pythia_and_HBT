#! /usr/bin/env bash

# Load header info
#source env.sh

# Load default variable values
#source defaults.sh

# Load OpenMP info
#source omp_env.sh
#export OMP_NUM_THREADS=$chosen_OMP_NUM_THREADS
#echo 'OMP_NUM_THREADS =' $OMP_NUM_THREADS

# Update any variables set from the command line
#for var in "$@"
#do
#    export "$var"
#done

declare -A boolVal=( ["true"]="1" ["false"]="0")

# make sure main results directory exists
if [ ! -d "$MAIN_RESULTS_DIRECTORY" ]
then
	mkdir $MAIN_RESULTS_DIRECTORY
	echo 'Created' $MAIN_RESULTS_DIRECTORY
fi

echo 'Doing some checks here:'
echo 'OMP_NUM_THREADS = '$OMP_NUM_THREADS
echo 'chosen_OMP_NUM_THREADS = '$chosen_OMP_NUM_THREADS
echo 'storeBjorkenCoordinates = '$storeBjorkenCoordinates

#===================
# Main calculation
#===================

echo 'Processing Nevents =' $Nevents $projectile'+'$target 'collisions at' $beamEnergy 'GeV'


# make sure current results directory exists
#nextCurrentResultsDirectoryName=`get_dirname $CURRENT_RESULTS_DIRECTORY '-' "true"`
#CURRENT_RESULTS_DIRECTORY=$nextCurrentResultsDirectoryName
#mkdir $CURRENT_RESULTS_DIRECTORY
#echo 'Created' $CURRENT_RESULTS_DIRECTORY
CURRENT_RESULTS_DIRECTORY=$MAIN_RESULTS_DIRECTORY

#=====================================
# Clean all working directories
clean_directory $HBT_DIRECTORY
clean_directory $HBT_EVENT_GEN_DIRECTORY
clean_directory $HBT_FITCF_DIRECTORY
clean_directory $HBT_SV_DIRECTORY
#=====================================


#=================================================
# Set centrality class from command-line argument
#=================================================
centralityCutString=$1	# should be passed in as argument directly


#=================================================
# process centrality class information
success=0
echo '  -- analyzing centrality class' $centralityCutString

centralityCut=(`echo $centralityCutString | sed 's/-/ /g' | sed 's/%//g'`)
thisCentrality="C"${centralityCut[0]}"_"${centralityCut[1]}


#=================================================
# if Pythia was minimum bias (default), do centrality selection in subsequent codes
# otherwise, just do whatever events have been produced
lowerLimit=${centralityCut[0]}
upperLimit=${centralityCut[1]}
if $centralitySelectionInPythia
then
	lowerLimit=0
	upperLimit=100
fi


#=================================================
# set names of sub-directories
HBT_RESULTS_DIRECTORY=$CURRENT_RESULTS_DIRECTORY/HBT_results
HBT_CEN_RESULTS_DIRECTORY=$HBT_RESULTS_DIRECTORY/$thisCentrality

# make sure HBT results directory exists
if [ -d "$HBT_CEN_RESULTS_DIRECTORY" ]
then
	rm -rf $HBT_CEN_RESULTS_DIRECTORY
fi
mkdir -p $HBT_CEN_RESULTS_DIRECTORY
#=================================================

collisionSystemStem=$projectile$target"_"`echo $beamEnergy`"GeV_Nev"$Nevents

#=====================================
#=====================================
#=====================================
# Run HBT_event_generator
if $runHBTEG
then
(

	cd $HBT_EVENT_GEN_DIRECTORY
	echo '     Now in '`pwd`

	# using OpenMP (leave a couple cores free)
	export OMP_NUM_THREADS=$chosen_OMP_NUM_THREADS

	if [ ! -d "./results" ]; then
		    mkdir results
	fi

	cp ../parameters.dat .

	#default: assume BE effects are turned off in Pythia
	chosen_BE_mode=0
	if [ "$BEeffects" == "on" ]; then
			chosen_BE_mode=1
	fi

	# time and run
	nohup time ./run_HBT_event_generator.e \
			BE_mode=$chosen_BE_mode \
			chosen_MCID=$chosenHBTparticle \
			store_Bjorken_coordinates="${boolVal[$storeBjorkenCoordinates]}" \
			1> HBT_event_generator.out \
			2> HBT_event_generator.err
	# N.B. - centralities determined either here or in Pythia, but not both

	# check and report whether run was successful
	runSuccess=`echo $?`
	check_success 'HBT_event_generator' $runSuccess

	# copy results
	cp HBT_event_generator.[oe]* ./results
	mkdir $HBT_CEN_RESULTS_DIRECTORY/CF_results
	cp -r ./results/* $HBT_CEN_RESULTS_DIRECTORY/CF_results

	readlink -f ./results/HBT_pipiCF.dat > $HBT_FITCF_DIRECTORY/catalogue.dat

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

	cd $HBT_FITCF_DIRECTORY
	echo '     Now in '`pwd`

	if [ ! -d "./results" ]
	then
		mkdir results
	fi

	cp ../parameters.dat .

	# time and run
	nohup time ./run_fit_correlation_function.e \
			1> fit_correlation_function.out \
			2> fit_correlation_function.err

	# check and report whether run was successful
	runSuccess=`echo $?`
	check_success 'fit_correlation_function' $runSuccess

	# copy results
	cp fit_correlation_function.[oe]* ./results
	mkdir $HBT_CEN_RESULTS_DIRECTORY/fit_results
	cp -r ./results/* $HBT_CEN_RESULTS_DIRECTORY/fit_results

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

	cd $HBT_SV_DIRECTORY
	echo '     Now in '`pwd`

	if [ ! -d "./results" ]
	then
		mkdir results
	fi

	cp ../parameters.dat .

	# time and run
	nohup time ./SV.e \
			chosen_MCID=$chosenHBTparticle \
			store_Bjorken_coordinates="${boolVal[$storeBjorkenCoordinates]}" \
			1> SV_record.out \
			2> SV_record.err

	# check and report whether run was successful
	runSuccess=`echo $?`
	check_success 'source_variances' $runSuccess

	# copy results
	cp SV_record.[oe]* ./results
	mkdir $HBT_CEN_RESULTS_DIRECTORY/SV_results
	cp -r ./results/* $HBT_CEN_RESULTS_DIRECTORY/SV_results

	#exit $runSuccess
)
fi

cp $HBT_DIRECTORY/parameters.dat $HBT_CEN_RESULTS_DIRECTORY

echo 'Finished everything!'


# End of file