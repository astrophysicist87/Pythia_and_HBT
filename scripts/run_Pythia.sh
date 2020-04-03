#! /usr/bin/env bash

#=============
# BEGIN SETUP
#=============
# Load header info
#source env.sh

# Load default variable values
#source defaults.sh

# Load OpenMP info
source omp_env.sh
#export OMP_NUM_THREADS=$chosen_OMP_NUM_THREADS
#echo 'OMP_NUM_THREADS =' $OMP_NUM_THREADS

# Update any variables set from the command line
#for var in "$@"
#do
#    export "$var"
#done

# Save the settings this job was run with (for future defaults)
#output_settings > settings.sh

#===========
# END SETUP
#===========

if [ -z ${chosen_OMP_NUM_THREADS+x} ]
then
	echo "chosen_OMP_NUM_THREADS is unset"
else
	echo "chosen_OMP_NUM_THREADS is set to '$chosen_OMP_NUM_THREADS'"
fi


# make sure main results directory exists
if [ ! -d "$MAIN_RESULTS_DIRECTORY" ]
then
	mkdir $MAIN_RESULTS_DIRECTORY
	echo 'Created' $MAIN_RESULTS_DIRECTORY
fi
CURRENT_RESULTS_DIRECTORY=$MAIN_RESULTS_DIRECTORY
if [ -z ${chosen_OMP_NUM_THREADS+x} ]
then
	echo "chosen_OMP_NUM_THREADS is unset"
else
	echo "chosen_OMP_NUM_THREADS is set to '$chosen_OMP_NUM_THREADS'"
fi


echo 'RUN_PYTHIA: Processing Nevents =' \
		$Nevents $projectile'+'$target \
		'collisions at' $beamEnergy 'GeV'

if [ -z ${chosen_OMP_NUM_THREADS+x} ]
then
	echo "chosen_OMP_NUM_THREADS is unset"
else
	echo "chosen_OMP_NUM_THREADS is set to '$chosen_OMP_NUM_THREADS'"
fi


echo 'Doing some checks inside:'
echo 'OMP_NUM_THREADS = '$OMP_NUM_THREADS
echo 'chosen_OMP_NUM_THREADS = '$chosen_OMP_NUM_THREADS
echo 'storeBjorkenCoordinates = '$storeBjorkenCoordinates
echo 'versionNumber = '$versionNumber
echo 'bMax = '$bMax


clean_directory $PYTHIA_DIRECTORY
PYTHIA_RESULTS_DIRECTORY=$CURRENT_RESULTS_DIRECTORY/Pythia_results

collisionSystemStem=$projectile$target"_"`echo $beamEnergy`"GeV_Nev"$Nevents

#=================
# Run Pythia here
(
	cd $PYTHIA_DIRECTORY
	echo '     Now in '`pwd`

	# make sure results directory exists
	if [ ! -d "$PYTHIA_RESULTS_DIRECTORY" ]
	then
		mkdir $PYTHIA_RESULTS_DIRECTORY
	fi

	# set some options for Pythia to run on
	rm main_BEeffects.cmnd 2>/dev/null
	# Turn on tracking of space-time information
	echo 'Fragmentation:setVertices =' $SetFragmentationVertices >> main_BEeffects.cmnd
	echo 'PartonVertex:setVertex =' $SetPartonVertices >> main_BEeffects.cmnd

	# turn on and set Bose-Einstein effects
	echo 'HadronLevel:BoseEinstein =' $BEeffects >> main_BEeffects.cmnd
	echo 'BoseEinstein:QRef =' $QRefValue >> main_BEeffects.cmnd

	# turn on/off other interesting mechanisms to test
	echo 'ColourReconnection:reconnect = ' $UseColorReconnection >> main_BEeffects.cmnd
	echo 'Ropewalk:RopeHadronization = ' $UseRopeHadronization >> main_BEeffects.cmnd
	echo 'Ropewalk:doShoving = ' $IncludeStringShoving >> main_BEeffects.cmnd
	echo 'Ropewalk:doFlavour = ' $IncludeFlavourRopesMechanism >> main_BEeffects.cmnd

	# New Pythia options, flags, parameters, etc.
	# which I've added myself (not generally compatible yet)
	# Comment out lines below this one if running on unmodified Pythia
	echo 'BoseEinstein:enhanceMode =' $BEEnhancementMode >> main_BEeffects.cmnd
	echo 'BoseEinstein:useInvariantSourceSize =' $useInvariantSourceSize >> main_BEeffects.cmnd
	echo 'BoseEinstein:useDistribution =' $useDistribution >> main_BEeffects.cmnd
	echo 'BoseEinstein:useRelativeDistance =' $useRelativeDistance >> main_BEeffects.cmnd
	echo 'BoseEinstein:useRestFrame =' $useRestFrame >> main_BEeffects.cmnd
	echo 'BoseEinstein:includePhaseSpace =' $includePhaseSpace >> main_BEeffects.cmnd
	echo 'BoseEinstein:linearInterpolateCDF =' $linearInterpolateCDF >> main_BEeffects.cmnd
	echo 'BoseEinstein:usePositiveShiftsForCompensation =' $usePositiveShiftsForCompensation >> main_BEeffects.cmnd
	echo 'BoseEinstein:computeBEEnhancementExactly =' $computeBEEnhancementExactly >> main_BEeffects.cmnd

	# Default is now to compute all events in Pythia from the get-go
	lowerLimit=0
	upperLimit=100

	# time and run
	./run_BEeffects_arbtryParticle_OpenMP.sh \
						$projectile $target $beamEnergy \
						$Nevents $chosenHBTparticle \
						$PYTHIA_RESULTS_DIRECTORY \
						$lowerLimit $upperLimit \
						$bMin $bMax $storeBjorkenCoordinates

	# check and report whether run was successful
	runSuccess=`echo $?`
	check_success 'Pythia' $runSuccess


	# if Pythia was minimum bias (default), do centrality selection in subsequent codes
	# otherwise, just do whatever events have been produced
	lowerLimit=${centralityCut[0]}
	upperLimit=${centralityCut[1]}
	if $centralitySelectionInPythia
	then
		lowerLimit=0
		upperLimit=100
	fi


	# Get the filenames which need to be processed
	recordOfOutputFilenames_Sxp=$PYTHIA_RESULTS_DIRECTORY/`echo $collisionSystemStem`"_S_x_p_filenames.dat"
	recordOfOutputFilename_mult=$PYTHIA_RESULTS_DIRECTORY/`echo $collisionSystemStem`"_total_N_filename.dat"
	\rm $HBT_EVENT_GEN_DIRECTORY/catalogue.dat
	\rm $HBT_SV_DIRECTORY/catalogue.dat
	for line in `cat $recordOfOutputFilenames_Sxp`
	do
		readlink -f $PYTHIA_RESULTS_DIRECTORY/$line >> $HBT_EVENT_GEN_DIRECTORY/catalogue.dat
		readlink -f $PYTHIA_RESULTS_DIRECTORY/$line >> $HBT_SV_DIRECTORY/catalogue.dat
	done
	# Set particle catalogue
	readlink -f $PYTHIA_RESULTS_DIRECTORY/HBT_particle_0.dat > $HBT_EVENT_GEN_DIRECTORY/particle_catalogue.dat
	readlink -f $PYTHIA_RESULTS_DIRECTORY/HBT_particle_0.dat > $HBT_FITCF_DIRECTORY/particle_catalogue.dat
	readlink -f $PYTHIA_RESULTS_DIRECTORY/HBT_particle_0.dat > $HBT_SV_DIRECTORY/particle_catalogue.dat
	# Set ensemble catalogue
	echo $projectile $target $beamEnergy $lowerLimit $upperLimit $Nevents > $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
	readlink -f $PYTHIA_RESULTS_DIRECTORY/`cat $recordOfOutputFilename_mult` >> $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
	echo $projectile $target $beamEnergy $lowerLimit $upperLimit $Nevents > $HBT_SV_DIRECTORY/ensemble_catalogue.dat
	readlink -f $PYTHIA_RESULTS_DIRECTORY/`cat $recordOfOutputFilename_mult` >> $HBT_SV_DIRECTORY/ensemble_catalogue.dat

	cp $PYTHIA_DIRECTORY/main_BEeffects.cmnd $PYTHIA_RESULTS_DIRECTORY

)

echo 'Finished everything!'


# End of file
