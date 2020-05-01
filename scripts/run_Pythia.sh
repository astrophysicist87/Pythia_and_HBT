#! /usr/bin/env bash
#-------------------

CWD=`pwd`
echo '| ------------------------------------------------------------'
echo '| - '`basename "$0"`': Executing this script in the following directory:'
echo '| | '$CWD

# convenient toggle
declare -A boolVal=( ["true"]="1" ["false"]="0")

# make sure main results directory exists
if [ ! -d "$MAIN_RESULTS_DIRECTORY" ]
then
	mkdir $MAIN_RESULTS_DIRECTORY
	echo '| - '`basename "$0"`': Created directory:' `realpath --relative-to="${PWD}" "$MAIN_RESULTS_DIRECTORY"`
fi
CURRENT_RESULTS_DIRECTORY=$MAIN_RESULTS_DIRECTORY


echo '| - '`basename "$0"`': Processing' \
		$Nevents $projectile'+'$target \
		'collisions at' $beamEnergy 'GeV'


clean_directory $PYTHIA_DIRECTORY
PYTHIA_RESULTS_DIRECTORY=$CURRENT_RESULTS_DIRECTORY/Pythia_results

collisionSystemStem=$projectile$target"_"`echo $beamEnergy`"GeV_Nev"$Nevents

#=================
# Run Pythia here
(
	echo '| - '`basename "$0"`':     Now entering '`realpath --relative-to="${PWD}" "$PYTHIA_DIRECTORY"`
	cd $PYTHIA_DIRECTORY

	cp $HOME_DIRECTORY/parameters.dat .

	#----------------------------------------------------
	# If false, no Pythia settings or code run, but still
	# export appropriate directories and file paths to
	# other directories for post-processing
	if $runPythia
	then

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
		./run_BEeffects_OpenMP.sh \
							$projectile $target \
							$PYTHIA_RESULTS_DIRECTORY \
							pythiaHBT::beam_energy=$beamEnergy \
							pythiaHBT::number_of_events=$Nevents \
							pythiaHBT::chosenParticleID=$chosenHBTparticle \
							pythiaHBT::lowerLimit=$lowerLimit \
							pythiaHBT::upperLimit=$upperLimit \
							pythiaHBT::bmin=$bMin \
							pythiaHBT::bmax=$bMax \
							pythiaHBT::output_Bjorken_variables="${boolVal[$storeBjorkenCoordinates]}"

		# check and report whether run was successful
		runSuccess=`echo $?`
		check_success 'Pythia' $runSuccess '-' '| | - '

	fi


	# Get the filenames which need to be processed
	recordOfOutputFilenames_Sxp=$PYTHIA_RESULTS_DIRECTORY/`echo $collisionSystemStem`"_S_x_p_filenames.dat"
	recordOfOutputFilename_mult=$PYTHIA_RESULTS_DIRECTORY/`echo $collisionSystemStem`"_total_N_filename.dat"
	\rm -f $HBT_EVENT_GEN_DIRECTORY/catalogue.dat
	\rm -f $HBT_SV_DIRECTORY/catalogue.dat
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
	#echo $projectile $target $beamEnergy $lowerLimit $upperLimit $Nevents > $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
	echo $projectile $target $beamEnergy $Nevents > $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
	readlink -f $PYTHIA_RESULTS_DIRECTORY/`cat $recordOfOutputFilename_mult` >> $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
	#echo $projectile $target $beamEnergy $lowerLimit $upperLimit $Nevents > $HBT_SV_DIRECTORY/ensemble_catalogue.dat
	echo $projectile $target $beamEnergy $Nevents > $HBT_SV_DIRECTORY/ensemble_catalogue.dat
	readlink -f $PYTHIA_RESULTS_DIRECTORY/`cat $recordOfOutputFilename_mult` >> $HBT_SV_DIRECTORY/ensemble_catalogue.dat

	cp $PYTHIA_DIRECTORY/main_BEeffects.cmnd $PYTHIA_RESULTS_DIRECTORY

)

echo '| - '`basename "$0"`': Finished everything!'
echo '| ------------------------------------------------------------'


# End of file
