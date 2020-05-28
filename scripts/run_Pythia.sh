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

# only command-line argument is dataset seed
datasetSeed=$1

echo '| - '`basename "$0"`': Processing' \
		$Nevents $projectile'+'$target \
		'collisions at' $beamEnergy 'GeV'


clean_directory $PYTHIA_DIRECTORY
datasetIndex="_${datasetSeed}"
#if [ "$datasetSeed" -lt 0 ]; then
#	datasetIndex=""
#fi
#PYTHIA_RESULTS_DIRECTORY=$CURRENT_RESULTS_DIRECTORY/Pythia_results
#PYTHIA_RESULTS_DIRECTORY=$CURRENT_RESULTS_DIRECTORY/Pythia_results/dataset_${datasetSeed}
PYTHIA_RESULTS_DIRECTORY=$CURRENT_RESULTS_DIRECTORY/Pythia_results/dataset${datasetIndex}

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
			mkdir -p $PYTHIA_RESULTS_DIRECTORY
		fi

		#--------------------------
		# set some options for Pythia to run on
		PYTHIA_CMND_FILE=$PYTHIA_RESULTS_DIRECTORY/main_BEeffects.cmnd
		rm $PYTHIA_CMND_FILE 2>/dev/null
		
		# Set random seed
		#echo 'Random:setSeed = on'                                                      >> $PYTHIA_CMND_FILE
		#echo 'Random:seed ='                              $[datasetSeed+1]              >> $PYTHIA_CMND_FILE
		#echo 'Random:seed = 0'                                                          >> $PYTHIA_CMND_FILE

		# Turn on tracking of space-time information
		echo 'Fragmentation:setVertices ='                $SetFragmentationVertices     >> $PYTHIA_CMND_FILE
		echo 'PartonVertex:setVertex ='                   $SetPartonVertices            >> $PYTHIA_CMND_FILE

		# turn on and set Bose-Einstein effects
		echo 'HadronLevel:BoseEinstein ='                 $BEeffects                    >> $PYTHIA_CMND_FILE
		echo 'BoseEinstein:QRef ='                        $QRefValue                    >> $PYTHIA_CMND_FILE

		# turn on/off other interesting mechanisms to test
		echo 'ColourReconnection:reconnect ='             $UseColorReconnection         >> $PYTHIA_CMND_FILE
		echo 'Ropewalk:RopeHadronization ='               $UseRopeHadronization         >> $PYTHIA_CMND_FILE
		echo 'Ropewalk:doShoving ='                       $IncludeStringShoving         >> $PYTHIA_CMND_FILE
		echo 'Ropewalk:doFlavour ='                       $IncludeFlavourRopesMechanism >> $PYTHIA_CMND_FILE

		# New Pythia options, flags, parameters, etc.
		# which I've added myself (not generally compatible yet)
		# Comment out lines below this one if running on unmodified Pythia
		echo 'BoseEinstein:enhanceMode ='                 $BEEnhancementMode            >> $PYTHIA_CMND_FILE
		echo 'BoseEinstein:useInvariantSourceSize ='      $useInvariantSourceSize       >> $PYTHIA_CMND_FILE
		echo 'BoseEinstein:useDistribution ='             $useDistribution              >> $PYTHIA_CMND_FILE
		echo 'BoseEinstein:useRelativeDistance ='         $useRelativeDistance          >> $PYTHIA_CMND_FILE
		echo 'BoseEinstein:useRestFrame ='                $useRestFrame                 >> $PYTHIA_CMND_FILE
		echo 'BoseEinstein:includePhaseSpace ='           $includePhaseSpace            >> $PYTHIA_CMND_FILE
		echo 'BoseEinstein:linearInterpolateCDF ='        $linearInterpolateCDF         >> $PYTHIA_CMND_FILE
		echo 'BoseEinstein:computeBEEnhancementExactly =' $computeBEEnhancementExactly  >> $PYTHIA_CMND_FILE
		echo 'BoseEinstein:shiftingSet ='                 $shiftingSet                  >> $PYTHIA_CMND_FILE
		echo 'BoseEinstein:compensationSet ='             $compensationSet              >> $PYTHIA_CMND_FILE
		echo 'BoseEinstein:compensationMode ='            $compensationMode             >> $PYTHIA_CMND_FILE

		# Default is now to compute all events in Pythia from the get-go
		lowerLimit=0
		upperLimit=100

		# time and run
		nice -n -10 ./run_BEeffects_OpenMP.sh \
			$projectile $target \
			$PYTHIA_RESULTS_DIRECTORY \
			pythiaHBT::beam_energy=$beamEnergy \
			pythiaHBT::number_of_events=$Nevents \
			pythiaHBT::chosenParticleID=$chosenHBTparticle \
			pythiaHBT::lowerLimit=$lowerLimit \
			pythiaHBT::upperLimit=$upperLimit \
			pythiaHBT::bmin=$bMin \
			pythiaHBT::bmax=$bMax \
			pythiaHBT::seed=$datasetSeed \
			pythiaHBT::output_Bjorken_variables="${boolVal[$storeBjorkenCoordinates]}"

		# check and report whether run was successful
		runSuccess=`echo $?`
		check_success 'Pythia' $runSuccess '-' '| | - '

	fi


	# Get the filenames which need to be processed
	recordOfOutputFilenames_Sxp=$PYTHIA_RESULTS_DIRECTORY/`echo $collisionSystemStem`"_S_x_p_filenames.dat"
	recordOfOutputFilename_mult=$PYTHIA_RESULTS_DIRECTORY/`echo $collisionSystemStem`"_total_N_filename.dat"

	rm -f $HBT_EVENT_GEN_DIRECTORY/catalogue.dat
	rm -f $HBT_SV_DIRECTORY/catalogue.dat

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
	echo $projectile $target $beamEnergy $Nevents > $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
	readlink -f $PYTHIA_RESULTS_DIRECTORY/`cat $recordOfOutputFilename_mult` >> $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
	echo $projectile $target $beamEnergy $Nevents > $HBT_SV_DIRECTORY/ensemble_catalogue.dat
	readlink -f $PYTHIA_RESULTS_DIRECTORY/`cat $recordOfOutputFilename_mult` >> $HBT_SV_DIRECTORY/ensemble_catalogue.dat

	#cp $PYTHIA_DIRECTORY/main_BEeffects.cmnd $PYTHIA_RESULTS_DIRECTORY

)

echo '| - '`basename "$0"`': Finished everything!'
echo '| ------------------------------------------------------------'


# End of file
