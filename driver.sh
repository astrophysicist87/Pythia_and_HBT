#! /usr/bin/env bash

# Load header info
source env.sh

# Load default variable values
source defaults.sh

# Update any variables set from the command line
for var in "$@"
do
    eval "$var"
done

# make sure main results directory exists
if [ ! -d "$MAIN_RESULTS_DIRECTORY" ]
then
	mkdir $MAIN_RESULTS_DIRECTORY
	echo 'Created' $MAIN_RESULTS_DIRECTORY
fi

#===================
# Main calculation
#===================

echo 'Processing Nevents =' $Nevents $projectile'+'$target 'collisions at' $beamEnergy 'GeV'

nCC=0
for centralityCutString in "0-100%"
do

	# make sure current results directory exists
	nextCurrentResultsDirectoryName=`get_dirname $CURRENT_RESULTS_DIRECTORY '-' "true"`
	CURRENT_RESULTS_DIRECTORY=$nextCurrentResultsDirectoryName
	mkdir $CURRENT_RESULTS_DIRECTORY
	echo 'Created' $CURRENT_RESULTS_DIRECTORY

	# set names of sub-directories
	PYTHIA_RESULTS_DIRECTORY=$CURRENT_RESULTS_DIRECTORY/Pythia_results
	HBT_RESULTS_DIRECTORY=$CURRENT_RESULTS_DIRECTORY/HBT_results

	# make sure HBT results directory exists
	if [ ! -d "$HBT_RESULTS_DIRECTORY" ]
	then
		mkdir $HBT_RESULTS_DIRECTORY
		echo 'Created' $HBT_RESULTS_DIRECTORY
	fi

	success=0
	echo '  -- analyzing centrality class' $centralityCutString

	centralityCut=(`echo $centralityCutString | sed 's/-/ /g' | sed 's/%//g'`)
	lowerLimit=${centralityCut[0]}
	upperLimit=${centralityCut[1]}

	collisionSystemStem=$projectile$target"_"`echo $beamEnergy`"GeV_Nev"$Nevents
	collisionSystemCentralityStem=$projectile$target"_"`echo $beamEnergy`"GeV_C"$lowerLimit"_"$upperLimit"_Nev"$Nevents

	#=====================================
	# Run Pythia (if desired)
	#PythiaSuccess=$(
	(

		cd $PYTHIA_DIRECTORY
		echo '     Now in '`pwd`

		if $runPythia && [ "$nCC" -eq "0" ]
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


			# time and run
			if $useParallel
			then
				./run_BEeffects_OpenMP.sh $projectile $target $beamEnergy \
									$Nevents $PYTHIA_RESULTS_DIRECTORY \
									$lowerLimit $upperLimit $ThermalOnly
			else
				./run_BEeffects.sh $projectile $target $beamEnergy \
									$Nevents $PYTHIA_RESULTS_DIRECTORY \
									$lowerLimit $upperLimit $ThermalOnly
			fi

			# check and report whether run was successful
			runSuccess=`echo $?`
			check_success 'Pythia' $runSuccess

			# copy results
			#cp $PYTHIA_RESULTS_DIRECTORY/* $RESULTS_DIRECTORY/
		fi

		# Get the filenames which need to be processed
		recordOfOutputFilenames_Sxp=$PYTHIA_RESULTS_DIRECTORY/`echo $collisionSystemStem`"_S_x_p_filenames.dat"
		recordOfOutputFilename_mult=$PYTHIA_RESULTS_DIRECTORY/`echo $collisionSystemStem`"_total_N_filename.dat"
		rm $HBT_EVENT_GEN_DIRECTORY/catalogue.dat
		rm $HBT_SV_DIRECTORY/catalogue.dat
		for line in `cat $recordOfOutputFilenames_Sxp`
		do
			readlink -f $PYTHIA_RESULTS_DIRECTORY/$line >> $HBT_EVENT_GEN_DIRECTORY/catalogue.dat
			readlink -f $PYTHIA_RESULTS_DIRECTORY/$line >> $HBT_SV_DIRECTORY/catalogue.dat
		done
		# Set particle catalogue
		readlink -f $PYTHIA_RESULTS_DIRECTORY/HBT_particle.dat > $HBT_EVENT_GEN_DIRECTORY/particle_catalogue.dat
		readlink -f $PYTHIA_RESULTS_DIRECTORY/HBT_particle.dat > $HBT_FITCF_DIRECTORY/particle_catalogue.dat
		readlink -f $PYTHIA_RESULTS_DIRECTORY/HBT_particle.dat > $HBT_SV_DIRECTORY/particle_catalogue.dat
		# Set ensemble catalogue
		#echo $projectile $target $beamEnergy 0 100 $Nevents > $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
		echo $projectile $target $beamEnergy $lowerLimit $upperLimit $Nevents > $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
		readlink -f $PYTHIA_RESULTS_DIRECTORY/`cat $recordOfOutputFilename_mult` >> $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
		echo $projectile $target $beamEnergy $lowerLimit $upperLimit $Nevents > $HBT_SV_DIRECTORY/ensemble_catalogue.dat
		readlink -f $PYTHIA_RESULTS_DIRECTORY/`cat $recordOfOutputFilename_mult` >> $HBT_SV_DIRECTORY/ensemble_catalogue.dat

		#exit $runSuccess
	)

	#=====================================
	# Run HBT_event_generator
	#HBTegSuccess=$(
	(

		cd $HBT_EVENT_GEN_DIRECTORY
		echo '     Now in '`pwd`

		# using OpenMP (leave a couple cores free)
		export OMP_NUM_THREADS=40

		if [ ! -d "./results" ]; then
			    mkdir results
		fi

		cp ../parameters.dat .

		#default: assume BE effects are turned off in Pythia
		chosen_BE_mode=0
		if [ "$BEeffects" == "on" ]; then
				# BE_mode value for enhancement mode = 0 or 1
				chosen_BE_mode=1
				# handle the special case
			    if [ "$BEEnhancementMode" == "2" ]; then
					chosen_BE_mode=2
				fi
		fi

		# time and run
		nohup time ./run_HBT_event_generator.e \
				centrality_minimum=$lowerLimit \
				centrality_maximum=$upperLimit \
				BE_mode=$chosen_BE_mode \
				1> HBT_event_generator.out \
				2> HBT_event_generator.err
		# N.B. - centralities now determined in Pythia

		# check and report whether run was successful
		runSuccess=`echo $?`
		check_success 'HBT_event_generator' $runSuccess

		# copy results
		cp HBT_event_generator.[oe]* ./results
		mkdir $HBT_RESULTS_DIRECTORY/CF_results
		cp -r ./results/* $HBT_RESULTS_DIRECTORY/CF_results

		readlink -f ./results/HBT_pipiCF.dat > $HBT_FITCF_DIRECTORY/catalogue.dat

		#exit $runSuccess
	)

	#=====================================
	# Run fit_correlation_function
	#fitCFSuccess=$(
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
		mkdir $HBT_RESULTS_DIRECTORY/fit_results
		cp -r ./results/* $HBT_RESULTS_DIRECTORY/fit_results

		#exit $runSuccess
	)

	#=====================================
	# Run SV.e
	#SVSuccess=$(
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
				1> SV_record.out \
				2> SV_record.err

		# check and report whether run was successful
		runSuccess=`echo $?`
		check_success 'source_variances' $runSuccess

		# copy results
		cp SV_record.[oe]* ./results
		mkdir $HBT_RESULTS_DIRECTORY/SV_results
		cp -r ./results/* $HBT_RESULTS_DIRECTORY/SV_results

		#exit $runSuccess
	)

	# only store results and clean up if run was successful
	#success=$[PythiaSuccess + HBTegSuccess + fitCFSuccess]

	# check if all codes were successful
	#if [ "$success" -eq "0" ]
	#then
		#add a few more files
		cp $PYTHIA_DIRECTORY/main_BEeffects.cmnd $PYTHIA_RESULTS_DIRECTORY
		cp $HBT_DIRECTORY/parameters.dat $HBT_RESULTS_DIRECTORY

		#typeStem=""
		#if [ "$ThermalOnly" == 'true' ]
		#then
		#	typeStem="_THERMAL"
		#fi

		zipFilename=$CURRENT_RESULTS_DIRECTORY".zip"

		zip -r $zipFilename $CURRENT_RESULTS_DIRECTORY

		# Clean-up HBT directories (but not Pythia results directory!!!)
		#rm -rf $HBT_EVENT_GEN_DIRECTORY/*HBT_event_generator.[oe]* $HBT_EVENT_GEN_DIRECTORY/results\
		#	   $HBT_FITCF_DIRECTORY/*fit_correlation_function.[oe]* $HBT_FITCF_DIRECTORY/results
	#fi

	# move on to next centrality class
	nCC=$[nCC+1]

done	# all centralities finished

echo 'Finished everything!'


# End of file
