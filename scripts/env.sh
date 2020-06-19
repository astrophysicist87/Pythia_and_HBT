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
export PYTHIA_DIRECTORY=$SOURCE_DIRECTORY/pythia${versionNumber}/pythiaHBT
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
echo 'export computeBEEnhancementExactly='$computeBEEnhancementExactly
echo
echo 'export shiftingSet='$shiftingSet
echo 'export compensationSet='$compensationSet
echo 'export compensationMode='$compensationMode
echo
echo 'export UseColorReconnection='$UseColorReconnection
echo 'export UseRopeHadronization='$UseRopeHadronization
echo 'export IncludeStringShoving='$IncludeStringShoving
echo 'export IncludeFlavourRopesMechanism='$IncludeFlavourRopesMechanism
echo 'export storeBjorkenCoordinates='$storeBjorkenCoordinates

}


output_this_job () {

declare -A boolVal=( ["true"]="1" ["false"]="0")

echo '--------------------------------------------------------------------------------'	# width is 80 spaces
echo '| Date              :' `date +%A" "%B" "%d" "%r" "%Z" "%Y`
echo '| Overview          :' $Nevents $projectile'+'$target 'collisions at \sqrt{s_{NN}} =' $beamEnergy 'GeV'
echo '--------------------------------------------------------------------------------'	# width is 80 spaces
echo '| -------------------------------------------'
echo '| Collision system  :' $projectile'+'$target
echo '| Beam energy (GeV) :' $beamEnergy
echo '| Number of events  :' $Nevents
printf "| ${eventClassSelectionMode^} classes:"
for eventClassCutString in "${class_ranges[@]}"
do
	if [ "$eventClassSelectionMode" == "centrality" ]
	then
	    printf " $eventClassCutString""%"
		#printf '%%'
	else
	    printf " $eventClassCutString"
	fi
done
printf "\n"
echo '| Particle IDs used in HBT:' $chosenHBTparticle
echo '|'
echo '| -------------------------------------------'
echo '| Chosen analyses :'
echo '|     * Pythia/Angantyr event generation    :' $runPythia
echo '|     * Compute HBT correlation function    :' $runHBTEG
echo '|     * Fit HBT correlation function        :' $runFitCF
echo '|     * Compute source variances            :' $runSV
echo '|'
echo '| -------------------------------------------'
echo '| Chosen Pythia/Angantyr settings           :'
echo '|     * Impact parameter range (fm)         : ('$bMin','$bMax')'
echo '|     * Set fragmentation vertices          :' $SetFragmentationVertices
echo '|     * Set parton vertices                 :' $SetPartonVertices
echo '|     * Print thermal particles only        :' $ThermalOnly
echo '|     * Include color reconnection          :' $UseColorReconnection
echo '|     * Include rope hadronization          :' $UseRopeHadronization
echo '|     * Include string shoving              :' $IncludeStringShoving
echo '|     * Include flavor ropes mechanism      :' $IncludeFlavourRopesMechanism
echo '|     * Store Bjorken coordinates           :' $storeBjorkenCoordinates
echo '|'
echo '|     * Bose-Einstein effects included      :' $BEeffects
echo '|     * Bose-Einstein mode                  :' $BEEnhancementMode
echo '|       if Bose-Einstein mode == 0          :'
echo '|          - QRef default value             :' $QRefValue
echo '|          - estimate source size from'
echo '|            space-time distribution        :' $useDistribution
echo '|            -- use space-time instead'
echo '|               of spatial separation       :' $useInvariantSourceSize
echo '|            -- use relative separations'
echo '|               instead of coordinates      :' $useRelativeDistance
echo '|            -- evaluate distances in PRF'
echo '|               instead of lab frame        :' $useRestFrame
echo '|       if Bose-Einstein mode == 1          :'
echo '|          - include phase space            :' $includePhaseSpace
echo '|          - include estimated pair density :' $linearInterpolateCDF
echo '|          - compute shift integral exactly :' $computeBEEnhancementExactly
echo '|          - shifting set definition        :' $shiftingSet
echo '|          - compensation set definition    :' $compensationSet
echo '|          - compensation mode              :' $compensationMode
echo '|          - compensation version           :' $compensationVersion
echo '|'
echo '| -------------------------------------------'
echo '| Chosen HBT analysis settings              :'
echo '|     * Number of KT points                 :' `cat parameters.dat | grep n_KT_pts | awk -F= '{print $2}'`
echo '|          - KT minimum                     :' `cat parameters.dat | grep KTmin | awk -F= '{print $2}'`
echo '|          - KT maximum                     :' `cat parameters.dat | grep KTmax | awk -F= '{print $2}'`
echo '|     * Number of KPhi points               :' `cat parameters.dat | grep n_Kphi_pts | awk -F= '{print $2}'`
echo '|     * Number of KL points                 :' `cat parameters.dat | grep n_KL_pts | awk -F= '{print $2}'`
echo '|          - KL minimum                     :' `cat parameters.dat | grep KLmin | awk -F= '{print $2}'`
echo '|          - KL maximum                     :' `cat parameters.dat | grep KLmax | awk -F= '{print $2}'`
echo '|     * Relative momentum q mode            :' `cat parameters.dat | grep q_mode | awk -F= '{print $2}'`
echo '|       if q mode == 0                      :'
echo '|          - Number of q points (out)       :' `cat parameters.dat | grep n_qo_pts | awk -F= '{print $2}'`
echo '|          - Number of q points (side)      :' `cat parameters.dat | grep n_qs_pts | awk -F= '{print $2}'`
echo '|          - Number of q points (long)      :' `cat parameters.dat | grep n_ql_pts | awk -F= '{print $2}'`
echo '|          - Bin width (out)                :' `cat parameters.dat | grep delta_qo | awk -F= '{print $2}'`
echo '|          - Bin width (side)               :' `cat parameters.dat | grep delta_qs | awk -F= '{print $2}'`
echo '|          - Bin width (long)               :' `cat parameters.dat | grep delta_ql | awk -F= '{print $2}'`
echo '|       if q mode == 1                      :'
echo '|          - Scalar type                    :' `cat parameters.dat | grep scalar_mode | awk -F= '{print $2}'`
echo '|          - Number of Q points             :' `cat parameters.dat | grep n_Q_pts | awk -F= '{print $2}'`
echo '|          - Bin width (scalar)             :' `cat parameters.dat | grep delta_Q | awk -F= '{print $2}'`
echo '|          - Points in qRP radial integral  :' `cat parameters.dat | grep n_qRP_pts | awk -F= '{print $2}'`
echo '|          - Points in qRP angular integral :' `cat parameters.dat | grep n_thq_pts | awk -F= '{print $2}'`
echo '|     * Method of computing correlator      :' `cat parameters.dat | grep method_mode | awk -F= '{print $2}'`
echo '|'
echo '| -------------------------------------------'
echo '| Chosen settings to fit correlator         :'
echo '|     * Fitting mode                        :' `cat parameters.dat | grep fit_mode | awk -F= '{print $2}'`
echo '|     * Include cross terms                 :' `cat parameters.dat | grep include_cross_terms | awk -F= '{print $2}'`
echo '|     * Use q-axis slices only              :' `cat parameters.dat | grep use_slices_only | awk -F= '{print $2}'`
echo '|'
echo '| -------------------------------------------'
}

export -f check_success
export -f clean_directory
export -f output_settings
export -f output_this_job

# End of file
