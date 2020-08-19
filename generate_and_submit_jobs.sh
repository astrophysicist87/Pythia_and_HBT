#! /usr/bin/env bash

(
	# info for this run
	NTHREADS=32
	#DIRECTORY=results_pp_7TeV_Nev30000000_V0A
	DIRECTORY=results_pp7TeV_PbPb2760GeV
	NICENESS=19

	# Set job specifications here
	#declare -a class_ranges=("1-11" "12-16" "17-22" "23-28" "29-34" "35-41" "42-51" "52-151")
	#declare -a class_ranges=("0-20%" "20-40%" "40-60%" "60-100%")
	declare -a class_ranges=("0-100%")
	declare -a specs=(
		##############
		#blixen
		##############
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=1000 BEeffects="on" BEEnhancementMode="0" useDistribution="on"'
		##############
		#egan
		##############
		#'projectile="p" target="p" beamEnergy="7000" Nevents=30000000 SetPartonVertices="on"'
		#'projectile="p" target="p" beamEnergy="7000" Nevents=1000000'
		##############
		#lema
		##############
		'projectile="p" target="Pb" beamEnergy="7000" Nevents=1000000 BEeffects="on"'
		'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="1" compensationSet="1" compensationMode="1" compensationVersion="0"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="1" compensationSet="2" compensationMode="1" compensationVersion="0"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="0" compensationSet="1" compensationMode="1" compensationVersion="0"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="1" compensationSet="1" compensationMode="1" compensationVersion="3"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="1" compensationSet="2" compensationMode="1" compensationVersion="3"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="0" compensationSet="1" compensationMode="1" compensationVersion="3"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="1" compensationSet="1" compensationMode="2" compensationVersion="0"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="1" compensationSet="2" compensationMode="2" compensationVersion="0"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="0" compensationSet="1" compensationMode="2" compensationVersion="0"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="1" compensationSet="1" compensationMode="2" compensationVersion="3"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="1" compensationSet="2" compensationMode="2" compensationVersion="3"'
		#'projectile="Pb" target="Pb" beamEnergy="2760" Nevents=10000 BEeffects="on" shiftingSet="0" compensationSet="1" compensationMode="2" compensationVersion="3"'
	)

	#-----------------------------------------------------
	# Shouldn't have to change anything below this point
	#-----------------------------------------------------
	# Make sure all scripts are executable!
	find . -name "*.sh" | xargs chmod 755

	# Compile source code
	./compile_all.sh $NTHREADS &> compile_all.out

	# Export job specifications
	export niceness=$NICENESS
	source scripts/env.sh
	source scripts/export_specs.sh
	SPECS_FULL_PATH=`readlink -f scripts/specs.sh`
	export_specs_array "${specs[@]}" > $SPECS_FULL_PATH
	export_class_range_array "${class_ranges[@]}" >> $SPECS_FULL_PATH

	# Generate jobs for these specifications
	DIRECTORY_FULL_PATH=`readlink -f $DIRECTORY`
	./scripts/generate_jobs.sh $NTHREADS $DIRECTORY_FULL_PATH

	# Submit generated jobs
	./scripts/submit_jobs.sh $DIRECTORY_FULL_PATH

) &> /dev/null &

# End of file
