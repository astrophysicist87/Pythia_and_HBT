#! /bin/bash
#-----------

declare -a class_ranges=("1-11" "12-16" "17-22" "23-28" "29-34" "35-41" "42-51" "52-151")

study_S_x_p_results_directory=study_S_x_p_results
mkdir -p $study_S_x_p_results_directory

PythiaHBTDirectory=/scratch/blixen/plumberg/Pythia_and_HBT
topResultsDirectory=results_pp_7TeV_Nev60000000/job-1/results
workingDirectory=$PythiaHBTDirectory/$topResultsDirectory/Pythia_results/dataset_0

for eventClassCutString in "${class_ranges[@]}"
do
    eventClassCut=(`echo $eventClassCutString | sed 's/-/ /g' | sed 's/%//g'`)
    classMin=${eventClassCut[0]}
    classMax=${eventClassCut[1]}
    for KTmin in $(seq 0 0.1 0.8)
    do
        KTmax=`echo "$KTmin+0.1" | bc`
        resultsFilename=S_x_p_N${classMin}_${classMax}_${KTmin}_${KTmax}.dat
        outputFilename=S_x_p_N${classMin}_${classMax}_${KTmin}_${KTmax}.out
        ./study_S_x_p $classMin $classMax $KTmin $KTmax \
                      $workingDirectory/pp_7000GeV_Nev60000000_mult.dat \
                      $workingDirectory/pp_7000GeV_Nev60000000_{0..599}.dat \
                      1> $study_S_x_p_results_directory/$resultsFilename \
                      2> $study_S_x_p_results_directory/$outputFilename &
    done
    wait
done