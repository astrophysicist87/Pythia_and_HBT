#! /bin/bash
#-----------

declare -a class_ranges=("1-11" "12-16" "17-22" "23-28" "29-34" "35-41" "42-51" "52-151")

study_S_x_p_results_directory=study_S_x_p_results
mkdir -p $study_S_x_p_results_directory

for eventClassCutString in "${class_ranges[@]}"
do
    qsub -F $eventClassCutString study_S_x_p_multClass.pbs
done
