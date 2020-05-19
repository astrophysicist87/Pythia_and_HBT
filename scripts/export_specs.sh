#! /usr/bin/env bash

function export_class_range_array () {

loc_class_ranges=("$@")

echo 'declare -a class_ranges=('
for CR in "${loc_class_ranges[@]}"
do
    echo \'$CR\'
done
echo ')'

}

function export_specs_array () {

loc_specs=("$@")

echo 'declare -a specs=('
for spec in "${loc_specs[@]}"
do
    echo \'$spec\'
done
echo ')'

}

export -f export_class_range_array
export -f export_specs_array
