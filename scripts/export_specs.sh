#! /usr/bin/env bash

function export_specs_array () {

loc_specs=("$@")

echo 'declare -a specs=('
for spec in "${loc_specs[@]}"
do
    echo \'$spec\'
done
echo ')'

}

export -f export_specs_array
