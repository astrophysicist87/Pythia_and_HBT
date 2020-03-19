#! /usr/bin/env bash

Rval=$1
aval=$2
stem=R`echo $Rval`_a`echo $aval`

RESULTSDIRECTORY=results_${stem}

rm -rf $RESULTSDIRECTORY
mkdir $RESULTSDIRECTORY
echo $RESULTSDIRECTORY > ./resultsDirectory.dat
cp ../parameters.dat
cp ../parameters.dat $RESULTSDIRECTORY
time ./run_checks.e RNG_R=$Rval RNG_a=$aval &> run_checks_${stem}.out &

sleep 3
