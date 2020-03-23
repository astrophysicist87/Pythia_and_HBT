#! /usr/bin/env bash

Rval=$1
aval=$2
stem=R`echo $Rval`_a`echo $aval`

RESULTSDIRECTORY=without_pair_density_Nev10000/results_${stem}
#RESULTSDIRECTORY=results_${stem}

rm -rf $RESULTSDIRECTORY
mkdir -p $RESULTSDIRECTORY
echo $RESULTSDIRECTORY > ./resultsDirectory.dat
cp ../parameters.dat .
cp ../parameters.dat $RESULTSDIRECTORY
time ./run_checks.e RNG_R=$Rval RNG_a=$aval &> $RESULTSDIRECTORY/run_checks_${stem}.out &

sleep 3
