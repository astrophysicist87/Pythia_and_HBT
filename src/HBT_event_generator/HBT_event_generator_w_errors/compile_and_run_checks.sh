#! /usr/bin/env bash

(
	# Command-line argument names and creates results directory
	RESULTSDIRECTORY=$1

	rm -rf $RESULTSDIRECTORY
	mkdir $RESULTSDIRECTORY
	echo $RESULTSDIRECTORY > ./resultsDirectory.dat
	cp ../parameters.dat $RESULTSDIRECTORY

	cd ../../../../test_shifter
	git pull && gmake distclean && gmake all

	cd -
	git pull && gmake distclean && gmake all && \
	cp ../parameters.dat . && time ./run_checks.e &> run_checks.out
) &> compile_and_run_checks.out &
