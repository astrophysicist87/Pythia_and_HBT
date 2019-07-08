#! /usr/bin/env bash

(
	rm -rf RESULTS
	rm compile_all.out driver.out submit.pbs.o*

	clean_directory () {
		rm $1/*.out $1/*.err $1/*.txt
		rm $1/*catalogue.dat
		rm $1/parameters.dat
		rm -rf $1/results
	}

	clean_directory pythia8235/examples

	clean_directory HBT_event_generator/HBT_event_generator_w_errors

	clean_directory HBT_event_generator/fit_correlation_function

	clean_directory HBT_event_generator/source_variances
) 2> /dev/null
