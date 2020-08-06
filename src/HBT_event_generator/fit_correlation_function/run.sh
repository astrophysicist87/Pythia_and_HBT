#! /usr/bin/env bash
#-------------------

# make sure results directory exists
#DIRECTORY=results
DIRECTORY=$1
if [ ! -d "$DIRECTORY" ]; then
	mkdir $DIRECTORY
fi

# time and run
nohup time ./run_fit_correlation_function.e "$@" \
				1> $DIRECTORY/fit_correlation_function.out \
				2> $DIRECTORY/fit_correlation_function.err
