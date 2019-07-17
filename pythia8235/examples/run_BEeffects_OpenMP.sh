#! /usr/bin/env bash

# make sure results directory exists
DIRECTORY="$5"
if [ ! -d "$DIRECTORY" ]; then
	mkdir -p $DIRECTORY
fi

echo 'Running ./main_BEeffects_OpenMP' $@

# time and run
nohup time ./main_BEeffects_OpenMP $@\
			1> $DIRECTORY/main_BEeffects_OpenMP.out\
			2> $DIRECTORY/main_BEeffects_OpenMP.err
