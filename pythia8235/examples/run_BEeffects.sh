#! /usr/bin/env bash

# make sure results directory exists
DIRECTORY="$5"
if [ ! -d "$DIRECTORY" ]; then
	mkdir -p $DIRECTORY
fi

echo 'Running ./main_BEeffects' $@

# time and run
nohup time ./main_BEeffects $@\
			1> $DIRECTORY/main_BEeffects.out\
			2> $DIRECTORY/main_BEeffects.err
