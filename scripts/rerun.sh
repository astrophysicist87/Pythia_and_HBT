#! /usr/bin/env bash

# Overwrite current defaults file with settings of previous run
mv settings.sh defaults.sh

# Turns Pythia off by default, but can be turned
# back on by including runPythia=true in "$@"
nohup ./driver.sh runPythia=false $@ &>> driver.out &
