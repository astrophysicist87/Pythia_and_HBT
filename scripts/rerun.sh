#! /usr/bin/env bash

# Overwrite current defaults file with settings of previous run
mv settings.sh scripts/defaults.sh

# Width is 80 spaces
echo '<<<==========================================================================>>>' >> driver.out
echo '<<<======================== RESUBMISSION STARTS HERE ========================>>>' >> driver.out
echo '<<<==========================================================================>>>' >> driver.out

# Turns Pythia off by default, but can be turned
# back on by including runPythia=true in "$@"
nohup ./driver.sh runPythia=false $@ &>> driver.out &
