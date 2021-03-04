#! /usr/bin/env bash

# Overwrite current defaults file with settings of previous run
mv settings.sh scripts/defaults.sh

# Width is 80 spaces
echo '<<<==========================================================================>>>' >> driver_pbs.out
echo '<<<======================== RESUBMISSION STARTS HERE ========================>>>' >> driver_pbs.out
echo '<<<==========================================================================>>>' >> driver_pbs.out

# Turns Pythia off by default, but can be turned
# back on by including runPythia=true in "$@"
nohup ./driver_pbs.sh runPythia=false $@ &>> driver_pbs.out &
