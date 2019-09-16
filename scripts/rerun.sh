#! /usr/bin/env bash

# Overwrite current defaults file with settings of previous run
mv settings.sh defaults.sh

nohup ./driver.sh $@ &>> driver.out &
