#! /usr/bin/env bash

for file in "$1"/job-*/submit.pbs
do
	cd `dirname $file`
	pwd
	qsub submit.pbs
	cd -
done

# End of file