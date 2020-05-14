#! /usr/bin/env bash

#for file in "$1"/job-*/submit*.pbs
for file in "$1"/job-*/submit.sh
do
	cd `dirname $file`
	pwd
	#qsub submit.pbs
	#qsub `basename $file`
	bash ./submit.sh
	cd -
done

# End of file
