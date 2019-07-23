#! /usr/bin/env bash

for file in "$1"/job-*/submit.sh
do
	cd `dirname $file`
	pwd
	./submit.sh
	cd -
done

# End of file
