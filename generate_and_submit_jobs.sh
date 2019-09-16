#! /usr/bin/env bash

(
	./generate_jobs.sh 1 this_is_a_test

	./submit_jobs.sh this_is_a_test
) &> /dev/null &

# End of file
