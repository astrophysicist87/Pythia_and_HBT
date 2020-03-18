#! /usr/bin/env bash

(
	rm -rf results
	mkdir results

	cd ../../../../test_shifter
	git pull && gmake distclean && gmake all

	cd -
	git pull && gmake distclean && gmake all && \
	cp ../parameters.dat . && time ./run_checks.e &> run_checks.out
) &> compile_and_run_checks.out &
