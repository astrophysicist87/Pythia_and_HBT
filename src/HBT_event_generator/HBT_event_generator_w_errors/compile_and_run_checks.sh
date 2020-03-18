#! /usr/bin/env bash

(
	cd ../../../../test_shifter
	git pull && gmake distclean && gmake all

	cd -
	git pull && gmake distclean && gmake all

	./run_checks.sh 5 0.02

	./run_checks.sh 5 0.2

	./run_checks.sh 5 2.0

	./run_checks.sh 10 0.01

	./run_checks.sh 10 0.1

	./run_checks.sh 10 1.0


) &> compile_and_run_checks.out &
