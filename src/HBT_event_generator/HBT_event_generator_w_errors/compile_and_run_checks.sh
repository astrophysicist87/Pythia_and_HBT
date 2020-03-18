#! /usr/bin/env bash

(
	cd ../../../../test_shifter
	git pull && gmake distclean && gmake all

	cd -
	git pull && gmake distclean && gmake all

	./run_checks.sh 5 0.2

	./run_checks.sh 5 0.5

	./run_checks.sh 5 0.75

	./run_checks.sh 5 1.0

	./run_checks.sh 5 1.5

	./run_checks.sh 5 2.0

	./run_checks.sh 5 5.0

	./run_checks.sh 5 7.5

	./run_checks.sh 5 10.0

	./run_checks.sh 5 15.0

	./run_checks.sh 5 20.0

	wait

	./run_checks.sh 10 0.2

	./run_checks.sh 10 0.5

	./run_checks.sh 10 0.75

	./run_checks.sh 10 1.0

	./run_checks.sh 10 1.5

	./run_checks.sh 10 2.0

	./run_checks.sh 10 5.0

	./run_checks.sh 10 7.5

	./run_checks.sh 10 10.0

	./run_checks.sh 10 15.0

	./run_checks.sh 10 20.0


) &> compile_and_run_checks.out &
