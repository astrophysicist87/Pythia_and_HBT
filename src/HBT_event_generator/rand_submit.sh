#! /usr/bin/env bash

(
	./clean_all.sh
	./compile_all.sh &> compile_all.out

	readlink -f ./HBT_particle.dat > ./HBT_event_generator_w_errors/particle_catalogue.dat
	cp parameters.dat ./HBT_event_generator_w_errors/

	readlink -f ./HBT_particle.dat > ./fit_correlation_function/particle_catalogue.dat
	cp parameters.dat ./fit_correlation_function/

	readlink -f ./HBT_particle.dat > ./source_variances/particle_catalogue.dat
	cp parameters.dat ./source_variances/
	
	cd HBT_event_generator_w_errors
	mkdir results
	export OMP_NUM_THREADS=12
	./run_HBT_event_generator.e \
		file_mode=0 \
		RNG_mult=1000 \
		RNG_Nev=100 \
		RNG_nLoops=1 \
		1> HBT_event_generator.out \
		2> HBT_event_generator.err

	cp -r results ../fit_correlation_function/
	cd ../fit_correlation_function/
	readlink -f results/HBT_pipiCF.dat > catalogue.dat
	./run_fit_correlation_function.e \
		1> fit_correlation_function.out \
		2> fit_correlation_function.err
	
	cd ../source_variances
	mkdir results
	./SV.e file_mode=0 \
		RNG_mult=1000 \
		RNG_Nev=1000 \
		RNG_nLoops=1\
		1> SV.out 2> SV.err

	# back to home directory
	cd ../

	resultsDirec=rand_submit_results_run3
	mkdir $resultsDirec
	cp parameters.dat $resultsDirec

	mv HBT_event_generator_w_errors/HBT_event_generator.* $resultsDirec/
	mv HBT_event_generator_w_errors/results $resultsDirec/HBTeg_results

	mv fit_correlation_function/fit_correlation_function.* $resultsDirec/
	mv fit_correlation_function/results $resultsDirec/fit_CF_results

	mv source_variances/SV.out $resultsDirec/
	mv source_variances/SV.err $resultsDirec/
	mv source_variances/results $resultsDirec/SV_results




) &

# End of file
