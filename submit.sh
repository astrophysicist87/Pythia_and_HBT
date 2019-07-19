#! /usr/bin/env bash

chosen_OMP_NUM_THREADS=12
echo 'chosen_OMP_NUM_THREADS='$chosen_OMP_NUM_THREADS > omp_env.sh

( ./compile_all.sh		\
	$chosen_OMP_NUM_THREADS	\
	&> compile_all.out	\
  && ./driver.sh		\
	useParallel=false	\
	projectile="p"		\
	target="p"		\
	beamEnergy="13000.0"	\
	Nevents=1000		\
	ThermalOnly='true'	\
	SetPartonVertices='off'	\
	&> driver.out ) &
