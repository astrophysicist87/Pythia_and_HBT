#! /usr/bin/env bash

chosen_OMP_NUM_THREADS=2
echo 'chosen_OMP_NUM_THREADS='$chosen_OMP_NUM_THREADS > omp_env.sh

( ./compile_all.sh $chosen_OMP_NUM_THREADS &> compile_all.out && ./driver.sh useParallel=true projectile="Pb" target="Pb" beamEnergy="2760.0" Nevents=10 ThermalOnly='true' SetPartonVertices='off' &> driver.out ) &
