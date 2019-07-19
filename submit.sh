#! /usr/bin/env bash

<<<<<<< HEAD
chosen_OMP_NUM_THREADS=2
echo 'chosen_OMP_NUM_THREADS='$chosen_OMP_NUM_THREADS > omp_env.sh

( ./compile_all.sh $chosen_OMP_NUM_THREADS &> compile_all.out && ./driver.sh useParallel=true projectile="Pb" target="Pb" beamEnergy="2760.0" Nevents=10 ThermalOnly='true' SetPartonVertices='off' &> driver.out ) &
=======
( ./compile_all.sh 12 &> compile_all.out && ./driver.sh &> driver.out ) &
>>>>>>> 591b1244348db8df19b43796c4cbabcd08133462
