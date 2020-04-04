#! /usr/bin/env bash

export chosen_Pythia_walltime='12:00:00'
export chosen_HBT_walltime_per_centrality='24:00:00'

generate_pbs () {

echo '#! /usr/bin/env bash

#PBS -l walltime='$1'
#PBS -l nodes=1:ppn='$2'
#PBS -A PAS0254
#PBS -j oe

cd "$PBS_O_WORKDIR" || exit $?

chosen_OMP_NUM_THREADS='$2
echo 'echo '\''export chosen_OMP_NUM_THREADS='\''$chosen_OMP_NUM_THREADS > omp_env.sh'

echo '
#./compile_all.sh $chosen_OMP_NUM_THREADS &> compile_all.out

if [ "$?" -eq "0" ]
then'
echo '	'${@:3}
echo 'fi
'
echo 'echo '\''Finished everything'\'

}

export -f generate_pbs
