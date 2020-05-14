#! /usr/bin/env bash

generate_pbs () {

echo '#! /usr/bin/env bash

#PBS -l walltime='$1'
#PBS -l nodes=1:ppn='$2'
#PBS -j oe

cd "$PBS_O_WORKDIR" || exit $?
'
#chosen_OMP_NUM_THREADS='$2
#echo 'echo '\''export chosen_OMP_NUM_THREADS='\''$chosen_OMP_NUM_THREADS > omp_env.sh'

echo '
if [ "$?" -eq "0" ]
then'
echo '	'${@:3}
echo 'fi
'
echo 'echo '\''Finished everything'\'

}

generate_sh () {

echo '#! /usr/bin/env bash
'
echo "$@"
echo 'echo '\''Finished everything'\'

}

export -f generate_pbs
export -f generate_sh
