#!/bin/bash

generate_slurm () {

echo '#!/bin/bash

#SBATCH --time='$1'
#SBATCH --nodes=1
#SBATCH --ntasks-per-node='$2'
#SBATCH --job-name=slurm_run
#SBATCH --partition=qgp
#SBATCH --output=slurm_run.o%j
##SBATCH --error=slurm_run.e%j
##SBATCH --mail-user=plumberg@illinois.edu
##SBATCH --mail-type=BEGIN,END


#cd ${SLURM_SUBMIT_DIR} || exit $?
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

setup_env () {

case "$1" in
	qgp)
		export def_nthreads=12
		export def_pythia_walltime='12:00:00'
		export def_HBT_walltime='12:00:00'
		export def_queuename='qgp'
		;;

	*)
	echo $"Usage: $0 {qgp}"
	exit 1
 
esac

}

export -f generate_slurm
export -f generate_sh
export -f setup_env
