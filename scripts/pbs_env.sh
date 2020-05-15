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

setup_env () {

case "$1" in
	mangi)
		export def_nthreads=128
		export def_pythia_walltime='12:00:00'
		export def_HBT_walltime='12:00:00'
		export def_queuename='amdsmall'
		;;

	mesabi)
		export def_nthreads=24
		export def_pythia_walltime='12:00:00'
		export def_HBT_walltime='12:00:00'
		export def_queuename='small'
		;;

	pitzer)
		export def_nthreads=40
		export def_pythia_walltime='12:00:00'
		export def_HBT_walltime='12:00:00'
		export def_account_string='PAS0254'
		;;

	vishnu)
		export def_nthreads=16
		export def_pythia_walltime='12:00:00'
		export def_HBT_walltime='12:00:00'
		;;

	*)
	echo $"Usage: $0 {mangi|mesabi|pitzer|vishnu}"
	exit 1
 
esac

}

export -f generate_pbs
export -f generate_sh
export -f setup_env
