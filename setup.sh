#! /usr/bin/env bash
#-------------------

#===================
# Load header info
source scripts/env.sh

# Compile Pythia source
cd $PYTHIA_DIRECTORY/..

# Process command-line options
setupMode="short"
while getopts ":fh" opt; do
  case ${opt} in
    f ) setupMode="full"
      ;;
	h ) echo "Usage: ./setup.sh [-fh] [number of OMP threads]"
      ;;
    \? ) echo "Usage: ./setup.sh [-fh] [number of OMP threads]"
      ;;
  esac
done
shift $((OPTIND -1))

# Decide whether to do short (default) or full compilation
if [ $setupMode = "full" ]
then
	( cd $HEPMC_DIRECTORY; ./setup_HepMC.sh; cd - )
	./setup_pythia.sh
else
	gmake
fi

# Go back
cd -

# Compile everything else
./compile_all.sh $@
