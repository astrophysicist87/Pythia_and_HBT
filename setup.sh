#! /usr/bin/env bash
#-------------------

# Compile Pythia source
cd src/pythia8243

# Process command-line options
setupMode="short"
while getopts ":f" opt; do
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
	./setup_pythia.sh
else
	gmake
fi

# Go back
cd -

# Compile everything else
./compile_all.sh $@
