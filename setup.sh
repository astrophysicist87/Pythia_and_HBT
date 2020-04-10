#! /usr/bin/env bash

cd src/pythia8243


setupMode="short"
while getopts ":f" opt; do
  case ${opt} in
    f ) setupMode="full"
      ;;
  esac
done
shift $((OPTIND -1))

if [ $setupMode = "full" ]
then
	./setup_pythia.sh
else
	gmake
fi

cd -

./compile_all.sh $@
