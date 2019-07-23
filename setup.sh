#! /usr/bin/env bash

cd src/pythia8235

./setup_pythia.sh

cd -

./compile_all.sh $@
