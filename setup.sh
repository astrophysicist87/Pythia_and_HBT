#! /usr/bin/env bash

cd src/pythia8243

./setup_pythia.sh

cd -

./compile_all.sh $@
