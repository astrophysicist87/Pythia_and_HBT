#! /usr/bin/env bash

cd src/pythia8243

gmake

cd -

./compile_all.sh $@
