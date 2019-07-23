#! /usr/bin/env bash

cd src/pythia8235

gmake

cd -

./compile_all.sh $@
