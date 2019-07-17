#! /usr/bin/env bash

cd pythia8235

gmake

cd ..

./compile_all.sh $@
