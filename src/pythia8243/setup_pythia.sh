#! /usr/bin/env bash

gmake distclean

./configure --with-hepmc2=../HepMC/build

gmake -j4
