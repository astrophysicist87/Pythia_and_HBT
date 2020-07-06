#! /usr/bin/env bash

tar xvzf HepMC-2.06.09.tar.gz

rm -rf build && mkdir build

cd build

cmake -DCMAKE_INSTALL_PREFIX=. \
      -Dmomentum:STRING=GEV \
      -Dlength:STRING=MM \
      ../HepMC-2.06.09

make
make test
make install
