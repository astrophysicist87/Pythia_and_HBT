#! /usr/bin/env bash

rm -rf HepMC-2.06.09 && tar xvzf HepMC-2.06.09.tar.gz

rm -rf build && mkdir build

cd build

cmake -DCMAKE_INSTALL_PREFIX=. \
      -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
      -Dmomentum:STRING=GEV \
      -Dlength:STRING=MM \
      ../HepMC-2.06.09

make
make test
make install
