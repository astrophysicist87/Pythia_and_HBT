#! /usr/bin/env bash

gmake distclean

./configure

gmake -j4
