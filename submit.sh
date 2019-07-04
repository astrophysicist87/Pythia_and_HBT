#! /usr/bin/env bash

( ./compile_all.sh 12 &> compile_all.out && ./driver.sh &> driver.out ) &
