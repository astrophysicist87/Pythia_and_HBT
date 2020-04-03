#! /usr/bin/env bash

variable=19870426

for var in "$@"
do
    export "$var"
done

./simple_test2.sh

echo $variable
