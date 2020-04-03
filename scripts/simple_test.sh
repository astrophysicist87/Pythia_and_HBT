#! /usr/bin/env bash

variable=19870426

for var in "$@"
do
    #eval "$var"
    "$var"
done

echo $variable
