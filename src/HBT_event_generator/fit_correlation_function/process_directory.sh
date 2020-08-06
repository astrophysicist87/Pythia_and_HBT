#! /usr/bin/env bash
#-------------------

directoryToProcess=$1
directoryToStoreResults=$2

mkdir -p $directoryToStoreResults

parametersFilepath=$directoryToProcess/parameters.dat
particleCatalogueFilepath=$directoryToProcess/particle_catalogue.dat
catalogueFilepath=$directoryToProcess/catalogue.dat

./run.sh $directoryToStoreResults \
         $parametersFilepath \
         $particleCatalogueFilepath \
         $catalogueFilepath &

#echo 'This is a test' > $directoryToStoreResults/test.out
