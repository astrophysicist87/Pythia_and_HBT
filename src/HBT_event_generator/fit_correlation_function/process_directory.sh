#! /usr/bin/env bash
#-------------------

directoryToProcess=$1
directoryToStoreResults=$2

rm -rf $directoryToStoreResults
mkdir -p $directoryToStoreResults

parametersFilepath=$directoryToProcess/parameters.dat
particleCatalogueFilepath=$directoryToProcess/particle_catalogue.dat
catalogueFilepath=$directoryToProcess/catalogue.dat

./run.sh $directoryToStoreResults \
         $parametersFilepath \
         $particleCatalogueFilepath \
         $catalogueFilepath \
         format_with_pairs=1 &

#echo 'This is a test' > $directoryToStoreResults/test.out
