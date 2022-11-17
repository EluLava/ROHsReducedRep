#! /bin/bash

#Create output directory
mkdir -p ./Slim_creating_inbreeding

launching_slim() {

	seed=$1

	slim -d seed=${seed} ./Cattle_creating_ped.slim

}

export -f launching_slim

#CLUSTER
parallel -j1 ::: launching_slim ::: 17390