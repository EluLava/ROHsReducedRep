#! /bin/bash

#This script runs on my COMPUTER

simulations=`basename -s '.txt' $(find ./Simulations/mating*) | awk '{split($0,a,"_");print a[3]}'`

slim_applying_ped() {

        #Getwd
        cwd=$(pwd)

        slim -d "simID='$1'" -d chr=$2 -d replicate=$3 -d "cwd='${cwd}'" ./1.4_Applying_pedigree_nonWF.slim

}

export -f slim_applying_ped

#Creating Output directory
mkdir -p ./TreeSequences

parallel --dryrun --env ::: slim_applying_ped ::: $(for sim in $simulations; do echo $sim; done) ::: $(echo {1..30}) > command1
parallel -j 6 :::: command1 :::+ $(echo {1..300})

rm ./command1