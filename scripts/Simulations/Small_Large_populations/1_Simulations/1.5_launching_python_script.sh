#! /bin/bash

#this script runs on my computer

#List the trees from 1.4
trees=$(find ./TreeSequences/ -name "*.tree")

msprime_processing() {

        #First argument is the tree from trees above
        tree=$1
        #Create a simID variable
        simID=$(basename -s '.tree' $tree | awk '{split($0,a,"_");print a[3]}')

        #Create a CHR variable
        CHR=$(basename -s '.tree' $tree | awk '{split($0,a,"_");print a[4]}')

        #Lauch the actual script
        python3 MSPRIME_Recap_mut_subset.py ${tree} ${simID} ${CHR}

}

export -f msprime_processing
parallel -j 5 ::: msprime_processing ::: $(for tree in $trees; do echo $tree; done)
