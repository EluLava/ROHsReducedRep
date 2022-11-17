#! /bin/bash

trees=$(find ./Slim_Applying_inbreeding/*.trees)

mkdir -p ./MSprime_Recap_Mut_VCF

msprime_processing() {

        #First argument is the tree from trees above
        tree=$1
        #CHR variable
        CHR=$(basename -s '.trees' ${tree} | awk '{split($0,a,"_");print a[3]}')

        #Create a SimID variable
        simID='Sim17390'

        #Lauch the actual script
        python3 ./MSPRIME_mut_subset.py ${tree} ${simID} ${CHR}

}

export -f msprime_processing
parallel -j29 ::: msprime_processing ::: $(for tree in $trees; do echo $tree; done)
