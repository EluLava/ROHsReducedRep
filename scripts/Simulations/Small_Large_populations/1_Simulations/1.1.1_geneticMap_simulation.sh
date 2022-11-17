#! /bin/bash

#THIS SCRIPT RUNS ON MY COMPUTER

#source `which env_parallel.bash`

fregene_function() {

        replicate=$1

        #Create directory if needed
        mkdir -p ./Fregene/Input
        mkdir -p ./Fregene/Output

        #HERE Input files must be moved to Input director (from Fregene downloaded folder)

        fregene -i ./Fregene/Input/in_100Mb.xml -recomb ./Fregene/Input/recomb_100Mb.xml -p ./Fregene/Output/Fregeneoutput${replicate}
        rm ./Fregene/Output/Fregeneoutput${replicate}

}

export -f fregene_function
parallel --env  ::: fregene_function ::: $(echo {1..300})

#Create directory for slim (script 1.2)
mkdir -p ./Simulations