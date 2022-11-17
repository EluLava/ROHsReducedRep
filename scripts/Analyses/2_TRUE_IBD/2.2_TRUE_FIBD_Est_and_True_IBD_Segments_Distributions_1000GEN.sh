#! /bin/bash

#Define the function
TRUE_IBD_F_ROHs_Dist() {

    #First argument passed to the function is the IBD segments file (per Sim)
    FILE=$1
    #Extract SimID from VCF filename
    Sim_ID=$(basename -s ".txt" ${FILE} | cut -d'_' -f4)
    
    #Open R for F and ROHs dist estimation
    R --vanilla <<EOF

        library(doParallel)
        library(foreach)
        cl <- 10
        registerDoParallel(cl)

        # read the file
        dta = read.table("${FILE}", header = T)

        ## F TRUE IBD

        #Get sum of genome within IBD seg per INDV
        FtrueIBD = aggregate((POS2 - POS1) ~ ID, data = dta, sum)
        #Divide SUM within IBD seg by genome length
        FtrueIBD[,2] = FtrueIBD[,2]/3000000000
        #Add SimID column
        FtrueIBD.2 = cbind(SimID=rep("${Sim_ID}",nrow(FtrueIBD)), FtrueIBD)
        #change columns names
        colnames(FtrueIBD.2) = c("SimID", "ID", "F_true_IBD_GEN1000")

        #save dataframe
        write.table(FtrueIBD.2, "./TRUE_IBD/1000GEN/stats/FTRUEIBD_${Sim_ID}_GEN1000.txt", quote = F, col.names = F, row.names = F)

        ## ROHs DIST

        #add one column for "length" of the IBD segment
        dta_dist = cbind(dta, LENGTH=vector(mode = "numeric", length = nrow(dta)))
        #Add length
        dta_dist[,6] = dta_dist[,5] - dta_dist[,4]

        #Create the output df
        dta_roh = as.data.frame(matrix(nrow=6,ncol=6))
        colnames(dta_roh) = c("SimID","WINDOWS.NB","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")

        #Fill SimuID, CALL_METH, REP, PERC
        dta_roh[,1] = rep("${Sim_ID}",nrow(dta_roh))
        dta_roh[,2] = rep("TRUE_IBD_1000",nrow(dta_roh))
        dta_roh[,3] = rep("1",nrow(dta_roh))
        dta_roh[,4] = c(1:6)

        #Loop through ROHs and Add the roh class
        foreachoutput = foreach(i=1:nrow(dta_dist), .combine=rbind) %dopar% {

            #Create output
            output = as.data.frame(matrix(nrow=1,ncol=2))
            colnames(output) = c("i","ROH_CLASS")
            output[1,1] = i

            #Define ROH class
            if(dta_dist[i,6] < 2000000){
                output[1,2] = 1

            } else if(dta_dist[i,6] >= 2000000 & dta_dist[i,6] < 4000000) {
                output[1,2] = 2

            } else if(dta_dist[i,6] >= 4000000 & dta_dist[i,6] < 6000000) {
                output[1,2] = 3

            } else if(dta_dist[i,6] >= 6000000 & dta_dist[i,6] < 10000000) {
                output[1,2] = 4

            } else if(dta_dist[i,6] >= 10000000 & dta_dist[i,6] < 16000000) {
                output[1,2] = 5

            } else if(dta_dist[i,6] >= 16000000) {
                output[1,2] = 6
            }

            #return the rownb and class
            return(output)
        }

        #Sort output of foreach per i
        foreachoutput.2 = foreachoutput[order(foreachoutput[,1]),]
        #Merge with dta_dist
        dta_dist.2 = cbind(dta_dist, ROH_CLASS=foreachoutput.2[,2])

        #For each ROH category estimate sum of lengths and mean sum of length for each individual
        x = aggregate(dta_dist.2[,6], by = list(INDIV=dta_dist.2[,2], CLASS=dta_dist.2[,7]), FUN = sum)
        y = aggregate(dta_dist.2[,6], by = list(INDIV=dta_dist.2[,2], CLASS=dta_dist.2[,7]), FUN = length)

        #Count the number of ROHs and and sum Length for each catgory
        for(roh_class in 1:6) {
            dta_roh[roh_class,6] = sum(x[3][x[2] == roh_class])/(nrow(FtrueIBD.2))
            dta_roh[roh_class,5] = sum(y[3][y[2] == roh_class])/(nrow(FtrueIBD.2))
        }

        #write output
        write.table(dta_roh, file = "./TRUE_IBD/1000GEN/stats/ROHs_Dist_${Sim_ID}_GEN1000.txt", quote = F, row.names = F, col.names = F)

EOF

}

export -f TRUE_IBD_F_ROHs_Dist

#Loop through simualations replicates
for f in $(find ./TRUE_IBD/1000GEN/stats/ -name "TRUE_IBD_SEGMENTS_*_GEN1000.txt"); do TRUE_IBD_F_ROHs_Dist ${f}; done

# Then just merge the SIM replicates into one dataframe

## F TRUE IBD
#Create header
echo -e "SimID\tindividuals\tFroh_TRUE_IBD_1000GEN" > ./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt
#Add all files
for file in $(find ./TRUE_IBD/1000GEN/stats/ -name "FTRUEIBD_*_GEN1000.txt"); do cat ${file} >> ./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt; done
#Rm per SimID file
rm ./TRUE_IBD/1000GEN/stats/FTRUEIBD_*_GEN1000.txt

## ROHs Dist
#Create header
echo -e "SimID\tWINDOWS.NB\tREP\tROHs_CLASS\tMean_NROH_per_ind\tMean_Tot_Kb_per_ind" > ./TRUE_IBD/1000GEN/stats/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt
#Add all files
for file in $(find ./TRUE_IBD/1000GEN/stats/ -name "ROHs_Dist_*_GEN1000.txt"); do cat ${file} >> ./TRUE_IBD/1000GEN/stats/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt; done
#Rm per SimID file
rm ./TRUE_IBD/1000GEN/stats/ROHs_Dist_*_GEN1000.txt

