#!/bin/bash

#Define a function for ROHs calling PLINK for WGS (and FHOM for SM)
WGS_ROHsCal_PLINK(){
        
        #First argument passed to the function is VCF file
	VCF=$1
        #Extract SimuID from VCFfile name
        SimID=$(basename -s ".vcf.gz" $VCF | cut -d '_' -f2)

        #Create output directory for ROHs called with PLINK if does not exists
        mkdir -p ./WGS_ROHs_PLINK
        #Same for FHOM as well
        mkdir -p ./WGS_FHOM

        #Change format from VCF to bed & Call ROHs      
        plink --vcf ${VCF} --chr-set 30 --make-bed --out ./WGS_ROHs_PLINK/BTMP_${SimID}

        #Call ROHs
        plink --bfile ./WGS_ROHs_PLINK/BTMP_${SimID} --homozyg --homozyg-window-snp 100 --homozyg-density 100 \
        --homozyg-snp 100 --homozyg-gap 1000 --homozyg-window-het 2 --homozyg-het 8 --homozyg-kb 100 --chr-set 30 --out ./WGS_ROHs_PLINK/WGS_${SimID}_ROH
        
        #plink FHOM
        plink --bfile ./WGS_ROHs_PLINK/BTMP_${SimID} --het --chr-set 30 --out ./WGS_FHOM/HET_WGS_${SimID}

        #rm BEDfile
        rm ./WGS_ROHs_PLINK/BTMP_${SimID}.*

        #Open R for F and ROHs distrbutions
        R --vanilla <<EOF

                #Read the roh indv plink file
                dtaroh_tmp = read.table("./WGS_ROHs_PLINK/WGS_${SimID}_ROH.hom.indiv", header = T)
                #Read the het plink file
                dtahet_tmp = read.table("./WGS_FHOM/HET_WGS_${SimID}.het", header = T)
                #Read the all roh plink file
                dtaroh_dist = read.table("./WGS_ROHs_PLINK/WGS_${SimID}_ROH.hom", header = T)[,-c(3,11:13)]

                ## FROH

                #Set the genome length IN KB
                genome_length = 3000000
                #Set individuals nb (use dtahet here in case some INDVs have zero ROHs and are not in PLINK .hom file, often the case with RADseq)
                nb_indiv=nrow(dtahet_tmp)

                #Set the individuals names vector
                individuals = c(1000, seq(1,nb_indiv-1))

                ###Calculate FROH
                Froh_tmp = apply(dtaroh_tmp, 1, function(x) as.numeric(as.character(x[5]))/genome_length)

                #Cbind INDV name, FROH, PLINK info such as total # of ROHs together in one dataframe
                F_sub = cbind(individuals, Froh_tmp, dtaroh_tmp[,4], dtaroh_tmp[,5])
                colnames(F_sub) = c("individuals", "Froh_WGS", "NROH_WGS", "SROH_WGS")

                #Adding SimID info
                WGS_dta = as.data.frame(cbind(rep("${SimID}",nrow(F_sub)), F_sub))
                colnames(WGS_dta) = c("SimuID", "individuals", "Froh_WGS", "NROH_WGS", "SROH_WGS")

                #Write the dataframe FROH
                write.table(WGS_dta, "./WGS_ROHs_PLINK/WGS_PLINK_${SimID}_FROH.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

                ## ROHs Dist

                #Create the output dataframe which will contain lengths per ROH category
                dta_roh = as.data.frame(matrix(nrow = 6, ncol = 6))
                colnames(dta_roh) = c("SimuID","WINDOWS.NB","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")

                #Fill SimuID, REP, and WINDOWS.NB
                dta_roh[,1] = rep("${SimID}",nrow(dta_roh))
                dta_roh[,2] = rep("WGS",nrow(dta_roh))
                dta_roh[,3] = rep("1",nrow(dta_roh)) #since WGS only 1
                dta_roh[,4] = c(1:6)

                #Add one column to plink df for ROH category
                dtaroh_dist = cbind(dtaroh_dist, rep(NA,nrow(dtaroh_dist)))
                colnames(dtaroh_dist)[10] = "ROH_CLASS"

                #Loop through ROHs in plink file and Add the roh class
                for(row in 1:nrow(dtaroh_dist)){

                        if(dtaroh_dist[row,8] < 2000){
                                dtaroh_dist[row,10] = 1

                        } else if(dtaroh_dist[row,8] >= 2000 & dtaroh_dist[row,8] < 4000) {
                                dtaroh_dist[row,10] = 2

                        } else if(dtaroh_dist[row,8] >= 4000 & dtaroh_dist[row,8] < 6000) {
                                dtaroh_dist[row,10] = 3

                        } else if(dtaroh_dist[row,8] >= 6000 & dtaroh_dist[row,8] < 10000) {
                                dtaroh_dist[row,10] = 4

                        } else if(dtaroh_dist[row,8] >= 10000 & dtaroh_dist[row,8] < 16000) {
                                dtaroh_dist[row,10] = 5

                        } else if(dtaroh_dist[row,8] >= 16000) {
                                dtaroh_dist[row,10] = 6
                        }
                }

                #For each ROH category estimate sum of lengths and mean sum of length for each individual
                x = aggregate(dtaroh_dist[,8], by = list(INDIV=dtaroh_dist[,2], CLASS=dtaroh_dist[,10]), FUN = sum)
                y = aggregate(dtaroh_dist[,8], by = list(INDIV=dtaroh_dist[,2], CLASS=dtaroh_dist[,10]), FUN = length)
                
                #Count the number of ROHs and sum Length for each catgory
                for(roh_class in 1:6) {
                        dta_roh[roh_class,6] = sum(x[3][x[2] == roh_class])/(sort(unique(dtaroh_dist[,2]))[length(unique(dtaroh_dist[,2]))-1])
                        dta_roh[roh_class,5] = sum(y[3][y[2] == roh_class])/(sort(unique(dtaroh_dist[,2]))[length(unique(dtaroh_dist[,2]))-1])
                
                }

                #write the dataframe
                write.table(dta_roh, file = "./WGS_ROHs_PLINK/WGS_PLINK_${SimID}_ROHsDist.txt", quote = F, row.names = F, col.names = F)

EOF

}

export -f WGS_ROHsCal_PLINK

#Run this function on all Simulation replicates VCFs
parallel -j10 ::: WGS_ROHsCal_PLINK ::: $(find ./Final_VCF -name "FINAL*.vcf.gz")

## Then just merge the SIM replicates into one dataframe

## FHOM
#Create header
echo -e "SimID\tFID\tIID\tO(HOM)_wgs\tE(HOM)_wgs\tN(NM)_wgs\tFHOM_wgs" > ./WGS_FHOM/WGS_FHOM.txt
#Add other files (and SimID as first column)
for file in $(find ./WGS_FHOM -name "HET_WGS_*.het"); do SimID=$(basename -s ".het" ${file} | cut -d '_' -f3); awk -v simid=${SimID} 'NR > 1 {print simid"\t"$0}' ${file} >> ./WGS_FHOM/WGS_FHOM.txt; done
#Rm per SIMID file
rm ./WGS_FHOM/WGS_*_HET.het

#Open R for merging
R --vanilla <<EOF

        #Read RAD-seq + WGS data
        dtasub = read.table("./WGS_FHOM/WGS_FHOM.txt", header = T)
        colnames(dtasub)[3] = "individuals"
        dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
        dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

        #merge both (without windows nb and REP for WGS)
        dtaFULL = merge(dtasub[,c(1,3:7)], dtaTRUE100, by = c("SimID","individuals"))
        dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("SimID","individuals"))
        #Save the data
        write.table(dtaFULL.1, "./Analyses/WGS_PLINK_FHOM.txt", quote = F, col.names = T, row.names = F)

EOF

## FROH
#Create header
echo -e "SimID\tindividuals\tFroh_WGS\tNROH_WGS\tSROH_WGS" > ./WGS_ROHs_PLINK/WGS_PLINK_FROH.txt
#Add all files
for file in $(find ./WGS_ROHs_PLINK -name "WGS_PLINK_*_FROH.txt"); do cat ${file} >> ./WGS_ROHs_PLINK/WGS_PLINK_FROH.txt; done
#Rm per SimID file
rm ./WGS_ROHs_PLINK/WGS_PLINK_*_FROH.txt

#Open R for merging
R --vanilla <<EOF

        #Read RAD-seq + WGS data
        dtasub = read.table("./WGS_ROHs_PLINK/WGS_PLINK_FROH.txt", header = T)
        dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
        dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

        #merge both (without windows nb and REP for WGS)
        dtaFULL = merge(dtasub, dtaTRUE100, by = c("SimID","individuals"))
        dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("SimID","individuals"))
        #Save the data
        write.table(dtaFULL.1, "./Analyses/WGS_PLINK_FROH.txt", quote = F, col.names = T, row.names = F)

EOF

## ROHs Dist
#Create header
echo -e "SimID\tWINDOWS.NB\tREP\tROHs_CLASS\tMean_NROH_per_ind\tMean_Tot_Kb_per_ind" > ./Analyses/WGS_PLINK_ROHsDistribution.txt
#Add all files
for file in $(find ./WGS_ROHs_PLINK/ -name "WGS_PLINK_*_ROHsDist.txt"); do cat ${file} >> ./Analyses/WGS_PLINK_ROHsDistribution.txt; done
#Rm per SimID file
rm ./WGS_ROHs_PLINK/WGS_PLINK_*_ROHsDist.txt


