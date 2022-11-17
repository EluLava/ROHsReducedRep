#!/bin/bash

#Create outut directory SAME FOR 1000 GEN jsust change name
mkdir -p ./TFPN_Rates_RZooRoH/TRUE_IBD/100GEN/bedfiles
mkdir -p ./TFPN_Rates_RZooRoH/TRUE_IBD/100GEN/Rates
mkdir -p ./TFPN_Rates_RZooRoH/WGS/bedfiles

#Create bedfile with regions of the genome within ROHs for WGS
WGS_bedfiles(){

    #First argument passed to the function is VCF file
    VCF=$1

    #From VCF extract SimID
    SimID=$(basename -s ".vcf.gz" ${VCF} | cut -d '_' -f2)

    R --vanilla <<EOF

    #Read the file with on ROH per line
    dta = read.table("./WGS_ROHs_RZooRoH/ROHs_WGS_${SimID}.hom", header = T)

    #Extract indv
    individuals = unique(dta[,1])

    #Iterate through individuals
    for(indiv in individuals){
    
        #Extract ROHs only for this individual and select only CHR, POS1 and POS2 column
        dta_ind = dta[dta[,1] == indiv, c(2,5,6)]
        #Order by CHR first then start positionx2
        dta_ind.2 = dta_ind[order(dta_ind[,1],dta_ind[,2]),]
        
        #change colnames
        colnames(dta_ind.2) = c("#chrom", "start", "stop")
    
        #PASS ALL TO INTEGER for bedtools
        dta_ind.2[,1] = as.integer(dta_ind.2[,1])
        dta_ind.2[,2] = as.integer(dta_ind.2[,2])
        dta_ind.2[,3] = as.integer(dta_ind.2[,3])
    
        #Write the bedfile
        write.table(dta_ind.2, paste0("./TFPN_Rates_RZooRoH/WGS/bedfiles/WGS_BED_${SimID}_Ind_", indiv, ".bed"), row.names = F, quote = F, col.names = T, sep = "\t")
    }

EOF

}

export -f WGS_bedfiles

#Apply the function to all simulation files
parallel -j10 ::: WGS_bedfiles ::: $(find ./Final_VCF/ -name "FINAL_*.vcf.gz")


###########################################################################
# TRUE IBD BEDFILES
#Create bedfile with regions of the genome which coalesce less than 100 (or 1000) GEN ago

#List SimIDs
SimuIDs=$(find ./Final_VCF/ -name "FINAL_*.vcf.gz" | xargs -I {} bash -c 'basename -s ".vcf.gz" {} | cut -d'_' -f2')

#Loop trhough simualtions
for sim in ${SimuIDs}
do

R --vanilla <<EOF
    dta = read.table("./TRUE_IBD/100GEN/stats/TRUE_IBD_SEGMENTS_${sim}_GEN100.txt", header = T)
    
    #Subset data for this SimID
    dtaSIM = dta[dta[,1] == "${sim}",]
    
    #List individuals
    individuals = unique(dtaSIM[,2])
    
    #Iterate through individuals
    for(indiv in individuals){

        #Extract ROHs only for this individual and select only CHR, POS1 and POS2 column
        dta_ind = dtaSIM[dtaSIM[,2] == indiv, c(3,4,5)]
        
        #change colnames
        colnames(dta_ind) = c("#chrom", "start", "stop")

        #PASS ALL TO INTEGER for bedtools
        dta_ind[,1] = as.integer(dta_ind[,1])
        dta_ind[,2] = as.integer(dta_ind[,2])
        dta_ind[,3] = as.integer(dta_ind[,3])
        
        #Write the bedfile
        write.table(dta_ind, paste0("./TFPN_Rates_RZooRoH/TRUE_IBD/100GEN/bedfiles/TRUEIBD_BED_${sim}_Ind_", indiv, ".bed"), row.names = F, quote = F, col.names = T, sep = "\t")
    }

EOF

done

###########################################################################

# RATES

#Create the function
INDIV_Rates_TRUE_IBD() {

#First argument passed to the function = TRUE IBD bedfile !
BEDFILE=$1

#Extract Simulation
SimuID=$(basename -s ".bed" ${BEDFILE} | cut -d'_' -f3)

INDIV=$(basename -s ".bed" ${BEDFILE} | cut -d'_' -f5)

#Get corresponding WGS RZooRoH file
WGS_FILE=./TFPN_Rates_RZooRoH/WGS/bedfiles/WGS_BED_${SimuID}_Ind_${INDIV}.bed

#Estimate the TP rate
TP=$(bedtools intersect -a ${BEDFILE} -b ${WGS_FILE} -g ./genome.bed | awk -v TP=0 -v genlength=3000000000 '{TP = TP + ($3 - $2)};END{print TP/genlength}')

#TN rate with complement function from bedtools, but only works with one single file --> I need to "merge" WGS and RAd file first, then complement the result
TN=$(cat ${BEDFILE} ${WGS_FILE} | sort -k1,1n -k2,2n | bedtools complement -i stdin -g ./genome.bed | awk -v TN=0 -v genlength=3000000000 '{TN = TN + ($3 - $2)};END{print TN/genlength}')

#Estimate the FP rate
FN=$(bedtools subtract -a ${BEDFILE} -b ${WGS_FILE} -g ./genome.bed | awk -v FN=0 -v genlength=3000000000 '{FN = FN + ($3 - $2)};END{print FN/genlength}')

#Estimate the TP rate
FP=$(bedtools subtract -a ${WGS_FILE} -b ${BEDFILE} -g ./genome.bed | awk -v FP=0 -v genlength=3000000000 '{FP = FP + ($3 - $2)};END{print FP/genlength}')

#Complete the TFPN Rates FINAL file
echo -e "${SimuID}\tRZooRoH\t${INDIV}\t${TP}\t${TN}\t${FP}\t${FN}" >> ./Analyses/WGS_RZooRoH_TFPN_Rates_GEN100.txt

}

export -f INDIV_Rates_TRUE_IBD

#Create the header for the file where we'll record TP, FP, TN, FN rates per indiv
echo -e "SimuID\tNB_RAD_FRAG\tINDIV\tTP\tTN\tFP\tFN" > ./Analyses/WGS_RZooRoH_TFPN_Rates_GEN100.txt
#Apply function to all BEDFILES
find ./TFPN_Rates_RZooRoH/TRUE_IBD/100GEN/bedfiles/ -name "TRUEIBD_BED_*_Ind_*.bed" | xargs -I {} bash -c 'INDIV_Rates_TRUE_IBD {}'





