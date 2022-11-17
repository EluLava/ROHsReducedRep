#!/bin/bash

LargeArray_bedfiles_and_TFPN_Rates_PLINK(){

    #First argument passed to the function is ROH file with one line per ROH --> .hom for PLINK
    ROHfile=$1

    #Extract SimuID
    SimID=$(basename -s ".hom" ${ROHfile} | cut -d'_' -f3)
    
    #Extract Replicate
    REP=$(basename -s ".hom" ${ROHfile} | cut -d'_' -f4)

    #Create directory if needed
    mkdir -p ./TFPN_Rates_PLINK/LargeArray/bedfiles
    mkdir -p ./TFPN_Rates_PLINK/LargeArray/Rates

    #Find the corresponding INDV file (for list of indvs)
    INDIV_FILE=$(find ./SNP_ARRAYs_PLINK/Large_ARRAY/ -name "PLINK_ROH_${SimID}_${REP}.hom.indiv")

    #Read the data in R to create one bedfile per individual
    R --vanilla <<EOF

        #Read the file with all indiv in case one has F of 0 --> 0 ROHs
        dta_indiv = read.table("${INDIV_FILE}", header = T)
    
        #Get vector of all individuals from which we can iterate
        individuals = unique(dta_indiv[,2])
    
        #We don't need dta_indiv anymore
        rm(dta_indiv)
    
        #Read the file with on ROH per line
        dta = read.table("${ROHfile}", header = T)
    
        #Iterate through individuals
        for(indiv in individuals){
    
            #Extract ROHs only for this individual and select only CHR, POS1 and POS2 column
            dta_ind = dta[dta[,2] == indiv, c(4,7,8)]
    
            #change colnames
            colnames(dta_ind) = c("#chrom", "start", "stop")

            #PASS ALL TO INTEGER for bedtools
            dta_ind[,1] = as.integer(dta_ind[,1])
            dta_ind[,2] = as.integer(dta_ind[,2])
            dta_ind[,3] = as.integer(dta_ind[,3])

            #Write the bedfile
            write.table(dta_ind, paste0("./TFPN_Rates_PLINK/LargeArray/bedfiles/PLINK_LargeArray_BED_${SimID}_${REP}_Ind_", indiv, ".bed"), row.names = F, quote = F, col.names = T, sep = "\t")
    
        }

EOF
    
    #List all the bedfiles for this replicate
    INDIV_BEDFILES=$(find ./TFPN_Rates_PLINK/LargeArray/bedfiles/ -name "PLINK_LargeArray_BED_${SimID}_${REP}_Ind_*.bed")

    #Create the header of the output file (one line with all rate epr indv)
    echo -e "SimuID\tNB_RAD_FRAG\tREP\tINDIV\tTP\tTN\tFP\tFN" > ./TFPN_Rates_PLINK/LargeArray/Rates/LargeArray_${SimID}_${REP}_TF_PN_Rates_GEN100.txt

    #Loop through these BEDFILES
    for BEDFILE in ${INDIV_BEDFILES}
    do

        #Extract individual
        INDIV=$(basename -s ".bed" ${BEDFILE} | cut -d'_' -f7)
        #Get corresponding WGS file
        TRUEIBD_FILE=./TFPN_Rates_PLINK/TRUE_IBD/100GEN/bedfiles/TRUEIBD_BED_${SimID}_Ind_${INDIV}.bed

        #Estimate the TP true positive rate
        TP=$(bedtools intersect -a ${TRUEIBD_FILE} -b ${BEDFILE} -g ./genome.bed | awk -v TP=0 -v genlength=3000000000 '{TP = TP + ($3 - $2)};END{print TP/genlength}')

        #TN true-negative rate with complement function from bedtools, but only works with one single file --> I need to "merge" WGS and RAd file first, then complement the result
        TN=$(cat ${TRUEIBD_FILE} ${BEDFILE} | sort -k1,1n -k2,2n | bedtools complement -i stdin -g ./genome.bed | awk -v TN=0 -v genlength=3000000000 '{TN = TN + ($3 - $2)};END{print TN/genlength}')
    
        #Estimate the FP rate
        FN=$(bedtools subtract -a ${TRUEIBD_FILE} -b ${BEDFILE} -g ./genome.bed | awk -v FN=0 -v genlength=3000000000 '{FN = FN + ($3 - $2)};END{print FN/genlength}')
    
        #Estimate the TP rate
        FP=$(bedtools subtract -a ${BEDFILE} -b ${TRUEIBD_FILE} -g ./genome.bed | awk -v FP=0 -v genlength=3000000000 '{FP = FP + ($3 - $2)};END{print FP/genlength}')

        #Complete the TFPN Rates FINAL file
        echo -e "${SimID}\tLargeArray\t${REP}\t${INDIV}\t${TP}\t${TN}\t${FP}\t${FN}" >> ./TFPN_Rates_PLINK/LargeArray/Rates/LargeArray_${SimID}_${REP}_TF_PN_Rates_GEN100.txt

        #rm the INDIV bedfile
        rm ${BEDFILE}

    done


}


export -f LargeArray_bedfiles_and_TFPN_Rates_PLINK

parallel -j10 ::: LargeArray_bedfiles_and_TFPN_Rates_PLINK ::: $(find ./SNP_ARRAYs_PLINK/Large_ARRAY/ -name "PLINK_ROH_*.hom")

#Then just merge all the replicate into one df

#List the files
FILES=./TFPN_Rates_PLINK/LargeArray/Rates/LargeArray_*_TF_PN_Rates_GEN100.txt
#Create header
echo -e "SimID\tNB_RAD_FRAG\tREP\tINDIV\tTP\tTN\tFP\tFN" > ./Analyses/LargeArray_PLINK_TFPN_Rates_GEN100.txt
#Add the files
for f in $FILES; do awk 'NR > 1 {print $0}' $f >> ./Analyses/LargeArray_PLINK_TFPN_Rates_GEN100.txt; done


