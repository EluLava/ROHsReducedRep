#!/bin/bash

SmallArray_bedfiles_and_TFPN_Rates_PLINK(){

    #First argument passed to the function is ROH file with one line per ROH --> .hom for PLINK
    ROHfile=$1

    #Extract SimuID
    SimID=$(basename -s ".hom" ${ROHfile} | cut -d'_' -f3)
    
    #Extract Replicate
    REP=$(basename -s ".hom" ${ROHfile} | cut -d'_' -f4)
    
    #List all the bedfiles for this replicate
    INDIV_BEDFILES=$(find ./TFPN_Rates_PLINK/SmallArray/bedfiles/ -name "PLINK_SmallArray_BED_${SimID}_${REP}_Ind_*.bed")

    #Create the header of the output file (one line with all rate epr indv)
    echo -e "SimuID\tNB_RAD_FRAG\tREP\tINDIV\tTP\tTN\tFP\tFN" > ./TFPN_Rates_PLINK/SmallArray/Rates/SmallArray_${SimID}_${REP}_TF_PN_Rates_GEN1000.txt

    #Loop through these BEDFILES
    for BEDFILE in ${INDIV_BEDFILES}
    do

        #Extract individual
        INDIV=$(basename -s ".bed" ${BEDFILE} | cut -d'_' -f7)
        #Get corresponding WGS file
        TRUEIBD_FILE=./TFPN_Rates_PLINK/TRUE_IBD/1000GEN/bedfiles/TRUEIBD_BED_${SimID}_Ind_${INDIV}.bed

        #Estimate the TP true positive rate
        TP=$(bedtools intersect -a ${TRUEIBD_FILE} -b ${BEDFILE} -g ./genome.bed | awk -v TP=0 -v genlength=3000000000 '{TP = TP + ($3 - $2)};END{print TP/genlength}')

        #TN true-negative rate with complement function from bedtools, but only works with one single file --> I need to "merge" WGS and RAd file first, then complement the result
        TN=$(cat ${TRUEIBD_FILE} ${BEDFILE} | sort -k1,1n -k2,2n | bedtools complement -i stdin -g ./genome.bed | awk -v TN=0 -v genlength=3000000000 '{TN = TN + ($3 - $2)};END{print TN/genlength}')
    
        #Estimate the FP rate
        FN=$(bedtools subtract -a ${TRUEIBD_FILE} -b ${BEDFILE} -g ./genome.bed | awk -v FN=0 -v genlength=3000000000 '{FN = FN + ($3 - $2)};END{print FN/genlength}')
    
        #Estimate the TP rate
        FP=$(bedtools subtract -a ${BEDFILE} -b ${TRUEIBD_FILE} -g ./genome.bed | awk -v FP=0 -v genlength=3000000000 '{FP = FP + ($3 - $2)};END{print FP/genlength}')

        #Complete the TFPN Rates FINAL file
        echo -e "${SimID}\tSmallArray\t${REP}\t${INDIV}\t${TP}\t${TN}\t${FP}\t${FN}" >> ./TFPN_Rates_PLINK/SmallArray/Rates/SmallArray_${SimID}_${REP}_TF_PN_Rates_GEN1000.txt

        #rm the INDIV bedfile
        rm ${BEDFILE}

    done


}


export -f SmallArray_bedfiles_and_TFPN_Rates_PLINK

parallel -j1 ::: SmallArray_bedfiles_and_TFPN_Rates_PLINK ::: $(find ./SNP_ARRAYs_PLINK/Small_ARRAY/ -name "PLINK_ROH_*.hom")


parallel -j10 ::: SmallArray_bedfiles_and_TFPN_Rates_PLINK ::: $(find ./SNP_ARRAYs_PLINK/Small_ARRAY/ -name "PLINK_ROH_*.hom")

#Then just merge all the replicate into one df

#List the files
FILES=./TFPN_Rates_PLINK/SmallArray/Rates/SmallArray_*_TF_PN_Rates_GEN1000.txt
#Create header
echo -e "SimID\tNB_RAD_FRAG\tREP\tINDIV\tTP\tTN\tFP\tFN" > ./Analyses/SmallArray_PLINK_TFPN_Rates_GEN1000.txt
#Add the files
for f in $FILES; do awk 'NR > 1 {print $0}' $f >> ./Analyses/SmallArray_PLINK_TFPN_Rates_GEN1000.txt; done


