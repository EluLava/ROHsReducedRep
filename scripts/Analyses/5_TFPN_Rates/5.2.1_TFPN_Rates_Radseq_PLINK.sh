#!/bin/bash

Radseq_bedfiles_and_TFPN_Rates_PLINK(){

    #First argument passed to the function is ROH file with one line per ROH --> .hom for PLINK
    ROHfile=$1

    #Extract SimuID
    SimID=$(basename -s ".hom" ${ROHfile} | cut -d'_' -f3)
    
    #Extract Replicate
    REP=$(basename -s ".hom" ${ROHfile} | cut -d'_' -f5)
    
    #Extract NB_WINDOWS
    RAD_WINDOW_NB=$(basename -s ".hom" ${ROHfile} | cut -d'_' -f4)

    #Create directory if needed
    mkdir -p ./TFPN_Rates_PLINK/RADseq/bedfiles
    mkdir -p ./TFPN_Rates_PLINK/RADseq/Rates/

    #Find the corresponding INDV file (for list of indvs)
    INDIV_FILE=$(find ./RADseq_PLINK/ -name "PLINK_ROH_${SimID}_${RAD_WINDOW_NB}_${REP}.hom.indiv")

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
            write.table(dta_ind, paste0("./TFPN_Rates_PLINK/RADseq/bedfiles/PLINK_Radseq_BED_${SimID}_${RAD_WINDOW_NB}_${REP}_Ind_", indiv, ".bed"), row.names = F, quote = F, col.names = T, sep = "\t")
    
        }

EOF
    
    #List all the bedfiles for this replicate
    INDIV_BEDFILES=$(find ./TFPN_Rates_PLINK/RADseq/bedfiles/ -name "PLINK_Radseq_BED_${SimID}_${RAD_WINDOW_NB}_${REP}_Ind_*.bed")

    #Create the header of the output file (one line with all rate epr indv) 100GEN
    echo -e "SimuID\tNB_RAD_FRAG\tREP\tINDIV\tTP\tTN\tFP\tFN" > ./TFPN_Rates_PLINK/RADseq/Rates/RAD_${SimID}_${RAD_WINDOW_NB}_${REP}_TF_PN_Rates_GEN100.txt

    #Create the header of the output file (one line with all rate epr indv) 1000GEN
    echo -e "SimuID\tNB_RAD_FRAG\tREP\tINDIV\tTP\tTN\tFP\tFN" > ./TFPN_Rates_PLINK/RADseq/Rates/RAD_${SimID}_${RAD_WINDOW_NB}_${REP}_TF_PN_Rates_GEN1000.txt

    #Loop through these BEDFILES
    for BEDFILE in ${INDIV_BEDFILES}
    do

        #Extract individual
        INDIV=$(basename -s ".bed" ${BEDFILE} | cut -d'_' -f8)

        ## 100 GEN

        #Get corresponding TRUE IBD file
        TRUE_IBD=./TFPN_Rates_PLINK/TRUE_IBD/100GEN/bedfiles/TRUEIBD_BED_${SimID}_Ind_${INDIV}.bed

        #Estimate the TP true positive rate
        TP=$(bedtools intersect -a ${TRUE_IBD} -b ${BEDFILE} -g ./genome.bed | awk -v TP=0 -v genlength=3000000000 '{TP = TP + ($3 - $2)};END{print TP/genlength}')

        #TN true-negative rate with complement function from bedtools, but only works with one single file --> I need to "merge" WGS and RAd file first, then complement the result
        TN=$(cat ${TRUE_IBD} ${BEDFILE} | sort -k1,1n -k2,2n | bedtools complement -i stdin -g ./genome.bed | awk -v TN=0 -v genlength=3000000000 '{TN = TN + ($3 - $2)};END{print TN/genlength}')
    
        #Estimate the FP rate
        FN=$(bedtools subtract -a ${TRUE_IBD} -b ${BEDFILE} -g ./genome.bed | awk -v FN=0 -v genlength=3000000000 '{FN = FN + ($3 - $2)};END{print FN/genlength}')
    
        #Estimate the TP rate
        FP=$(bedtools subtract -a ${BEDFILE} -b ${TRUE_IBD} -g ./genome.bed | awk -v FP=0 -v genlength=3000000000 '{FP = FP + ($3 - $2)};END{print FP/genlength}')

        #Complete the TFPN Rates FINAL file
        echo -e "${SimID}\t${RAD_WINDOW_NB}\t${REP}\t${INDIV}\t${TP}\t${TN}\t${FP}\t${FN}" >> ./TFPN_Rates_PLINK/RADseq/Rates/RAD_${SimID}_${RAD_WINDOW_NB}_${REP}_TF_PN_Rates_GEN100.txt

        ## 1000 GEN

        #Get corresponding TRUE IBD file
        TRUE_IBD=./TFPN_Rates_PLINK/TRUE_IBD/1000GEN/bedfiles/TRUEIBD_BED_${SimID}_Ind_${INDIV}.bed

        #Estimate the TP true positive rate
        TP=$(bedtools intersect -a ${TRUE_IBD} -b ${BEDFILE} -g ./genome.bed | awk -v TP=0 -v genlength=3000000000 '{TP = TP + ($3 - $2)};END{print TP/genlength}')

        #TN true-negative rate with complement function from bedtools, but only works with one single file --> I need to "merge" WGS and RAd file first, then complement the result
        TN=$(cat ${TRUE_IBD} ${BEDFILE} | sort -k1,1n -k2,2n | bedtools complement -i stdin -g ./genome.bed | awk -v TN=0 -v genlength=3000000000 '{TN = TN + ($3 - $2)};END{print TN/genlength}')

        #Estimate the FP rate
        FN=$(bedtools subtract -a ${TRUE_IBD} -b ${BEDFILE} -g ./genome.bed | awk -v FN=0 -v genlength=3000000000 '{FN = FN + ($3 - $2)};END{print FN/genlength}')

        #Estimate the TP rate
        FP=$(bedtools subtract -a ${BEDFILE} -b ${TRUE_IBD} -g ./genome.bed | awk -v FP=0 -v genlength=3000000000 '{FP = FP + ($3 - $2)};END{print FP/genlength}')

        #Complete the TFPN Rates FINAL file
        echo -e "${SimID}\t${RAD_WINDOW_NB}\t${REP}\t${INDIV}\t${TP}\t${TN}\t${FP}\t${FN}" >> ./TFPN_Rates_PLINK/RADseq/Rates/RAD_${SimID}_${RAD_WINDOW_NB}_${REP}_TF_PN_Rates_GEN1000.txt

        #rm the INDIV bedfile
        rm ${BEDFILE}

    done


}


export -f Radseq_bedfiles_and_TFPN_Rates_PLINK

parallel -j20 ::: Radseq_bedfiles_and_TFPN_Rates_PLINK ::: $(find ./RADseq_PLINK/ -name "PLINK_ROH_*.hom")

#Then just merge all the replicate into one df

## 100 GEN

#List the files
FILES=./TFPN_Rates_PLINK/RADseq/Rates/RAD_*_TF_PN_Rates_GEN100.txt
#Create header
echo -e "SimID\tNB_RAD_FRAG\tREP\tINDIV\tTP\tTN\tFP\tFN" > ./Analyses/RADSEQ_PLINK_TFPN_Rates_GEN100.txt
#Add the files
for f in $FILES; do awk 'NR > 1 {print $0}' $f >> ./Analyses/RADSEQ_PLINK_TFPN_Rates_GEN100.txt; done

## 1000 GEN

#List the files
FILES=./TFPN_Rates_PLINK/RADseq/Rates/RAD_*_TF_PN_Rates_GEN1000.txt
#Create header
echo -e "SimID\tNB_RAD_FRAG\tREP\tINDIV\tTP\tTN\tFP\tFN" > ./Analyses/RADSEQ_PLINK_TFPN_Rates_GEN1000.txt
#Add the files
for f in $FILES; do awk 'NR > 1 {print $0}' $f >> ./Analyses/RADSEQ_PLINK_TFPN_Rates_GEN1000.txt; done


