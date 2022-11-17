#!/bin/bash

#create the directory
mkdir -p ./SNP_ARRAYs/Large_ARRAY/

#Check if filtering on MAF etc. was already done with the small array ! if yes no need, if no, do it --> Check for the tmp file created in 3.2.1
if [! -f ./SNP_ARRAYs/tmp_all_frq_estimated]; then
    #Estimate Allelic frequencies
    #Define the function
    PRESUBSAMPLING_ARRAY() {
    
        #First argument passed to the function is VCF WGS file
        VCF=$1
        #Extract the simulation ID
        SimID=$(basename -s ".vcf.gz" ${VCF} | cut -d'_' -f2)
        
        #Create directory if does not exist
        mkdir -p ./SNP_ARRAYs/Small_ARRAY
        
        #Filter on MAF
        bcftools view -q 0.05:minor -O z -o ./Final_VCF/MAF_${SimID}.vcf.gz ${VCF}
        
        #Create a dataframe with only CHR & POS
        bcftools query --format '%CHROM  %POS\n' ./Final_VCF/MAF_${SimID}.vcf.gz > ./SNP_ARRAYs/SNP_info_maf0.05_${SimID}.txt
        
        #Estimate alleles frequencies for this VCF
        vcftools --gzvcf ./Final_VCF/MAF_${SimID}.vcf.gz --freq2 --out ./SNP_ARRAYs/All_Freq_EST_maf0.05_${SimID}
        
        #rm the filtered VCF
        rm ./Final_VCF/MAF_${SimID}.vcf.gz
        
        #Open R to merge ALL-freq & SNP_info & add MAF column
        R --vanilla <<EOF
        
            #read the SNP INFO
            dta_snp_info = read.table("./SNP_ARRAYs/SNP_info_maf0.05_${SimID}.txt", header = T)
            #Read the All freq
            dta_all_frq = read.table("./SNP_ARRAYs/All_Freq_EST_maf0.05_${SimID}.frq", header = F, skip = 1)[,c(1,2,5,6)]
            
            #Update colnames of dta_frq for merging
            colnames(dta_all_frq) = c("CHR", "POS", "ALL_Frq_A", "ALL_Frq_T")
            colnames(dta_snp_info) = c("CHR", "POS")
            
            #Create new column for MAF
            dta_all_frq = cbind(dta_all_frq, MAF=vector(length = nrow(dta_all_frq)))
            
            #Get the minimum all.frq for each SNP in the MAF column
            dta_all_frq[,5] = apply(dta_all_frq, 1, FUN = function(x) {min(c(x[3],x[4]))})
            
            #merge both df
            dta_full = merge(dta_snp_info, dta_all_frq[,c(1,2,5)], by = c("CHR", "POS"))
            
            #write the df
            write.table(dta_full, "./SNP_ARRAYs/SNP_info_MAF_${SimID}.txt", row.names = F, col.names = T, quote = F, sep = '\t')
        
EOF
        
        #rm SNP-info file
        rm ./SNP_ARRAYs/SNP_info_maf0.05_${SimID}.txt

        #Create a tmp file so I know this script was run and won't need to bererun for the large array !
        echo "DONE" > ./SNP_ARRAYs/tmp_all_frq_estimated
    
    }
    
    export -f PRESUBSAMPLING_ARRAY
    
    #Launch on VCF files
    parallel -j10 ::: PRESUBSAMPLING_ARRAY ::: $(find ./Final_VCF/ -name "FINAL_*.vcf.gz")


fi


########################################################################################################################################################################
######################################################################## SELECTING THE SNPs  ###########################################################################
########################################################################################################################################################################


#Create a second function will will handle the python code
SUBSAMPLING_ARRAY() {

        file=$1
        SimuID=$(basename -s ".vcf.gz" ${file} | cut -d'_' -f2)

        python3 - << EOF

import numpy as np
import pandas as pd
import random as rd
import os
import sys

#define a function
def SNP_Arraying(SimID):

    #Since we take SNP with higher MAF, we only make 1 replicate !
    REP=1

    #Create the start and stop positions of the windows
    #start = 1, stop = 2999997001, num = 1000000
    start_windows=np.linspace(start = 1, stop = 2999997001, num = 1000000)
    stop_windows=np.linspace(start = 2501, stop = 2999999501, num = 1000000)

    #Read the SNP info file
    snp_info = pd.read_csv('./SNP_ARRAYs/SNP_info_MAF_{}.txt'.format(SimID), header='infer', sep='\t')

    #Add a EXT_POS column: CHR*POS to be able to loop through genome as one big CHR
    snp_info = snp_info.assign(EXT_POS = 'NaN')

    #Fill the new column
    snp_info['EXT_POS'] = ((snp_info['CHR'] -  1)*100000000) + snp_info['POS']

    #Create names for my df with CHR & POS
    snp_pos_names = ['CHR', 'POS']
    #Create an empty panda.df with CHR & SNP_positions that we will fill while looping through windows
    snp_pos = pd.DataFrame(columns = snp_pos_names)

    #Loop through the windows
    for wind in range(0,(len(start_windows)-1)):

        #subsample the df
        sub = snp_info[(snp_info['EXT_POS'] > start_windows[wind]) & (snp_info['EXT_POS'] < stop_windows[wind])]

        #Check if the subsample dataframe is emtpy and go to next iteartion if it is
        if sub.empty:
            continue

        #If the df is NOT emtpy
        else:
            #Sample the row with higher MAF in sub
            selected_row = sub.loc[sub['MAF'].idxmax()]
            snp_pos = snp_pos.append(pd.Series([selected_row.loc['CHR'], selected_row.loc['POS']], index = snp_pos.columns), ignore_index = True)
    
    #Pass CHR & POS to integer
    snp_pos['CHR'] = snp_pos['CHR'].astype(int)
    snp_pos['POS'] = snp_pos['POS'].astype(int)
    #add # before CHR to mimick bedfile
    snp_pos.rename(columns = {'CHR':'#CHR'}, inplace = True)

    #Save the subsampled CHR and POS in a file
    snp_pos.to_csv('./SNP_ARRAYs/Large_ARRAY/SNP_CHR_POS_{}_{}.csv'.format(SimID, REP), index = False, header = False, sep = "\t")

#LAUNCH THE FUNCTION WITH SimID
SNP_Arraying('${SimuID}')

EOF

}

export -f SUBSAMPLING_ARRAY

#List the files
FILES=$(find ./Final_VCF/ -name "FINAL_*.vcf.gz")

#Loop trough VCF files
for f in ${FILES}; do SUBSAMPLING_ARRAY ${f}; done
