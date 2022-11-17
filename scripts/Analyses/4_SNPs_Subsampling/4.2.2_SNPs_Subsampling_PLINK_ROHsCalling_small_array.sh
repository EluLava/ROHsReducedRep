#!/bin/bash

#define function for SNPs subsampling and ROHs calling with PLINK
SmallArray_ROHsCalling_PLINK() {

    #SNP_CHR_POS_${SimID}_${REP}.csv from script 3.2.1 is the first argument passed to the function
    FILE=$1

    #Extract SimuID and REP from file name
    SimID=$(basename -s ".csv" ${FILE} | cut -d '_' -f4)
    REP=$(basename -s ".csv" ${FILE} | cut -d'_' -f5)

    #If PLINK directory does not exist, create it to store the results of ROHs calling with PLINK
    mkdir -p ./SNP_ARRAYs_PLINK/Small_ARRAY/SNPs_DENSITY

    #Subsample SNPs to create a new VCF
    bcftools view -R ./SNP_ARRAYs_PLINK/Small_ARRAY/SNP_CHR_POS_${SimID}_${REP}.csv -O z -o ./SNP_ARRAYs_PLINK/Small_ARRAY/TMPVCF_SNPs_${SimID}_${REP}.vcf.gz ./Final_VCF/FINAL_${SimID}.vcf.gz

    #Estimate SNP_DENSITY
    vcftools --gzvcf ./SNP_ARRAYs_PLINK/Small_ARRAY/TMPVCF_SNPs_${SimID}_${REP}.vcf.gz --SNPdensity 100000 --out ./SNP_ARRAYs_PLINK/Small_ARRAY/SNPs_DENSITY/DENSITY_${SimID}_${REP}

    #Average SNP densities per windows --> Average SNP density per KB
    SNPdens=$(awk '{ sum += $4 } END { if (NR > 0) print sum / NR }' ./SNP_ARRAYs_PLINK/Small_ARRAY/SNPs_DENSITY/DENSITY_${SimID}_${REP}.snpden)

    #from VCF to BED
    plink --vcf ./SNP_ARRAYs_PLINK/Small_ARRAY/TMPVCF_SNPs_${SimID}_${REP}.vcf.gz --make-bed --chr-set 30 --out ./SNP_ARRAYs_PLINK/Small_ARRAY/BTMP_${SimID}_${REP}

    #If SNP dens (integer format) lower or equal to 30 , simply use 30
    if [[ `printf "%.0f" $(echo "${SNPdens}" | bc)` -le 30 ]]; then

        #Call ROHs
        plink --bfile ./SNP_ARRAYs_PLINK/Small_ARRAY/BTMP_${SimID}_${REP} --homozyg --homozyg-window-snp 30 --homozyg-density 50 \
        --homozyg-snp 30 --homozyg-gap 1000 --homozyg-window-het 0 --homozyg-het 0 --homozyg-kb 100 --chr-set 30 --out ./SNP_ARRAYs_PLINK/Small_ARRAY/PLINK_ROH_${SimID}_${REP}
  
    #else, if greater than 100, use 100
    elif [[ `printf "%.0f" $(echo "${SNPdens}" | bc)` -ge 100 ]]; then

        #Call ROHs
        plink --bfile ./SNP_ARRAYs_PLINK/Small_ARRAY/BTMP_${SimID}_${REP} --homozyg --homozyg-window-snp 100 --homozyg-density 50 \
        --homozyg-snp 100 --homozyg-gap 1000 --homozyg-window-het 0 --homozyg-het 0 --homozyg-kb 100 --chr-set 30 --out ./SNP_ARRAYs_PLINK/Small_ARRAY/PLINK_ROH_${SimID}_${REP}

    #Finally, else (--> if between 10 and 100) use params (--homozyg-window-snp + --homozyg-snp) according to SNP dens
    else

        plink --bfile ./SNP_ARRAYs_PLINK/Small_ARRAY/BTMP_${SimID}_${REP} --homozyg --homozyg-window-snp `printf "%.0f" $(echo "${SNPdens}*100" | bc)` --homozyg-density 50 \
        --homozyg-snp `printf "%.0f" $(echo "${SNPdens}*100" | bc)` --homozyg-gap 1000 --homozyg-window-het 0 --homozyg-het 0 --homozyg-kb 100 --chr-set 30 --out ./SNP_ARRAYs_PLINK/Small_ARRAY/PLINK_ROH_${SimID}_${REP}

    fi

    #Calculate Het
    plink --bfile ./SNP_ARRAYs_PLINK/Small_ARRAY/BTMP_${SimID}_${REP} --het --chr-set 30 --out ./SNP_ARRAYs_PLINK/Small_ARRAY/HET_${SimID}_${REP}
    
    #Remove TMPVCF file
    rm ./SNP_ARRAYs_PLINK/Small_ARRAY/TMPVCF_SNPs_${SimID}_${REP}.vcf.gz
    #Remove BED, BIM; BAM ETc.
    rm ./SNP_ARRAYs_PLINK/Small_ARRAY/BTMP_${SimID}_${REP}.*
    
    #Open R for F and ROHs distributions estimations
    R --vanilla <<EOF

        library(doParallel)
        library(foreach)
        cl <- 5
        registerDoParallel(cl)

        #Read the indv roh plink file
        dtaroh = read.table("./SNP_ARRAYs_PLINK/Small_ARRAY/PLINK_ROH_${SimID}_${REP}.hom.indiv", header = T)
        #Read the het plink file
        dtahet = read.table("./SNP_ARRAYs_PLINK/Small_ARRAY/HET_${SimID}_${REP}.het", header = T)
        #Read the rohs dist file
        dta_dist = read.table("./SNP_ARRAYs_PLINK/Small_ARRAY/PLINK_ROH_${SimID}_${REP}.hom", header = T)[,-c(3,11:13)]

        #Set the genome length IN KB
        genome_length = 3000000
        #Set individuals nb (again het here becasue some indvs might have 0 ROHs with RAD)
        nb_indiv=nrow(dtahet)

        #Set the individuals vector
        individuals = c(1000, seq(1,nb_indiv-1))

        ## FROH

        #Estimate FROH
        Froh = apply(dtaroh, 1, function(x) as.numeric(as.character(x[5]))/genome_length)

        #First output dataframe
        F_sub = cbind(individuals, Froh, dtaroh[,4], dtaroh[,5])
        colnames(F_sub) = c("individuals", "FROH_sub", "NROH_sub", "SROH_sub")

        #Merging SimID, Froh_sub, Number of windows and SNPs & REP
        SUB_dta = as.data.frame(cbind(rep("${SimID}",nrow(F_sub)), F_sub, rep("SMALL_ARRAY",nrow(F_sub)), rep("${REP}",nrow(F_sub))))
        colnames(SUB_dta) = c("SimID", "individuals", "FROH_sub", "NROH_sub", "SROH_sub", "WINDOWS.NB", "REP")

        #Writing the dataframe
        write.table(SUB_dta, "./SNP_ARRAYs_PLINK/Small_ARRAY/FROH_PLINK_SmallArray_${SimID}_${REP}.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

        ## ROHs Dist

        #Create the output dataframe which will contain lengths per ROH category
        dta_roh = as.data.frame(matrix(nrow = 6, ncol = 6))
        colnames(dta_roh) = c("SimuID","WINDOWS.NB","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")

        #Fill SimuID, REP, and WINDOWS.NB
        dta_roh[,1] = rep("${SimID}",nrow(dta_roh))
        dta_roh[,2] = rep("SMALL_ARRAY",nrow(dta_roh))
        dta_roh[,3] = rep("${REP}",nrow(dta_roh))
        dta_roh[,4] = c(1:6)

        if(nrow(dta_dist) == 0){

            #just put 0 in output d
            dta_roh[,5] = 0
            dta_roh[,6] = 0

        } else {

            #Loop through ROHs and Add the roh class
            foreachoutput = foreach(i=1:nrow(dta_dist), .combine=rbind) %dopar% {
                #Create output
                output = as.data.frame(matrix(nrow=1,ncol=2))
                colnames(output) = c("i","ROH_CLASS")
                output[1,1] = i

                #Define ROH class
                if(dta_dist[i,8] < 2000){
                    output[1,2] = 1
                } else if(dta_dist[i,8] >= 2000 & dta_dist[i,8] < 4000) {
                    output[1,2] = 2
                } else if(dta_dist[i,8] >= 4000 & dta_dist[i,8] < 6000) {
                    output[1,2] = 3
                } else if(dta_dist[i,8] >= 6000 & dta_dist[i,8] < 10000) {
                    output[1,2] = 4
                } else if(dta_dist[i,8] >= 10000 & dta_dist[i,8] < 16000) {
                    output[1,2] = 5
                } else if(dta_dist[i,8] >= 16000) {
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
            x = aggregate(dta_dist.2[,8], by = list(INDIV=dta_dist.2[,2], CLASS=dta_dist.2[,10]), FUN = sum)
            y = aggregate(dta_dist.2[,8], by = list(INDIV=dta_dist.2[,2], CLASS=dta_dist.2[,10]), FUN = length)

            #Count the number of ROHs and and sum Length for each catgory
            for(roh_class in 1:6) {
                #if ROHs calss was detected
                if(roh_class %in% unique(x[,2])){
                    dta_roh[roh_class,6] = sum(x[3][x[2] == roh_class])/(length(individuals))
                    dta_roh[roh_class,5] = sum(y[3][y[2] == roh_class])/(length(individuals))
                #else simply 0
                } else {
                    dta_roh[roh_class,6] = 0
                    dta_roh[roh_class,5] = 0
                }
            }
        }

        #write the output
        write.table(dta_roh, file = "./SNP_ARRAYs_PLINK/Small_ARRAY/ROHsDist_PLINK_SmallArray_${SimID}_${REP}_ROHsDist.txt", quote = F, row.names = F, col.names = F)
EOF

}

export -f SmallArray_ROHsCalling_PLINK

parallel -j10 ::: SmallArray_ROHsCalling_PLINK ::: $(find ./SNP_ARRAYs_PLINK/Small_ARRAY -name "SNP_CHR_POS_*.csv")

### Then just fuse all Simulations in one df

#Create df in case does not exists
mkdir -p ./Analyses

## FROH

#Create header
echo -e "SimID\tindividuals\tFROH_sub\tNROH_sub\tSROH_sub\tWINDOWS.NB\tREP" > ./SNP_ARRAYs_PLINK/Small_ARRAY/SmallArray_PLINK_FROH.txt
#Add all files
for file in $(find ./SNP_ARRAYs_PLINK/Small_ARRAY/ -name "FROH_PLINK_SmallArray_*.txt"); do cat ${file} >> ./SNP_ARRAYs_PLINK/Small_ARRAY/SmallArray_PLINK_FROH.txt; done
#Rm per SimID file
rm ./SNP_ARRAYs_PLINK/Small_ARRAY/FROH_PLINK_SmallArray_*.txt

#Merge with WGS (R)

#Open R for merging
R --vanilla <<EOF

    #Read Small array + True IBD data
    dtasub = read.table("./SNP_ARRAYs_PLINK/Small_ARRAY/SmallArray_PLINK_FROH.txt", header = T)
    dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
    dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

    #merge both (without windows nb and REP for WGS)
    dtaFULL = merge(dtasub, dtaTRUE100, by = c("SimID","individuals"))
    dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("SimID","individuals"))

    #Save the data
    write.table(dtaFULL.1, "./Analyses/SmallArray_PLINK_FROH.txt", quote = F, col.names = T, row.names = F)

EOF

## ROHs Dist

#Create header
echo -e "SimuID\tSUB_TECH\tREP\tROHs_CLASS\tMean_NROH_per_ind\tMean_Tot_Kb_per_ind" > ./Analyses/SmallArray_PLINK_ROHsDistributions.txt

#Loop through file to fill the data
for file in $(find ./SNP_ARRAYs_PLINK/Small_ARRAY/ -name "ROHsDist_PLINK_SmallArray_*.txt"); do cat ${file} >> ./Analyses/SmallArray_PLINK_ROHsDistributions.txt; done

#rm the sub files
rm ./SNP_ARRAYs_PLINK/Small_ARRAY/ROHsDist_PLINK_SmallArray_*.txt



