#!/bin/bash

#Define function
SNPDENSITY_RAD_Subsampling_ROH_calling_PLINK() {

    #First argument passed to the function is the WGS VCF
    VCF=$1
    #Extract SimuID from file name
    SimID=$(basename -s ".vcf.gz" ${VCF} | cut -d '_' -f2)

    #Second argument: # of RAD windows we'd like to subsample (≠ percentage of genome subsampled !)
    RAD_WINDOW_NB=$2
    
    #Third argument is replicate
    REP=$3

    #If PLINK directory does not exist, create it to store the results of ROHs calling with PLINK
    mkdir -p ./SUPP_ANALYSES/SNPdensity_PLINK/bedfiles
    mkdir -p ./SUPP_ANALYSES/SNPdensity_PLINK/SNPs_DENSITY
    mkdir -p ./SUPP_ANALYSES/SNPdensity_PLINK/FROH
    
    #Randomly create BED file with RAD-like fragments (500bp)
    bedtools random -l 500 -n ${RAD_WINDOW_NB} -g ./genome.bed | sort -k1V,1 -k2n,2 > ./SUPP_ANALYSES/SNPdensity_PLINK/bedfiles/sorted_${SimID}_${RAD_WINDOW_NB}_${REP}.bed
    
    #Subset SNPs
    vcftools --gzvcf ${VCF} --bed ./SUPP_ANALYSES/SNPdensity_PLINK/bedfiles/sorted_${SimID}_${RAD_WINDOW_NB}_${REP}.bed --out ./SUPP_ANALYSES/SNPdensity_PLINK/TMP_${SimID}_${RAD_WINDOW_NB}_${REP} --recode
    
    #Estimate SNP_DENSITY
    vcftools --vcf ./SUPP_ANALYSES/SNPdensity_PLINK/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}.recode.vcf --SNPdensity 100000 --out ./SUPP_ANALYSES/SNPdensity_PLINK/SNPs_DENSITY/DENSITY_${SimID}_${RAD_WINDOW_NB}_${REP}

    #Average SNP densities per windows --> Average SNP density per KB
    SNPdens=$(awk '{ sum += $4 } END { if (NR > 0) print sum / NR }' ./SUPP_ANALYSES/SNPdensity_PLINK/SNPs_DENSITY/DENSITY_${SimID}_${RAD_WINDOW_NB}_${REP}.snpden)

    #From VCF to bed
    plink --vcf ./SUPP_ANALYSES/SNPdensity_PLINK/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}.recode.vcf --make-bed --chr-set 30 --out ./SUPP_ANALYSES/SNPdensity_PLINK/BTMP_${SimID}_${RAD_WINDOW_NB}_${REP}

    #If SNP dens (integer format) lower or equal to 30 , simply use 30
    if [[ `printf "%.0f" $(echo "${SNPdens}" | bc)` -le 30 ]]; then

        #Call ROHs
        plink --bfile ./SUPP_ANALYSES/SNPdensity_PLINK/BTMP_${SimID}_${RAD_WINDOW_NB}_${REP} --homozyg --homozyg-window-snp 30 --homozyg-density 50 \
        --homozyg-snp 30 --homozyg-gap 1000 --homozyg-window-het 1 --homozyg-het 2 --homozyg-kb 100 --chr-set 30 --out ./SUPP_ANALYSES/SNPdensity_PLINK/PLINK_ROH_${SimID}_${RAD_WINDOW_NB}_${REP}
  
    #else, if greater or equal to 100
    elif [[ `printf "%.0f" $(echo "${SNPdens}" | bc)` -ge 100 ]]; then

        #Call ROHs
        plink --bfile ./SUPP_ANALYSES/SNPdensity_PLINK/BTMP_${SimID}_${RAD_WINDOW_NB}_${REP} --homozyg --homozyg-window-snp 100 --homozyg-density 100 \
        --homozyg-snp 100 --homozyg-gap 1000 --homozyg-window-het 1 --homozyg-het 2 --homozyg-kb 100 --chr-set 30 --out ./SUPP_ANALYSES/SNPdensity_PLINK/PLINK_ROH_${SimID}_${RAD_WINDOW_NB}_${REP}

    #Finally, else (--> if between 30 and 100) use params (--homozyg-window-snp + --homozyg-snp) according to SNP dens
    else

        plink --bfile ./SUPP_ANALYSES/SNPdensity_PLINK/BTMP_${SimID}_${RAD_WINDOW_NB}_${REP} --homozyg --homozyg-window-snp `printf "%.0f" $(echo "${SNPdens}*100" | bc)` --homozyg-density 50 \
        --homozyg-snp `printf "%.0f" $(echo "${SNPdens}*100" | bc)` --homozyg-gap 1000 --homozyg-window-het 1 --homozyg-het 2 --homozyg-kb 100 --chr-set 30 --out ./SUPP_ANALYSES/SNPdensity_PLINK/PLINK_ROH_${SimID}_${RAD_WINDOW_NB}_${REP}


    fi

    #Remove TMPVCF file
    rm ./SUPP_ANALYSES/SNPdensity_PLINK/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}.recode.vcf
    #Remove BED, BIM; BAM ETc.
    rm ./SUPP_ANALYSES/SNPdensity_PLINK/BTMP_${SimID}_${RAD_WINDOW_NB}_${REP}.*
    #Rm bedfile as well
    rm ./SUPP_ANALYSES/SNPdensity_PLINK/bedfiles/sorted_${SimID}_${RAD_WINDOW_NB}_${REP}.bed
    
    #Open R for FROH estimation and comparison with WGS estimations
    R --vanilla <<EOF

        #Read the indv roh plink file
        dtaroh = read.table("./SUPP_ANALYSES/SNPdensity_PLINK/PLINK_ROH_${SimID}_${RAD_WINDOW_NB}_${REP}.hom.indiv", header = T)

        #read the dta TRUE IBD
        dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
        #read the dta TRUE IBD
        dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

        #Set the genome length IN KB
        genome_length = 3000000
        #Set individuals nb (TRUE IBD here because some indvs might have 0 ROHs with RAD)
        nb_indiv=length(unique((dtaTRUE100[dtaTRUE100[,1] == "${SimID}",2])))

        #Set the individuals vector
        individuals = dtaroh[,2]

        ### FROH

        #Estimate FROH
        Froh = apply(dtaroh, 1, function(x) as.numeric(as.character(x[5]))/genome_length)

        #First output dataframe
        F_sub = as.data.frame(cbind(individuals, as.numeric(as.character(Froh))))
        colnames(F_sub) = c("individuals", "FROH_sub")

        #Merging SimID, Froh_sub, Number of windows & REP
        SUB_dta = as.data.frame(cbind(rep("${SimID}",nrow(F_sub)), F_sub, rep("${RAD_WINDOW_NB}",nrow(F_sub)), rep("${REP}",nrow(F_sub))))
        colnames(SUB_dta) = c("SimID", "individuals", "FROH_sub", "WINDOWS.NB", "REP")

        ### MERGING with TRUE IBD

        #merge by individuals
        dtaFULL = merge(SUB_dta, dtaTRUE100, by = c("individuals", "SimID"), all.x = T)
        dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("individuals", "SimID"), all.x = T)

        #Pass all F to numeric (just in case)
        dtaFULL.1[,3] = as.numeric(as.character(dtaFULL.1[,3]))
        dtaFULL.1[,4] = as.numeric(as.character(dtaFULL.1[,4]))
        dtaFULL.1[,5] = as.numeric(as.character(dtaFULL.1[,5]))

        #Estimating correlation between FROH SUB and FROH TRUE IBD 100 GEN
        COR100GEN = cor(dtaFULL.1[,3], dtaFULL.1[,6])
        
        #Building model lm(FROH RAD ~ FROH TRUE IBD 100GEN)to extract slope + intercept
        mod100GEN = lm(dtaFULL.1[,3] ~ dtaFULL.1[,6])
        slope100GEN = mod100GEN[[1]][2]
        intercept100GEN = mod100GEN[[1]][1]

        #Estimating correlation between FROH SUB and FROH TRUE IBD 1000 GEN
        COR1000GEN = cor(dtaFULL.1[,3], dtaFULL.1[,7])

        #Building model lm(FROH RAD ~ FROH TRUE IBD 1000GEN)to extract slope + intercept
        mod1000GEN = lm(dtaFULL.1[,3] ~ dtaFULL.1[,7])
        slope1000GEN = mod1000GEN[[1]][2]
        intercept1000GEN = mod1000GEN[[1]][1]

        ### SNP DENSITY

        #read the data
        snpdens = read.table("./SUPP_ANALYSES/SNPdensity_PLINK/SNPs_DENSITY/DENSITY_${SimID}_${RAD_WINDOW_NB}_${REP}.snpden", header = T)
        #Estimate mean NB of SNPs per MB (in KB in the file)
        meanSNPsnbMb = mean(snpdens[,4])*1000

        ### DIFFERENCE

        #Estimating DIFFERENCE between FROH sub and TRUE FROH
        Diff100GEN = dtaFULL.1[,3] - dtaFULL.1[,6]

        #Estimating DIFFERENCE between FROH sub and TRUE FROH
        Diff1000GEN = dtaFULL.1[,3] - dtaFULL.1[,7]

        ### OUTPUTs

        #create output
        output.1 = as.data.frame(matrix(nrow = 1, ncol = 10, c("${SimID}", "${RAD_WINDOW_NB}", "${REP}", COR100GEN, slope100GEN, intercept100GEN, COR1000GEN, slope1000GEN, intercept1000GEN, meanSNPsnbMb)))

        #Writing the dataframe
        write.table(output.1, "./SUPP_ANALYSES/SNPdensity_PLINK/FROH/DENSITY_AND_COR_${SimID}_${RAD_WINDOW_NB}_${REP}_summary_stats.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

        #create output
        output.2 = as.data.frame(matrix(nrow = length(Diff100GEN), ncol = 6, cbind(rep("${SimID}", length(Diff100GEN)), rep("${RAD_WINDOW_NB}", length(Diff100GEN)), rep("${REP}", length(Diff100GEN)), Diff100GEN, Diff1000GEN, rep(meanSNPsnbMb, length(Diff100GEN)))))

        #Writing the dataframe
        write.table(output.2, "./SUPP_ANALYSES/SNPdensity_PLINK/FROH/DIFFERENCE_AND_COR_${SimID}_${RAD_WINDOW_NB}_${REP}_summary_stats.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

EOF


    #Rm PLINK files
    #rm ./SUPP_ANALYSES/SNPdensity_PLINK/PLINK_ROH_${SimID}_${RAD_WINDOW_NB}_${REP}.*

}

export -f SNPDENSITY_RAD_Subsampling_ROH_calling_PLINK

#Launch jobs SMALL POPULATION
parallel -j20 ::: SNPDENSITY_RAD_Subsampling_ROH_calling_PLINK ::: $(find ./Final_VCF -name "FINAL*.vcf.gz") ::: $(seq 100000 2000 600000) ::: {1..100}
#Launch jobs LARGE POPULATION
#parallel -j20 ::: SNPDENSITY_RAD_Subsampling_ROH_calling_PLINK ::: $(find ./Final_VCF -name "FINAL*.vcf.gz") ::: $(seq 10000 200 60000) ::: {1..100}

## Then just merge the replicates into one dataframe
#Create header
echo -e "SimID\tWINDOWS.NB\tREP\tCOR_FROHSUB_TRUEIBD100\tslope_FROHSUB_TRUEIBD100\tintercept_FROHSUB_TRUEIBD100\tCOR_FROHSUB_TRUEIBD1000\tslope_FROHSUB_TRUEIBD1000\tintercept_FROHSUB_TRUEIBD1000\tNBSNPsperMb" > ./Analyses/COR_SLOPE_INTERCEPT_from_SNPDensity_RADseq_PLINK.txt
#Add all files
for file in $(find ./SUPP_ANALYSES/SNPdensity_PLINK/FROH/ -name "DENSITY_AND_COR_*_summary_stats.txt"); do cat ${file} >> ./Analyses/COR_SLOPE_INTERCEPT_from_SNPDensity_RADseq_PLINK.txt; done
#Rm per REPPLICATE file
#rm ./SUPP_ANALYSES/SNPdensity_PLINK/FROH/DENSITY_AND_COR_*_summary_stats.txt

## Then just merge the replicates into one dataframe
#Create header
echo -e "SimID\tWINDOWS.NB\tREP\tDIFF_FROHSUB_TRUEIBD100\tDIFF_FROHSUB_TRUEIBD1000\tNBSNPsperMb" > ./Analyses/DIFF_FROH_from_SNPDensity_RADseq_PLINK.txt
#Add all files
for file in $(find ./SUPP_ANALYSES/SNPdensity_PLINK/FROH/ -name "DIFFERENCE_AND_COR_*_summary_stats.txt"); do cat ${file} >> ./Analyses/DIFF_FROH_from_SNPDensity_RADseq_PLINK.txt; done
#Rm per REPPLICATE file
#rm ./SUPP_ANALYSES/SNPdensity_PLINK/FROH/DIFFERENCE_AND_COR_*_summary_stats.txt
