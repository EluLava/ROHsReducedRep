#!/bin/bash

#Define function
SNPDENSITY_RAD_Subsampling_ROH_calling_RZooRoH() {

    #First argument passed to the function is the WGS VCF
    VCF=$1
    #Extract SimuID from file name
    SimID=$(basename -s ".vcf.gz" ${VCF} | cut -d '_' -f2)

    #Second argument: # of RAD windows we'd like to subsample (≠ percentage of genome subsampled !)
    RAD_WINDOW_NB=$2
    
    #Third argument is replicate
    REP=$3

    #If RZooRoH directory does not exist, create it to store the results of ROHs calling with RZooRoH
    mkdir -p ./SUPP_ANALYSES/SNPdensity_RZooRoH/bedfiles
    mkdir -p ./SUPP_ANALYSES/SNPdensity_RZooRoH/SNPs_DENSITY
    mkdir -p ./SUPP_ANALYSES/SNPdensity_RZooRoH/FROH
    
    #Randomly create BED file with RAD-like fragments (500bp)
    bedtools random -l 500 -n ${RAD_WINDOW_NB} -g ./genome.bed | sort -k1V,1 -k2n,2 > ./SUPP_ANALYSES/SNPdensity_RZooRoH/bedfiles/sorted_${SimID}_${RAD_WINDOW_NB}_${REP}.bed
    
    #Subset SNPs
    vcftools --gzvcf ${VCF} --bed ./SUPP_ANALYSES/SNPdensity_RZooRoH/bedfiles/sorted_${SimID}_${RAD_WINDOW_NB}_${REP}.bed --out ./SUPP_ANALYSES/SNPdensity_RZooRoH/TMP_${SimID}_${RAD_WINDOW_NB}_${REP} --recode
    
    #Estimate SNP_DENSITY
    vcftools --vcf ./SUPP_ANALYSES/SNPdensity_RZooRoH/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}.recode.vcf --SNPdensity 100000 --out ./SUPP_ANALYSES/SNPdensity_RZooRoH/SNPs_DENSITY/DENSITY_${SimID}_${RAD_WINDOW_NB}_${REP}
        
    #convert the file in Oxford format
    plink --vcf ./SUPP_ANALYSES/SNPdensity_RZooRoH/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}.recode.vcf --recode oxford --autosome --chr-set 30 --out ./SUPP_ANALYSES/SNPdensity_RZooRoH/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}_OX

    #Open R for FROH estimation and comparison with WGS estimations
    R --vanilla <<EOF

    library(RZooRoH)
    library(doParallel)
    library(foreach)
    cl <- 10
    registerDoParallel(cl)

    # read the OX file
    data_Rohs = zoodata("./SUPP_ANALYSES/SNPdensity_RZooRoH/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}_OX.gen", zformat = "gp")

    #Create 3 HBD  and 1 non-HBD class.es model
    Mod3R <- zoomodel(K=5, krates=c(10,100,1000,10000,10000), err = 5e-5)
    
    #Run the model
    loc_mod3r <- zoorun(Mod3R, data_Rohs, localhbd = TRUE, nT = cl)
    
    #Extract table with all important information
    loc_table_3R <- loc_mod3r@hbdseg

    #Write this file for TFPN Rates later (not needed for this exra subsampling)
    #write.table(loc_table_3R, "./SUPP_ANALYSES/SNPdensity_RZooRoH/ROHs_SNPdens_${SimID}_${RAD_WINDOW_NB}_${REP}.hom", quote = F, col.names = T, row.names = F)
    

    ### FROH

    #create df
    F3R = 1-loc_mod3r@realized[,5]
    individuals = c(1000, seq(1,(loc_mod3r@nind - 1)))
    df_3R = data.frame(individuals, F3R)
    #Add SimID (fer merging with WGS)
    SUB_dta = cbind(df_3R, SimID=rep("${SimID}", nrow(df_3R)))
    
    ### MERGING with TRUE IBD

    #read the dta TRUE IBD
    dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
    #read the dta TRUE IBD
    dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

    #merge by individuals
    dtaFULL = merge(SUB_dta, dtaTRUE100, by = c("individuals", "SimID"), all.x = T)
    dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("individuals", "SimID"), all.x = T)

    #Pass all F to numeric (just in case)
    dtaFULL.1[,3] = as.numeric(as.character(dtaFULL.1[,3]))
    dtaFULL.1[,4] = as.numeric(as.character(dtaFULL.1[,4]))
    dtaFULL.1[,5] = as.numeric(as.character(dtaFULL.1[,5]))


    #Estimating correlation between FROH SUB and FROH TRUE IBD 100 GEN
    COR100GEN = cor(dtaFULL.1[,3], dtaFULL.1[,4])
    
    #Building model lm(FROH RAD ~ FROH TRUE IBD 100GEN)to extract slope + intercept
    mod100GEN = lm(dtaFULL.1[,3] ~ dtaFULL.1[,4])
    slope100GEN = mod100GEN[[1]][2]
    intercept100GEN = mod100GEN[[1]][1]

    #Estimating correlation between FROH SUB and FROH TRUE IBD 1000 GEN
    COR1000GEN = cor(dtaFULL.1[,3], dtaFULL.1[,5])
    
    #Building model lm(FROH RAD ~ FROH TRUE IBD 1000GEN)to extract slope + intercept
    mod1000GEN = lm(dtaFULL.1[,3] ~ dtaFULL.1[,5])
    slope1000GEN = mod1000GEN[[1]][2]
    intercept1000GEN = mod1000GEN[[1]][1]


    ### SNP DENSITY

    #read the data
    snpdens = read.table("./SUPP_ANALYSES/SNPdensity_RZooRoH/SNPs_DENSITY/DENSITY_${SimID}_${RAD_WINDOW_NB}_${REP}.snpden", header = T)
    #Estimate mean NB of SNPs per MB (in KB in the file)
    meanSNPsnbMb = mean(snpdens[,4])*1000

    ### DIFFERENCE

    #Estimating DIFFERENCE between FROH sub and TRUE FROH
    Diff100GEN = dtaFULL.1[,3] - dtaFULL.1[,4]

    #Estimating DIFFERENCE between FROH sub and TRUE FROH
    Diff1000GEN = dtaFULL.1[,3] - dtaFULL.1[,5]

    ### OUTPUTs

    #create output
    output.1 = as.data.frame(matrix(nrow = 1, ncol = 10, c("${SimID}", "${RAD_WINDOW_NB}", "${REP}", COR100GEN, slope100GEN, intercept100GEN, COR1000GEN, slope1000GEN, intercept1000GEN, meanSNPsnbMb)))

    #Writing the dataframe
    write.table(output.1, "./SUPP_ANALYSES/SNPdensity_RZooRoH/FROH/DENSITY_AND_COR_${SimID}_${RAD_WINDOW_NB}_${REP}_summary_stats.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

    #create output
    output.2 = as.data.frame(matrix(nrow = length(Diff100GEN), ncol = 6, cbind(rep("${SimID}", length(Diff100GEN)), rep("${RAD_WINDOW_NB}", length(Diff100GEN)), rep("${REP}", length(Diff100GEN)), Diff100GEN, Diff1000GEN, rep(meanSNPsnbMb, length(Diff100GEN)))))

    #Writing the dataframe
    write.table(output.2, "./SUPP_ANALYSES/SNPdensity_RZooRoH/FROH/DIFFERENCE_AND_COR_${SimID}_${RAD_WINDOW_NB}_${REP}_summary_stats.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

EOF

    #Rm OX files
    rm ./SUPP_ANALYSES/SNPdensity_RZooRoH/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}_OX.*

}

export -f SNPDENSITY_RAD_Subsampling_ROH_calling_RZooRoH

#Launch jobs SMALL POPULATION
parallel -j20 ::: SNPDENSITY_RAD_Subsampling_ROH_calling_RZooRoH ::: $(find ./Final_VCF -name "FINAL*.vcf.gz") ::: $(seq 400 300 31000) ::: {1..100}
#Launch jobs LARGE POPULATION
#parallel -j20 ::: SNPDENSITY_RAD_Subsampling_ROH_calling_RZooRoH ::: $(find ./Final_VCF -name "FINAL*.vcf.gz") ::: $(seq 100 30 3100) ::: {1..100}

## Then just merge the replicates into one dataframe
#Create header
echo -e "SimID\tWINDOWS.NB\tREP\tCOR_FROHSUB_TRUEIBD100\tslope_FROHSUB_TRUEIBD100\tintercept_FROHSUB_TRUEIBD100\tCOR_FROHSUB_TRUEIBD1000\tslope_FROHSUB_TRUEIBD1000\tintercept_FROHSUB_TRUEIBD1000\tNBSNPsperMb" > ./Analyses/COR_SLOPE_INTERCEPT_from_SNPDensity_RADseq_RZooRoH.txt
#Add all files
for file in $(find ./SUPP_ANALYSES/SNPdensity_RZooRoH/FROH/ -name "DENSITY_AND_COR_*_summary_stats.txt"); do cat ${file} >> ./Analyses/COR_SLOPE_INTERCEPT_from_SNPDensity_RADseq_RZooRoH.txt; done
#Rm per REPPLICATE file
rm ./SUPP_ANALYSES/SNPdensity_RZooRoH/FROH/DENSITY_AND_COR_*_summary_stats.txt

## Then just merge the replicates into one dataframe
#Create header
echo -e "SimID\tWINDOWS.NB\tREP\tDIFF_FROHSUB_TRUEIBD100\tDIFF_FROHSUB_TRUEIBD1000\tNBSNPsperMb" > ./Analyses/DIFF_FROH_from_SNPDensity_RADseq_RZooRoH.txt
#Add all files
for file in $(find ./SUPP_ANALYSES/SNPdensity_RZooRoH/FROH/ -name "DIFFERENCE_AND_COR_*_summary_stats.txt"); do cat ${file} >> ./Analyses/DIFF_FROH_from_SNPDensity_RADseq_RZooRoH.txt; done
#Rm per REPPLICATE file
rm ./SUPP_ANALYSES/SNPdensity_RZooRoH/FROH/DIFFERENCE_AND_COR_*_summary_stats.txt


