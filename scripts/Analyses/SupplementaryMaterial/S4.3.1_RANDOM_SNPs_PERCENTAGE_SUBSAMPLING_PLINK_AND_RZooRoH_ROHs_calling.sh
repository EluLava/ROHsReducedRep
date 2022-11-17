#Define function
SUBPERC_Subsampling_ROH_calling_RZooRoH() {

    #First argument passed to the function is the WGS VCF
    VCF=$1
    #Extract SimuID from file name
    SimID=$(basename -s ".vcf.gz" ${VCF} | cut -d '_' -f2)

    #Second argument: percentage of SNPs (≠ genome) subsampled
    PERCENTAGE=$2
    
    #Third argument is replicate
    REP=$3

    #If RZOOROH directory does not exist, create it to store the results of ROHs calling with PLINK
    mkdir -p ./SUPP_ANALYSES/RANDOM_PERC/RZooRoH
    mkdir -p ./SUPP_ANALYSES/RANDOM_PERC/PLINK
    
    #Extract SNPs list IF DOES NOT EXISTs yet
    if [ ! -f "./SUPP_ANALYSES/RANDOM_PERC/Full_SNPsList_${SimID}.tsv" ]; then
        bcftools query -f'%CHROM\t%POS\n' ${VCF} > ./SUPP_ANALYSES/RANDOM_PERC/Full_SNPsList_${SimID}.tsv
    fi

    # Subsample random SNPs
    shuf -n $(($PERCENTAGE*`wc -l ./SUPP_ANALYSES/RANDOM_PERC/Full_SNPsList_${SimID}.tsv | cut -d' ' -f1`/100)) ./SUPP_ANALYSES/RANDOM_PERC/Full_SNPsList_${SimID}.tsv | sort -k1,1n -k2,2n > ./SUPP_ANALYSES/RANDOM_PERC/SNPsList_${SimID}_${PERCENTAGE}_${REP}.tsv
    
    # Then subsample the SNPs
    bcftools view -R ./SUPP_ANALYSES/RANDOM_PERC/SNPsList_${SimID}_${PERCENTAGE}_${REP}.tsv -O z -o ./SUPP_ANALYSES/RANDOM_PERC/TMPVCF_${SimID}_${PERCENTAGE}_${REP}.vcf.gz ${VCF}

    # Estimate SNP density for PLINK parameters

    #Estimate SNP_DENSITY
    vcftools --gzvcf ./SUPP_ANALYSES/RANDOM_PERC/TMPVCF_${SimID}_${PERCENTAGE}_${REP}.vcf.gz --SNPdensity 100000 --out ./SUPP_ANALYSES/RANDOM_PERC/DENSITY_${SimID}_${PERCENTAGE}_${REP}

    #Average SNP densities per windows --> Average SNP density per KB
    SNPdens=$(awk '{ sum += $4 } END { if (NR > 0) print sum / NR }' ./SUPP_ANALYSES/RANDOM_PERC/DENSITY_${SimID}_${PERCENTAGE}_${REP}.snpden)

    #rm SNP density file
    rm ./SUPP_ANALYSES/RANDOM_PERC/DENSITY_${SimID}_${PERCENTAGE}_${REP}.*
    #rm SNP list file
    rm ./SUPP_ANALYSES/RANDOM_PERC/SNPsList_${SimID}_${PERCENTAGE}_${REP}.tsv

    ## LEt's start by PLINK ROHs calling

    #Number of heterozygous sites allowed in one ROH will depend on PERCENTAGE ! --> minimum 2, max 9 --> from 0% to 20% = 2 SNPs, then each 10% increase, one more het SNP allowed in a ROH, max 8
    #If PERC smaller or equal to 20, het = 2
    if [[ `printf "%.0f" $(echo "${PERCENTAGE}" | bc)` -le 20 ]]; then

        hetallowed=2

    #If higher or equal to 80
    elif [[ `printf "%.0f" $(echo "${PERCENTAGE}" | bc)` -ge 80 ]]; then

        hetallowed=8

    #Else simply percentage / 10
    else

        hetallowed=$((${PERCENTAGE}/10))

    fi
    
    #Convert from vcf to bed
    plink --vcf ./SUPP_ANALYSES/RANDOM_PERC/TMPVCF_${SimID}_${PERCENTAGE}_${REP}.vcf.gz --make-bed --autosome --chr-set 30 --out ./SUPP_ANALYSES/RANDOM_PERC/BEDTMP_${SimID}_${PERCENTAGE}_${REP}

    #If SNP dens (integer format) lower or equal to 30 , simply use 30
    if [[ `printf "%.0f" $(echo "${SNPdens}" | bc)` -le 30 ]]; then

        #Call ROHs
        plink --bfile ./SUPP_ANALYSES/RANDOM_PERC/BEDTMP_${SimID}_${PERCENTAGE}_${REP} --homozyg --homozyg-window-snp 30 --homozyg-density 50 \
        --homozyg-snp 30 --homozyg-gap 1000 --homozyg-window-het 1 --homozyg-het `echo ${hetallowed}` --homozyg-kb 100 --chr-set 30 --out ./SUPP_ANALYSES/RANDOM_PERC/PLINK/PLINK_ROH_${SimID}_${PERCENTAGE}_${REP}

    #else, if greater or equal to 100
    elif [[ `printf "%.0f" $(echo "${SNPdens}" | bc)` -ge 100 ]]; then

        #Call ROHs
        plink --bfile ./SUPP_ANALYSES/RANDOM_PERC/BEDTMP_${SimID}_${PERCENTAGE}_${REP} --homozyg --homozyg-window-snp 100 --homozyg-density 100 \
        --homozyg-snp 100 --homozyg-gap 1000 --homozyg-window-het 1 --homozyg-het `echo ${hetallowed}` --homozyg-kb 100 --chr-set 30 --out ./SUPP_ANALYSES/RANDOM_PERC/PLINK/PLINK_ROH_${SimID}_${PERCENTAGE}_${REP}

    #Finally, else (--> if between 30 and 100) use params (--homozyg-window-snp + --homozyg-snp) according to SNP dens
    else

        plink --bfile ./SUPP_ANALYSES/RANDOM_PERC/BEDTMP_${SimID}_${PERCENTAGE}_${REP} --homozyg --homozyg-window-snp `printf "%.0f" $(echo "${SNPdens}*100" | bc)` --homozyg-density 50 \
        --homozyg-snp `printf "%.0f" $(echo "${SNPdens}*100" | bc)` --homozyg-gap 1000 --homozyg-window-het 1 --homozyg-het `echo ${hetallowed}` --homozyg-kb 100 --chr-set 30 --out ./SUPP_ANALYSES/RANDOM_PERC/PLINK/PLINK_ROH_${SimID}_${PERCENTAGE}_${REP}

    fi

    ## Then RZooRoH HBD

    #convert the file in Oxford format
    plink --vcf ./SUPP_ANALYSES/RANDOM_PERC/TMPVCF_${SimID}_${PERCENTAGE}_${REP}.vcf.gz --recode oxford --autosome --chr-set 30 --out ./SUPP_ANALYSES/RANDOM_PERC/TMP_${SimID}_${PERCENTAGE}_${REP}_OX

    #rm VCF
    rm ./SUPP_ANALYSES/RANDOM_PERC/TMPVCF_${SimID}_${PERCENTAGE}_${REP}.vcf.gz
    #rm BED file
    rm ./SUPP_ANALYSES/RANDOM_PERC/BEDTMP_${SimID}_${PERCENTAGE}_${REP}.*

    #Open R to do the modelling
    R --vanilla <<EOF

    library(RZooRoH)
    library(doParallel)
    library(foreach)
    cl <- 10
    registerDoParallel(cl)

    ## RZooRoH

    # read the OX file
    data_Rohs = zoodata("./SUPP_ANALYSES/RANDOM_PERC/TMP_${SimID}_${PERCENTAGE}_${REP}_OX.gen", zformat = "gp")

    #Create 3 HBD  and 1 non-HBD class.es model
    Mod3R <- zoomodel(K=5, krates=c(10,100,1000,10000,10000), err = 5e-5)

    #Run the model
    loc_mod3r <- zoorun(Mod3R, data_Rohs, localhbd = TRUE, nT = cl)

    #Extract table with all important information
    loc_table_3R <- loc_mod3r@hbdseg

    #Change indv names so they match simulations + plink ouputs
    loc_table_3R[,1][loc_table_3R[,1] == 1] = 1000
    loc_table_3R[,1][loc_table_3R[,1] != 1000] = loc_table_3R[,1][loc_table_3R[,1] != 1000] - 1
    
    #Write this file for TFPN Rates later
    write.table(loc_table_3R, "./SUPP_ANALYSES/RANDOM_PERC/RZooRoH/PERCTMP_${SimID}_${PERCENTAGE}_${REP}.hom", quote = F, col.names = T, row.names = F)    

    ## FROH

    #create df
    F3R = 1-loc_mod3r@realized[,5]
    individuals = c(1000, seq(1,(loc_mod3r@nind - 1)))
    df_3R = data.frame(individuals, F3R)

    # cbind into output df
    dt_3R <- cbind(MODEL = rep("4R", nrow(df_3R)), SimID=rep("${SimID}" , nrow(df_3R)), df_3R, SNPs.PERCENTAGE=rep("${PERCENTAGE}", nrow(df_3R)), REP=rep("${REP}" , nrow(df_3R)))
    colnames(dt_3R) <- c("Model", "SimID" , "individuals" , "FROH_sub" , "SNPs.PERCENTAGE" , "REP")

    #Save the data
    write.table(dt_3R, "./SUPP_ANALYSES/RANDOM_PERC/RZooRoH/PERCTMP_${SimID}_${PERCENTAGE}_${REP}_summary_stats.txt", col.names = F , row.names = F , quote = F)


    ## ROHs DIST

    #Create ROHs DIST table
    dta_dist <- cbind(rep("4R", nrow(loc_table_3R)), rep("${SimID}" , nrow(loc_table_3R)), rep("${PERCENTAGE}" , nrow(loc_table_3R)), rep("${REP}" , nrow(loc_table_3R)), loc_table_3R)
    colnames(dta_dist) <- c("Model", "SimID" , "SNPs.PERCENTAGE" ,"REP" , "individuals" , "chrom" , "start_snp", "end_snp" , "start_pos" , "end_pos" , "number_snp" , "length" , "HBD")

    #Create the output df
    dta_roh = as.data.frame(matrix(nrow=6,ncol=6))
    colnames(dta_roh) = c("SimID","SNPs.PERCENTAGE","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")

    #Fill SimuID, CALL_METH, REP, PERC
    dta_roh[,1] = rep("${SimID}",nrow(dta_roh))
    dta_roh[,2] = rep("${PERCENTAGE}",nrow(dta_roh))
    dta_roh[,3] = rep("${REP}",nrow(dta_roh))
    dta_roh[,4] = c(1:6)

    if(  == 0){

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
            if(dta_dist[i,12] < 2000000){
                output[1,2] = 1
            } else if(dta_dist[i,12] >= 2000000 & dta_dist[i,12] < 4000000) {
                output[1,2] = 2
            } else if(dta_dist[i,12] >= 4000000 & dta_dist[i,12] < 6000000) {
                output[1,2] = 3
            } else if(dta_dist[i,12] >= 6000000 & dta_dist[i,12] < 10000000) {
                output[1,2] = 4
            } else if(dta_dist[i,12] >= 10000000 & dta_dist[i,12] < 16000000) {
                output[1,2] = 5
            } else if(dta_dist[i,12] >= 16000000) {
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
        x = aggregate(dta_dist.2[,12], by = list(INDIV=dta_dist.2[,5], CLASS=dta_dist.2[,14]), FUN = sum)
        y = aggregate(dta_dist.2[,12], by = list(INDIV=dta_dist.2[,5], CLASS=dta_dist.2[,14]), FUN = length)


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

    #write output
    write.table(dta_roh, file = "./SUPP_ANALYSES/RANDOM_PERC/RZooRoH/PERCTMP_${SimID}_${PERCENTAGE}_${REP}_ROHsDist.txt", quote = F, row.names = F, col.names = F)


    ## PLINK

    #Read the indv roh plink file
    dtaroh = read.table("./SUPP_ANALYSES/RANDOM_PERC/PLINK/PLINK_ROH_${SimID}_${PERCENTAGE}_${REP}.hom.indiv", header = T)
    #Read the rohs dist file
    dtaroh_dist = read.table("./SUPP_ANALYSES/RANDOM_PERC/PLINK/PLINK_ROH_${SimID}_${PERCENTAGE}_${REP}.hom", header = T)[,-c(3,11:13)]

    #Set the genome length IN KB
    genome_length = 3000000
    
    #Set individuals names (again dtaroh here because some indvs might have 0 ROHs)
    individuals= dtaroh[,2]

    ## FROH

    #Estimate FROH
    Froh = apply(dtaroh, 1, function(x) as.numeric(as.character(x[5]))/genome_length)

    #First output dataframe
    F_sub = cbind(individuals, Froh, dtaroh[,4], dtaroh[,5])
    colnames(F_sub) = c("individuals", "FROH_sub", "NROH_sub", "SROH_sub")

    #Merging SimID, Froh_sub, Number of windows and SNPs & REP
    SUB_dta = as.data.frame(cbind(rep("${SimID}",nrow(F_sub)), F_sub, rep("${PERCENTAGE}",nrow(F_sub)),rep("${REP}",nrow(F_sub))))
    colnames(SUB_dta) = c("SimID", "individuals", "FROH_sub", "NROH_sub", "SROH_sub", "SNPs.PERCENTAGE", "REP")

    #Writing the dataframe
    write.table(SUB_dta, "./SUPP_ANALYSES/RANDOM_PERC/PLINK/PERCTMP_${SimID}_${PERCENTAGE}_${REP}_summary_stats.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

    ## ROHs Dist

    #Create the output dataframe which will contain lengths per ROH category
    dta_roh = as.data.frame(matrix(nrow = 6, ncol = 6))
    colnames(dta_roh) = c("SimuID","SNPs.PERCENTAGE","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")

    #Fill SimuID, REP, and WINDOWS.NB
    dta_roh[,1] = rep("${SimID}",nrow(dta_roh))
    dta_roh[,2] = rep("${PERCENTAGE}",nrow(dta_roh))
    dta_roh[,3] = rep("${REP}",nrow(dta_roh))
    dta_roh[,4] = c(1:6)

    #If no ROHs called all is 0 !
    if(nrow(dtaroh_dist) == 0) {

        #just put 0 in output d
        dta_roh[,5] = 0
        dta_roh[,6] = 0

    } else {

        #Loop through ROHs and Add the roh class
        foreachoutput = foreach(i=1:nrow(dtaroh_dist), .combine=rbind) %dopar% {
        
            #Create output
            output = as.data.frame(matrix(nrow=1,ncol=2))
            colnames(output) = c("i","ROH_CLASS")
            output[1,1] = i
            
            #Define ROH class
            if(dtaroh_dist[i,8] < 2000){
                output[1,2] = 1
            } else if(dtaroh_dist[i,8] >= 2000 & dtaroh_dist[i,8] < 4000) {
                output[1,2] = 2
            } else if(dtaroh_dist[i,8] >= 4000 & dtaroh_dist[i,8] < 6000) {
                output[1,2] = 3
            } else if(dtaroh_dist[i,8] >= 6000 & dtaroh_dist[i,8] < 10000) {
                output[1,2] = 4
            } else if(dtaroh_dist[i,8] >= 10000 & dtaroh_dist[i,8] < 16000) {
                output[1,2] = 5
            } else if(dtaroh_dist[i,8] >= 16000) {
                output[1,2] = 6
            }

            #return the rownb and class
            return(output)
        }

        #Sort output of foreach per i
        foreachoutput.2 = foreachoutput[order(foreachoutput[,1]),]

        #Merge with dtaroh_dist
        dta_dist.2 = cbind(dtaroh_dist, ROH_CLASS=foreachoutput.2[,2])
        
        #For each ROH category estimate sum of lengths and mean sum of length for each individual
        x = aggregate(dta_dist.2[,8], by = list(INDIV=dta_dist.2[,2], CLASS=dta_dist.2[,10]), FUN = sum)
        y = aggregate(dta_dist.2[,8], by = list(INDIV=dta_dist.2[,2], CLASS=dta_dist.2[,10]), FUN = length)

        #Count the number of ROHs and and sum Length for each catgory
        for(roh_class in 1:6) {
            #if ROHs class was detected
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
    write.table(dta_roh, file = "./SUPP_ANALYSES/RANDOM_PERC/PLINK/PERCTMP_${SimID}_${PERCENTAGE}_${REP}_ROHsDist.txt", quote = F, row.names = F, col.names = F)

EOF

    #rm TMP OX files
    rm ./SUPP_ANALYSES/RANDOM_PERC/TMP_${SimID}_${PERCENTAGE}_${REP}_OX.*

}

export -f SUBPERC_Subsampling_ROH_calling_RZooRoH

#Launch jobs
parallel -j20 ::: SUBPERC_Subsampling_ROH_calling_RZooRoH ::: $(find ./Final_VCF -name "FINAL*.vcf.gz") ::: {2,3,4,5,10,20,30,40,50,60,70,80,90} ::: {1..100}

## Then just merge the SIM replicates into one dataframe

## RZooRoH

## FROH

#Create header
echo -e "Model\nSimID\nindividuals\nFROH_sub\nSNPs.PERCENTAGE\nREP") > ./SUPP_ANALYSES/RANDOM_PERC/RZooRoH/RandomPerc_RZooRoH_FROH.txt
#Add all files
for file in $(find ./SUPP_ANALYSES/RANDOM_PERC/RZooRoH/ -name "PERCTMP_*_summary_stats.txt"); do cat ${file} >> ./SUPP_ANALYSES/RANDOM_PERC/RZooRoH/RandomPerc_RZooRoH_FROH.txt; done
#Rm per SimID file
rm ./SUPP_ANALYSES/RANDOM_PERC/RZooRoH/PERCTMP_*_summary_stats.txt

#Merge with TRUE IBD (R)

#Open R for merging
R --vanilla <<EOF

    #Read RAD-seq + TRUE IBD data
    dtasub = read.table("./SUPP_ANALYSES/RANDOM_PERC/RZooRoH/RandomPerc_RZooRoH_FROH.txt", header = T)
    dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
    dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

    #merge both (without windows nb and REP for WGS)
    dtaFULL = merge(dtasub, dtaTRUE100, by = c("SimID","individuals"))
    dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("SimID","individuals"))
    #Save the data
    write.table(dtaFULL.1, "./Analyses/SM_RANDOMPERC_RZooRoH_FROH.txt", quote = F, col.names = T, row.names = F)

EOF

## ROHs Dist

#Create header
echo -e "SimID\tSNPs.PERCENTAGE\tREP\tROHs_CLASS\tMean_NROH_per_ind\tMean_Tot_Kb_per_ind" > ./Analyses/SM_RANDOMPERC_RZooRoH_ROHsDistribution.txt
#Add all files
for file in $(find ./SUPP_ANALYSES/RANDOM_PERC/RZooRoH/ -name "PERCTMP_*_ROHsDist.txt"); do cat ${file} >> ./Analyses/SM_RANDOMPERC_RZooRoH_ROHsDistribution.txt; done
#Rm per SimID file
rm ./SUPP_ANALYSES/RANDOM_PERC/RZooRoH/PERCTMP_*_ROHsDist.txt


## PLINK

## FROH

#Create header
echo -e "SimID\tindividuals\tFROH_sub\tNROH_sub\tSROH_sub\SNPs.PERCENTAGE\tSNPs_NB\tREP" > ./SUPP_ANALYSES/RANDOM_PERC/PLINK/RandomPerc_PLINK_FROH.txt
#Add all files
for file in $(find ./SUPP_ANALYSES/RANDOM_PERC/PLINK/ -name "PERCTMP_*_summary_stats.txt"); do cat ${file} >> ./SUPP_ANALYSES/RANDOM_PERC/PLINK/RandomPerc_PLINK_FROH.txt; done
#Rm per SimID file
rm ./SUPP_ANALYSES/RANDOM_PERC/PLINK/PERCTMP_*_summary_stats.txt

#Merge with TRUE IBD (R)

#Open R for merging
R --vanilla <<EOF

    #Read RAD-seq + TRUE IBD data
    dtasub = read.table("./SUPP_ANALYSES/RANDOM_PERC/PLINK/RandomPerc_PLINK_FROH.txt", header = T)
    dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
    dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

    #merge both (without windows nb and REP for WGS)
    dtaFULL = merge(dtasub, dtaTRUE100, by = c("SimID","individuals"))
    dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("SimID","individuals"))
    #Save the data
    write.table(dtaFULL.1, "./Analyses/SM_RANDOMPERC_PLINK_FROH.txt", quote = F, col.names = T, row.names = F)

EOF

## ROHs Dist

#Create header
echo -e "SimID\tSNPs.PERCENTAGE\tREP\tROHs_CLASS\tMean_NROH_per_ind\tMean_Tot_Kb_per_ind" > ./Analyses/SM_RANDOMPERC_PLINK_ROHsDistribution.txt
#Add all files
for file in $(find ./SUPP_ANALYSES/RANDOM_PERC/PLINK/ -name "PERCTMP_*_ROHsDist.txt"); do cat ${file} >> ./Analyses/SM_RANDOMPERC_PLINK_ROHsDistribution.txt; done
#Rm per SimID file
rm ./SUPP_ANALYSES/RANDOM_PERC/PLINK/PERCTMP_*_ROHsDist.txt



