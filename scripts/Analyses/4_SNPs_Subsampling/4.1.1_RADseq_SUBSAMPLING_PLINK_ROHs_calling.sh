#!/bin/bash

#Define function
RAD_Subsampling_ROH_calling_PLINK() {

    #First argument passed to the function is the WGS VCF
    VCF=$1
    #Extract SimuID from file name
    SimID=$(basename -s ".vcf.gz" ${VCF} | cut -d '_' -f2)

    #Second argument: # of RAD windows we'd like to subsample (≠ percentage of genome subsampled !)
    RAD_WINDOW_NB=$2
    
    #Third argument is replicate
    REP=$3

    #If PLINK directory does not exist, create it to store the results of ROHs calling with PLINK
    mkdir -p ./RADseq_PLINK/bedfiles
    mkdir -p ./RADseq_PLINK/SNPs_DENSITY
    
    #Randomly create BED file with RAD-like fragments (500bp)
    bedtools random -l 500 -n ${RAD_WINDOW_NB} -g ./genome.bed | sort -k1V,1 -k2n,2 > ./RADseq_PLINK/bedfiles/sorted_${SimID}_${RAD_WINDOW_NB}_${REP}.bed
    
    #Subset SNPs
    vcftools --gzvcf ${VCF} --bed ./RADseq_PLINK/bedfiles/sorted_${SimID}_${RAD_WINDOW_NB}_${REP}.bed --out ./RADseq_PLINK/TMP_${SimID}_${RAD_WINDOW_NB}_${REP} --recode
    
    #Record number of subsampled SNPs
    SNPs_NB=$(grep -v '#' ./RADseq_PLINK/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}.recode.vcf | wc -l)
    
    #Save SNP_NUMBERs in bedfile name
    mv ./RADseq_PLINK/bedfiles/sorted_${SimID}_${RAD_WINDOW_NB}_${REP}.bed ./RADseq_PLINK/bedfiles/sorted_${SimID}_${RAD_WINDOW_NB}_${SNPs_NB}_${REP}.bed
    
    #Estimate SNP_DENSITY
    vcftools --vcf ./RADseq_PLINK/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}.recode.vcf --SNPdensity 100000 --out ./RADseq_PLINK/SNPs_DENSITY/DENSITY_${SimID}_${RAD_WINDOW_NB}_${REP}

    #Average SNP densities per windows --> Average SNP density per KB
    SNPdens=$(awk '{ sum += $4 } END { if (NR > 0) print sum / NR }' ./RADseq_PLINK/SNPs_DENSITY/DENSITY_${SimID}_${RAD_WINDOW_NB}_${REP}.snpden)

    #From VCF to bed
    plink --vcf ./RADseq_PLINK/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}.recode.vcf --make-bed --chr-set 30 --out ./RADseq_PLINK/BTMP_${SimID}_${RAD_WINDOW_NB}_${REP}

    #If SNP dens (integer format) lower or equal to 30 , simply use 30
    if [[ `printf "%.0f" $(echo "${SNPdens}" | bc)` -le 30 ]]; then

        #Call ROHs
        plink --bfile ./RADseq_PLINK/BTMP_${SimID}_${RAD_WINDOW_NB}_${REP} --homozyg --homozyg-window-snp 30 --homozyg-density 50 \
        --homozyg-snp 30 --homozyg-gap 1000 --homozyg-window-het 1 --homozyg-het 2 --homozyg-kb 100 --chr-set 30 --out ./RADseq_PLINK/PLINK_ROH_${SimID}_${RAD_WINDOW_NB}_${REP}
  
    #else, if greater or equal to 100
    elif [[ `printf "%.0f" $(echo "${SNPdens}" | bc)` -ge 100 ]]; then

        #Call ROHs
        plink --bfile ./RADseq_PLINK/BTMP_${SimID}_${RAD_WINDOW_NB}_${REP} --homozyg --homozyg-window-snp 100 --homozyg-density 50 \
        --homozyg-snp 100 --homozyg-gap 1000 --homozyg-window-het 1 --homozyg-het 2 --homozyg-kb 100 --chr-set 30 --out ./RADseq_PLINK/PLINK_ROH_${SimID}_${RAD_WINDOW_NB}_${REP}

    #Finally, else (--> if between 30 and 100) use params (--homozyg-window-snp + --homozyg-snp) according to SNP dens
    else

        plink --bfile ./RADseq_PLINK/BTMP_${SimID}_${RAD_WINDOW_NB}_${REP} --homozyg --homozyg-window-snp `printf "%.0f" $(echo "${SNPdens}*100" | bc)` --homozyg-density 50 \
        --homozyg-snp `printf "%.0f" $(echo "${SNPdens}*100" | bc)` --homozyg-gap 1000 --homozyg-window-het 1 --homozyg-het 2 --homozyg-kb 100 --chr-set 30 --out ./RADseq_PLINK/PLINK_ROH_${SimID}_${RAD_WINDOW_NB}_${REP}


    fi

    #Calculate Het
    plink --bfile ./RADseq_PLINK/BTMP_${SimID}_${RAD_WINDOW_NB}_${REP} --het --chr-set 30 --out ./RADseq_PLINK/HET_${SimID}_${RAD_WINDOW_NB}_${REP}

    #Remove TMPVCF file
    rm ./RADseq_PLINK/TMP_${SimID}_${RAD_WINDOW_NB}_${REP}.*
    #Remove BED, BIM; BAM ETc.
    rm ./RADseq_PLINK/BTMP_${SimID}_${RAD_WINDOW_NB}_${REP}.*
    
    #Open R for F and ROHs distributions estimations
    R --vanilla <<EOF

        library(doParallel)
        library(foreach)
        cl <- 3
        registerDoParallel(cl)

        #Read the indv roh plink file
        dtaroh = read.table("./RADseq_PLINK/PLINK_ROH_${SimID}_${RAD_WINDOW_NB}_${REP}.hom.indiv", header = T)
        #Read the het plink file
        dtahet = read.table("./RADseq_PLINK/HET_${SimID}_${RAD_WINDOW_NB}_${REP}.het", header = T)
        #Read the rohs dist file
        dta_dist = read.table("./RADseq_PLINK/PLINK_ROH_${SimID}_${RAD_WINDOW_NB}_${REP}.hom", header = T)[,-c(3,11:13)]

        #Set the genome length IN KB
        genome_length = 3000000
        #Set individuals nb (again het here becasue some indvs might have 0 ROHs with RAD)
        nb_indiv=nrow(dtahet)

        #Set the individuals vector
        individuals = dtaroh[,2]

        ## FROH

        #Estimate FROH
        Froh = apply(dtaroh, 1, function(x) as.numeric(as.character(x[5]))/genome_length)

        #First output dataframe
        F_sub = cbind(individuals, Froh, dtaroh[,4], dtaroh[,5])
        colnames(F_sub) = c("individuals", "FROH_sub", "NROH_sub", "SROH_sub")

        #Merging SimID, Froh_sub, Number of windows and SNPs & REP
        SUB_dta = as.data.frame(cbind(rep("${SimID}",nrow(F_sub)), F_sub, rep("${RAD_WINDOW_NB}",nrow(F_sub)), rep("${SNPs_NB}",nrow(F_sub)),rep("${REP}",nrow(F_sub))))
        colnames(SUB_dta) = c("SimID", "individuals", "FROH_sub", "NROH_sub", "SROH_sub", "WINDOWS.NB", "SNPs_NB", "REP")

        #Writing the dataframe
        write.table(SUB_dta, "./RADseq_PLINK/RADTMP_${SimID}_${RAD_WINDOW_NB}_${REP}_summary_stats.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

        ## ROHs Dist

        #Create the output dataframe which will contain lengths per ROH category
        dta_roh = as.data.frame(matrix(nrow = 6, ncol = 6))
        colnames(dta_roh) = c("SimuID","WINDOWS.NB","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")

        #Fill SimuID, REP, and WINDOWS.NB
        dta_roh[,1] = rep("${SimID}",nrow(dta_roh))
        dta_roh[,2] = rep("${RAD_WINDOW_NB}",nrow(dta_roh))
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
        write.table(dta_roh, file = "./RADseq_PLINK/RADTMP_${SimID}_${RAD_WINDOW_NB}_${REP}_ROHsDist.txt", quote = F, row.names = F, col.names = F)
EOF

}

export -f RAD_Subsampling_ROH_calling_PLINK

#Launch jobs
#SmallPop
parallel -j20 ::: RAD_Subsampling_ROH_calling_PLINK ::: $(find ./Final_VCF -name "FINAL*.vcf.gz") ::: {180000,195000,210000,225000,240000,255000,270000,295000,450000,600000} ::: {1..100}
#LargePop
#parallel -j20 ::: RAD_Subsampling_ROH_calling_PLINK ::: $(find ./Final_VCF -name "FINAL*.vcf.gz") ::: {10000,15000,20000,25000,30000,35000,40000,45000,50000,55000,60000} ::: {1..100}

## Then just merge the SIM replicates into one dataframe

## FHOM
#Create header
echo -e "SimID\tWINDOWS.NB\tFID\tIID\tO(HOM)_sub\tE(HOM)_sub\tN(NM)_sub\tFHOM_sub" > ./RADseq_PLINK/RADSEQ_HET.txt
#Add all files (and SimID as first column)
for file in $(find ./RADseq_PLINK -name "HET_*.het"); do SimID=$(basename -s ".het" ${file} | cut -d '_' -f2); NBRADFRAG=$(basename -s ".het" ${file} | cut -d '_' -f3); awk -v simid=${SimID} -v nbradfrag=${NBRADFRAG} 'NR > 1 {print simid"\t"nbradfrag"\t"$0}' ${file} >> ./RADseq_PLINK/RADSEQ_HET.txt; done
#Rm per SIMID file
rm ./RADseq_PLINK/HET_*

#Merge with TRUE IBD (R)

#Create directory if needed
mkdir -p ./Analyses

#Open R for merging
R --vanilla <<EOF

    #Read RAD-seq + TRUE IBD data
    dtasub = read.table("./RADseq_PLINK/RADSEQ_HET.txt", header = T)
    colnames(dtasub)[4] = "individuals"
    dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
    dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

    #merge both (without windows nb and REP for WGS)
    dtaFULL = merge(dtasub, dtaTRUE100, by = c("SimID","individuals"))
    dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("SimID","individuals"))
    
    #Save the data
    write.table(dtaFULL.1, "./Analyses/RADSEQ_PLINK_FHOM.txt", quote = F, col.names = T, row.names = F)

EOF


## FROH

#Create header
echo -e "SimID\tindividuals\tFROH_sub\tNROH_sub\tSROH_sub\tWINDOWS.NB\tSNPs_NB\tREP" > ./RADseq_PLINK/RADSEQ_PLINK_FROH.txt
#Add all files
for file in $(find ./RADseq_PLINK -name "RADTMP_*_summary_stats.txt"); do cat ${file} >> ./RADseq_PLINK/RADSEQ_PLINK_FROH.txt; done
#Rm per SimID file
rm ./RADseq_PLINK/RADTMP_*_summary_stats.txt

#Merge with TRUE IBD (R)

#Open R for merging
R --vanilla <<EOF

    #Read RAD-seq + TRUE IBD data
    dtasub = read.table("./RADseq_PLINK/RADSEQ_PLINK_FROH.txt", header = T)
    dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
    dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

    #merge both (without windows nb and REP for WGS)
    dtaFULL = merge(dtasub, dtaTRUE100, by = c("SimID","individuals"))
    dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("SimID","individuals"))

    #Save the data
    write.table(dtaFULL.1, "./Analyses/RADSEQ_PLINK_FROH.txt", quote = F, col.names = T, row.names = F)

EOF

## ROHs Dist

#Create header
echo -e "SimID\tWINDOWS.NB\tREP\tROHs_CLASS\tMean_NROH_per_ind\tMean_Tot_Kb_per_ind" > ./Analyses/RADSEQ_PLINK_ROHsDistribution.txt
#Add all files
for file in $(find ./RADseq_PLINK/ -name "RADTMP_*_ROHsDist.txt"); do cat ${file} >> ./Analyses/RADSEQ_PLINK_ROHsDistribution.txt; done
#Rm per SimID file
rm ./RADseq_PLINK/RADTMP_*_ROHsDist.txt

