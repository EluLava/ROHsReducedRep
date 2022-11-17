#!/bin/bash

#Define function
RAD_Subsampling_ROH_calling_RZooRoH() {

	#First argument passed to the function is VCF WGS file
	VCF=$1
	#Extract SimID from VCF filename
	SimID=$(basename -s ".vcf.gz" ${VCF} | cut -d'_' -f2)
	
	#Second argument is # of 500bp fragments for RD (≠ percentage genome subsampled)
	NB_RAD_FRAG=$2
	#Thirs argument replicate
	REP=$3

	#If RZooRoH directory does not exist, create it to store the results of ROHs calling with RZooRoH
	mkdir -p ./RADseq_RZooRoH/bedfiles
	mkdir -p ./RADseq_RZooRoH/SNPs_DENSITY

	#Create bedfile
	bedtools random -l 500 -n ${NB_RAD_FRAG}  -g ./genome.bed | sort -k1V,1 -k2n,2 > ./RADseq_RZooRoH/bedfiles/sorted_${SimID}_${NB_RAD_FRAG}_${REP}.bed

	#Subsample SNPs
	vcftools --gzvcf ${VCF} --bed ./RADseq_RZooRoH/bedfiles/sorted_${SimID}_${NB_RAD_FRAG}_${REP}.bed --out ./RADseq_RZooRoH/TMP_${SimID}_${NB_RAD_FRAG}_${REP} --recode

	#Record number of subsampled SNPs
	SNPs_NB=$(grep -v '#' ./RADseq_RZooRoH/TMP_${SimID}_${NB_RAD_FRAG}_${REP}.recode.vcf | wc -l)

	#Get SNP_NUMBERs in bedfile name
	mv ./RADseq_RZooRoH/bedfiles/sorted_${SimID}_${NB_RAD_FRAG}_${REP}.bed ./RADseq_RZooRoH/bedfiles/sorted_${SimID}_${NB_RAD_FRAG}_${SNPs_NB}_${REP}.bed
	
	#Estimate SNP_DENSITY
	vcftools --vcf ./RADseq_RZooRoH/TMP_${SimID}_${NB_RAD_FRAG}_${REP}.recode.vcf --SNPdensity 100000 --out ./RADseq_RZooRoH/SNPs_DENSITY/DENSITY_${SimID}_${NB_RAD_FRAG}_${REP}
	
	#SNPs distributions among windows
	bedtools coverage -a ./RADseq_RZooRoH/bedfiles/sorted_${SimID}_${NB_RAD_FRAG}_${SNPs_NB}_${REP}.bed -b ./RADseq_RZooRoH/TMP_${SimID}_${NB_RAD_FRAG}_${REP}.recode.vcf -counts > ./RADseq_RZooRoH/SNPs_DENSITY/Coverage_${SimID}_${NB_RAD_FRAG}_${REP}.txt

	#convert the file in Oxford format
	plink --vcf ./RADseq_RZooRoH/TMP_${SimID}_${NB_RAD_FRAG}_${REP}.recode.vcf --recode oxford --autosome --chr-set 30 --out ./RADseq_RZooRoH/TMP_${SimID}_${NB_RAD_FRAG}_${REP}_OX

	#Open R to do the modelling
	R --vanilla <<EOF

		library(RZooRoH)
		library(doParallel)
		library(foreach)
		cl <- 10
		registerDoParallel(cl)

		# read the OX file
		data_Rohs = zoodata("./RADseq_RZooRoH/TMP_${SimID}_${NB_RAD_FRAG}_${REP}_OX.gen", zformat = "gp")

		#Create 3 HBD  and 1 non-HBD class.es model
		Mod3R <- zoomodel(K=5, krates=c(10,100,1000,10000,10000), err = 2.5e-8)

		#Run the model
		loc_mod3r <- zoorun(Mod3R, data_Rohs, localhbd = TRUE, nT = 10)

		#Extract table with all important information
		loc_table_3R <- loc_mod3r@hbdseg

		#Change indv names so they match simulations + plink ouputs
		loc_table_3R[,1][loc_table_3R[,1] == 1] = 1000
		loc_table_3R[,1][loc_table_3R[,1] != 1000] = loc_table_3R[,1][loc_table_3R[,1] != 1000] - 1

		#Write this file for TFPN Rates later
		write.table(loc_table_3R, "./RADseq_RZooRoH/ROHs_RADseq_${SimID}_${NB_RAD_FRAG}_${REP}.hom", quote = F, col.names = T, row.names = F)

		## FROH

		#create df
		F3R = 1-loc_mod3r@realized[,5]
		individuals = c(1000, seq(1,(loc_mod3r@nind - 1)))
		df_3R = data.frame(individuals, F3R)

		# cbind into output df
		dt_3R <- cbind(MODEL = rep("4R", nrow(df_3R)), SimID=rep("${SimID}" , nrow(df_3R)), df_3R, WINDOWS.NB=rep("${NB_RAD_FRAG}", nrow(df_3R)), REP=rep("${REP}" , nrow(df_3R)))
		colnames(dt_3R) <- c("Model", "SimID" , "individuals" , "FROH_sub" , "WINDOWS.NB" , "REP")

		#Save the data
		write.table(dt_3R, "./RADseq_RZooRoH/RADTMP_${SimID}_${NB_RAD_FRAG}_${REP}_summary_stats.txt", col.names = F , row.names = F , quote = F)

		## ROHs DIST

		#Create ROHs DIST table
		dta_dist <- cbind(rep("4R", nrow(loc_table_3R)), rep("${SimID}" , nrow(loc_table_3R)), rep("${NB_RAD_FRAG}" , nrow(loc_table_3R)), rep("${REP}" , nrow(loc_table_3R)), loc_table_3R)
		colnames(dta_dist) <- c("Model", "SimID" , "WINDOWS.NB" ,"REP" , "individuals" , "chrom" , "start_snp", "end_snp" , "start_pos" , "end_pos" , "number_snp" , "length" , "HBD")

		#Create the output df
		dta_roh = as.data.frame(matrix(nrow=6,ncol=6))
		colnames(dta_roh) = c("SimID","WINDOWS.NB","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")

		#Fill SimuID, CALL_METH, REP, PERC
		dta_roh[,1] = rep("${SimID}",nrow(dta_roh))
		dta_roh[,2] = rep("${NB_RAD_FRAG}",nrow(dta_roh))
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
	write.table(dta_roh, file = "./RADseq_RZooRoH/RADTMP_${SimID}_${NB_RAD_FRAG}_${REP}_ROHsDist.txt", quote = F, row.names = F, col.names = F)

EOF

		#rm TMP VCFs & OX files
		rm ./RADseq_RZooRoH/TMP_${SimID}_${NB_RAD_FRAG}_${REP}.*
		rm ./RADseq_RZooRoH/TMP_${SimID}_${NB_RAD_FRAG}_${REP}_OX.*

}

export -f RAD_Subsampling_ROH_calling_RZooRoH

#THIS IS FOR THE SMALL POPULATION
parallel -j10 ::: RAD_Subsampling_ROH_calling_RZooRoH ::: $(find ./Final_VCF -name "FINAL*.vcf.gz") ::: {1000,3000,5000,7000,9000,10000,20000,30000,40000,50000,60000} ::: {1..100}

#The next line is the corresponding line for the large population !
#parallel -j10 ::: RAD_Subsampling_ROH_calling_RZooRoH ::: $(find ./Final_VCF -name "FINAL*.vcf.gz") ::: {250,500,750,1000,2500,5000,7500,10000,25000,50000,75000,100000,250000,500000,750000,1000000} ::: {1..100}


## Then just merge the SIM replicates into one dataframe

#in case
mkdir -p ./Analyses
mkdir -p ./RADseq_RZooRoH/stats


## FROH

#Create header
echo -e "Model\tSimID\tindividuals\tFROH_sub\tWINDOWS.NB\tREP" > ./RADseq_RZooRoH/stats/Radseq_RZooRoH_FROH.txt
#Add all files
for file in $(find ./RADseq_RZooRoH -name "RADTMP_*_summary_stats.txt"); do cat ${file} >> ./RADseq_RZooRoH/stats/Radseq_RZooRoH_FROH.txt; done
#Rm per SimID file
rm ./RADseq_RZooRoH/RADTMP_*_summary_stats.txt

#Open R for merging
R --vanilla <<EOF

	#Read RAD-seq + WGS data
	dtasub = read.table("./RADseq_RZooRoH/stats/Radseq_RZooRoH_FROH.txt", header = T)
	dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
	dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

	#merge both (without windows nb and REP for WGS)
	dtaFULL = merge(dtasub, dtaTRUE100, by = c("SimID","individuals"))
	dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("SimID","individuals"))
	#Save the data
	write.table(dtaFULL.1, "./Analyses/RADSEQ_RZooRoH_FROH.txt", quote = F, col.names = T, row.names = F)

EOF

## ROHs Dist

#Create header
echo -e "SimID\tWINDOWS.NB\tREP\tROHs_CLASS\tMean_NROH_per_ind\tMean_Tot_Kb_per_ind" > ./Analyses/RADSEQ_RZooRoH_ROHsDistribution.txt
#Add all files
for file in $(find ./RADseq_RZooRoH/ -name "RADTMP_*_ROHsDist.txt"); do cat ${file} >> ./Analyses/RADSEQ_RZooRoH_ROHsDistribution.txt; done
#Rm per SimID file
rm ./RADseq_RZooRoH/RADTMP_*_ROHsDist.txt
