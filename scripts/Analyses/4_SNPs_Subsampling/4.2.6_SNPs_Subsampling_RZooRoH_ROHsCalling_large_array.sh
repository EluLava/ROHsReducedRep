#!/bin/bash

SUB_ROH_CALL_FROH_ROHsDist_RZooRoH() {

	#SNP_CHR_POS_${SimID}_${REP}.csv from script 3.2.4 is the first argument passed to the function
	FILE=$1

	#Extract SimID & REP (1) info
	SimID=$(basename -s ".csv" $FILE | cut -d'_' -f4)
	REP=$(basename -s ".csv" $FILE | cut -d'_' -f5)
	
	#Subsample SNPs to create a new VCF
	bcftools view -R ./SNP_ARRAYs/Large_ARRAY/SNP_CHR_POS_${SimID}_${REP}.csv -O z -o ./SNP_ARRAYs/Large_ARRAY/VCF_SNPs_${SimID}_${REP}.vcf.gz ./Final_VCF/FINAL_${SimID}.vcf.gz

	#From VCF to OX
	plink --vcf ./SNP_ARRAYs/Large_ARRAY/VCF_SNPs_${SimID}_${REP}.vcf.gz --recode oxford --autosome --chr-set 30 --out ./SNP_ARRAYs/Large_ARRAY/TMP_${SimID}_${REP}_OX

	#Open R to do the modelling
	R --vanilla <<EOF

		library(RZooRoH)
		library(doParallel)
		library(foreach)
		cl <- 10
		registerDoParallel(cl)

		# read the OX file
		data_Rohs = zoodata("./SNP_ARRAYs/Large_ARRAY/TMP_${SimID}_${REP}_OX.gen", zformat = "gp")

		#Create 3 HBD  and 1 non-HBD class.es model
		Mod3R <- zoomodel(K=5, krates=c(10,100,1000,10000,10000), err = 5e-5)

		#Run the model
		loc_mod3r <- zoorun(Mod3R, data_Rohs, localhbd = TRUE, nT = 10)

		#Extract table with all important information
		loc_table_3R <- loc_mod3r@hbdseg

		#Change indv names so they match simulations + plink ouputs
		loc_table_3R[,1][loc_table_3R[,1] == 1] = 1000
		loc_table_3R[,1][loc_table_3R[,1] != 1000] = loc_table_3R[,1][loc_table_3R[,1] != 1000] - 1

		#Write this file fro TFPN Rates later
		write.table(loc_table_3R, "./SNP_ARRAYs/Large_ARRAY/ROHs_LargeArray_${SimID}_${REP}.hom", quote = F, col.names = T, row.names = F)


		## FROH

		#create df
		F3R = 1-loc_mod3r@realized[,5]
		individuals = c(1000, seq(1,(loc_mod3r@nind - 1)))
		df_3R = data.frame(individuals, F3R)

		# cbind into output df
		dt_3R <- cbind(MODEL = rep("4R", nrow(df_3R)), SimID=rep("${SimID}" , nrow(df_3R)), df_3R, REP=rep("${REP}" , nrow(df_3R)))
		colnames(dt_3R) <- c("Model", "SimID" , "individuals" , "FROH_LargeArray" , "REP")

		#Save the data
		write.table(dt_3R, "./SNP_ARRAYs/Large_ARRAY/FROH_RZooRoH_LargeArray_${SimID}_${REP}.txt", col.names = F , row.names = F , quote = F)

		## ROHs DIST

		#Create ROHs DIST table
		dta_dist <- cbind(rep("3R", nrow(loc_table_3R)), rep("${SimID}" , nrow(loc_table_3R)), rep("LARGE_ARRAY" , nrow(loc_table_3R)), rep("${REP}" , nrow(loc_table_3R)), loc_table_3R)
		colnames(dta_dist) <- c("Model", "SimID" , "WINDOWS.NB" ,"REP" , "individuals" , "chrom" , "start_snp", "end_snp" , "start_pos" , "end_pos" , "number_snp" , "length" , "HBD")

		#Create the output df
		dta_roh = as.data.frame(matrix(nrow=6,ncol=6))
		colnames(dta_roh) = c("SimID","WINDOWS.NB","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")

		#Fill SimuID, CALL_METH, REP, PERC
		dta_roh[,1] = rep("${SimID}",nrow(dta_roh))
		dta_roh[,2] = rep("LARGE_ARRAY",nrow(dta_roh))
		dta_roh[,3] = rep("${REP}",nrow(dta_roh))
		dta_roh[,4] = c(1:6)

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
			dta_roh[roh_class,6] = sum(x[3][x[2] == roh_class])/(sort(unique(dta_dist.2[,5]))[length(unique(dta_dist.2[,5]))-1])
			dta_roh[roh_class,5] = sum(y[3][y[2] == roh_class])/(sort(unique(dta_dist.2[,5]))[length(unique(dta_dist.2[,5]))-1])
		}

		#write output
		write.table(dta_roh, file = "./SNP_ARRAYs/Large_ARRAY/ROHsDist_RZooRoH_LargeArray_${SimID}_${REP}.txt", quote = F, row.names = F, col.names = F)

EOF

	#rm TMP VCFs & OX files
	rm ./SNP_ARRAYs/Large_ARRAY/TMP_${SimID}_${REP}_OX.*
	rm ./SNP_ARRAYs/Large_ARRAY/VCF_SNPs_${SimID}_${REP}.*

}

export -f SUB_ROH_CALL_FROH_ROHsDist_RZooRoH

#launch the function on all the Simulations replicates
parallel -j10 ::: SUB_ROH_CALL_FROH_ROHsDist_RZooRoH ::: $(find ./SNP_ARRAYs/Large_ARRAY/ -name "SNP_CHR_POS_*.csv")


### Then just fuse all Simulations in one df

## FROH

#Create df in case does not exists
mkdir -p ./Analyses

#Create header
echo -e "Model\tSimID\tindividuals\tFROH_LargeArray\tREP" > ./SNP_ARRAYs/Large_ARRAY/LargeArray_RZooRoH_FROH.txt

#Loop through file to fill the data
for file in $(find ./SNP_ARRAYs/Large_ARRAY/ -name "FROH_RZooRoH_LargeArray_*.txt"); do cat ${file} >> ./SNP_ARRAYs/Large_ARRAY/LargeArray_RZooRoH_FROH.txt; done

#rm the sub files
rm ./SNP_ARRAYs/Large_ARRAY/FROH_RZooRoH_LargeArray_*.txt

#Open R for merging
R --vanilla <<EOF

	#Read ARRAY + TRUE IBD data
	dtasub = read.table("./SNP_ARRAYs/Large_ARRAY/LargeArray_RZooRoH_FROH.txt", header = T)
	dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
	dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

	#merge both (without windows nb and REP for WGS)
	dtaFULL = merge(dtasub, dtaTRUE100, by = c("SimID","individuals"))
	dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("SimID","individuals"))
	#Save the data
	write.table(dtaFULL.1, "./Analyses/LargeArray_RZooRoH_FROH.txt", quote = F, col.names = T, row.names = F)

EOF


## ROHs Dist

#Create header
echo -e "SimuID\tSUB_TECH\tREP\tROHs_CLASS\tMean_NROH_per_ind\tMean_Tot_Kb_per_ind" > ./Analyses/LargeArray_RZooRoH_ROHsDistributions.txt

#Loop through file to fill the data
for file in $(find ./SNP_ARRAYs/Large_ARRAY/ -name "ROHsDist_RZooRoH_LargeArray_*.txt"); do cat ${file} >> ./Analyses/LargeArray_RZooRoH_ROHsDistributions.txt; done

#rm the sub files
rm ./SNP_ARRAYs/Large_ARRAY/ROHsDist_RZooRoH_LargeArray_*.txt


