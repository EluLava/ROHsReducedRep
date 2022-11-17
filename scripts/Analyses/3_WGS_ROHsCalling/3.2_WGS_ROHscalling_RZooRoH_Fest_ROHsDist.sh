#! /bin/bash

#Define the function
WGS_ROHsCal_RZooRoH() {

	#First argument passed to the function is the WGS VCF
	VCF=$1
	#Extract SimID from VCF filename
	Sim_ID=$(basename -s ".vcf.gz" ${VCF} | cut -d'_' -f2)
	
	#Create output directory for ROHs called with RZooRoH if does not exists
	mkdir -p ./WGS_ROHs_RZooRoH
	
	#convert the file in Oxford format
	plink --vcf ${VCF} --recode oxford --autosome --chr-set 30 --out ./WGS_ROHs_RZooRoH/WGS_${Sim_ID}_OX
	
	#Open R for ROHs calling etc.
	R --vanilla <<EOF

		library(RZooRoH)
		library(doParallel)
		library(foreach)
		cl <- 5
		registerDoParallel(cl)

		# read the OX file
		data_Rohs = zoodata("./WGS_ROHs_RZooRoH/WGS_${Sim_ID}_OX.gen", zformat = "gp")

		#Create 3 HBD  and 1 non-HBD class.es model
		Mod3R <- zoomodel(K=5, krates=c(10,100,1000,10000,10000), err = 5e-5)

		#Run the model
		loc_mod3r <- zoorun(Mod3R, data_Rohs, localhbd = TRUE, nT = cl)

		#Extract table with all important information
		loc_table_3R <- loc_mod3r@hbdseg

		#Change indv names so they match simulations + plink ouputs
		loc_table_3R[,1][loc_table_3R[,1] == 1] = 1000
		loc_table_3R[,1][loc_table_3R[,1] != 1000] = loc_table_3R[,1][loc_table_3R[,1] != 1000] - 1

		#Write this file fro TFPN Rates later
		write.table(loc_table_3R, "./WGS_ROHs_RZooRoH/ROHs_WGS_${Sim_ID}.hom", quote = F, col.names = T, row.names = F)

		## FROH

		#create df
		F3R = 1-loc_mod3r@realized[,5]
		individuals = c(1000, seq(1,(loc_mod3r@nind - 1)))
		df_3R = data.frame(individuals, F3R)

		# cbind into output df
		dt_3R <- cbind(MODEL = rep("4R", nrow(df_3R)), SimID=rep("${Sim_ID}" , nrow(df_3R)), df_3R, WINDOWS.NB=rep("WGS", nrow(df_3R)), REP=rep("1" , nrow(df_3R)))
		colnames(dt_3R) <- c("Model", "SimID" , "individuals" , "Froh_WGS" , "WINDOWS.NB" , "REP")

		#Save the data
		write.table(dt_3R, "./WGS_ROHs_RZooRoH/WGS_RZooRoH_${Sim_ID}_summary_stats.txt", col.names = F , row.names = F , quote = F)

		## ROHs DIST

		#Create ROHs DIST table
		dta_dist <- cbind(rep("4R", nrow(loc_table_3R)), rep("${Sim_ID}" , nrow(loc_table_3R)), rep("WGS" , nrow(loc_table_3R)), rep("1" , nrow(loc_table_3R)), loc_table_3R)
		colnames(dta_dist) <- c("Model", "SimID" , "WINDOWS.NB" ,"REP" , "individuals" , "chrom" , "start_snp", "end_snp" , "start_pos" , "end_pos" , "number_snp" , "length" , "HBD")

		#Create the output df
		dta_roh = as.data.frame(matrix(nrow=6,ncol=6))
		colnames(dta_roh) = c("SimID","WINDOWS.NB","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")

		#Fill SimuID, CALL_METH, REP, PERC
		dta_roh[,1] = rep("${Sim_ID}",nrow(dta_roh))
		dta_roh[,2] = rep("WGS",nrow(dta_roh))
		dta_roh[,3] = rep("1",nrow(dta_roh))
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
			dta_roh[roh_class,6] = sum(x[3][x[2] == roh_class])/(length(individuals))
			dta_roh[roh_class,5] = sum(y[3][y[2] == roh_class])/(length(individuals))
		}

		#write output
		write.table(dta_roh, file = "./WGS_ROHs_RZooRoH/WGS_RZooRoH_${Sim_ID}_ROHsDist.txt", quote = F, row.names = F, col.names = F)

EOF

	#rm OX files
	rm ./WGS_ROHs_RZooRoH/WGS_${Sim_ID}_OX.*

}

export -f WGS_ROHsCal_RZooRoH

#Loop through simualations replicates
for f in $(find ./Final_VCF -name "FINAL*.vcf.gz"); do WGS_ROHsCal_RZooRoH ${f}; done

# Then just merge the SIM replicates into one dataframe

mkdir -p ./Analyses/

## FROH
#Create header
echo -e "Model\tSimID\tindividuals\tFroh_WGS\tWINDOWS.NB\tREP" > ./WGS_ROHs_RZooRoH/WGS_RZooRoH_FROH.txt
#Add all files
for file in $(find ./WGS_ROHs_RZooRoH -name "WGS_RZooRoH_*_summary_stats.txt"); do cat ${file} >> ./WGS_ROHs_RZooRoH/WGS_RZooRoH_FROH.txt; done
#Rm per SimID file
rm ./WGS_ROHs_RZooRoH/WGS_RZooRoH_*_summary_stats.txt

#Open R for merging
R --vanilla <<EOF

	#Read RAD-seq + WGS data
	dtasub = read.table("./WGS_ROHs_RZooRoH/WGS_RZooRoH_FROH.txt", header = T)
	dtaTRUE100 = read.table("./TRUE_IBD/100GEN/stats/F_TRUE_IBD_100GEN.txt", header = T)
	dtaTRUE1000 = read.table("./TRUE_IBD/1000GEN/stats/F_TRUE_IBD_1000GEN.txt", header = T)

	#merge both (without windows nb and REP for WGS)
	dtaFULL = merge(dtasub, dtaTRUE100, by = c("SimID","individuals"))
	dtaFULL.1 = merge(dtaFULL, dtaTRUE1000, by = c("SimID","individuals"))
	#Save the data
	write.table(dtaFULL.1, "./Analyses/WGS_RZooRoH_FROH.txt", quote = F, col.names = T, row.names = F)

EOF

## ROHs Dist
#Create header
echo -e "SimID\tWINDOWS.NB\tREP\tROHs_CLASS\tMean_NROH_per_ind\tMean_Tot_Kb_per_ind" > ./Analyses/WGS_RZooRoH_ROHsDistribution.txt
#Add all files
for file in $(find ./WGS_ROHs_RZooRoH/ -name "WGS_RZooRoH_*_ROHsDist.txt"); do cat ${file} >> ./Analyses/WGS_RZooRoH_ROHsDistribution.txt; done
#Rm per SimID file
rm ./WGS_ROHs_RZooRoH/WGS_RZooRoH_*_ROHsDist.txt

