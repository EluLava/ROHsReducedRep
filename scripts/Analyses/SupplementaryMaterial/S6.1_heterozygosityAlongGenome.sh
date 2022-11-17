VCF_HetAlongGenome() {

	#FIrst argument is VCF file
	FILE=$1
	#Get the SmuID
	SimuID=$(basename -s ".vcf.gz" $FILE | cut -d '_' -f2)

	#Create directory if needed
	mkdir -p ./hetero_along_genome/

	R --vanilla <<EOF

		#library
		library("vcfR")
		#Read the VCF
		dta = read.vcfR("${FILE}")

		#Extract GT elements
		gt = extract_gt_tidy(dta, format_fields = c("GT"))[,c(1,2,3)]
		gt = as.data.frame(gt)
		gt = cbind(vector(length = nrow(gt)), gt)
		colnames(gt) = c("CHR", "POS", "INDIV", "HET")

		#individuals vector
		individuals = unique(gt[,3])

		#Extracting the positions from the vcf
		gt[,1] = rep(dta@fix[,1],length(individuals))
		gt[,2] = rep(dta@fix[,2],length(individuals))

		#CHanging the structure
		gt[,1] = factor(gt[,1], levels = seq(1,30))
		gt[,2] = as.numeric(gt[,2])
		gt[,3] = factor(gt[,3], levels = individuals)
		gt[,4] = as.character(gt[,4])

		#REMOVE VCF
		remove(dta)

		#Calculating heterozygosity dosage
		gt[gt == "0|0" | gt == "1|1"] <- 0
		gt[gt == "0|1" | gt == "1|0"] <- 1

		#Change HET in numeric
		gt[,4] = as.numeric(gt[,4])

		#empty dataframe
		data_wind = as.data.frame(matrix(ncol = 5,nrow = 0))
		colnames(data_wind) = c("SimuID", "CHR", "POS", "INDIV", "HET")
		data_wind[,3] = as.numeric(data_wind[,3])
		data_wind[,5] = as.numeric(data_wind[,5])

		print("Starting sliding windows")

		for (chr in levels(gt[,1])) {
  			#Tmp data frame with only CHR
  			gt_tmp = gt[gt[,1] == chr,]
  			#Last position in the contig
  			max_pos = max(gt_tmp[,2])
  			#non-overlapping sliding-windows each 50kb, last one smaller until the last pos
  			windows = c(seq(0,max_pos,50000), max_pos)
  
  			# LOOPING THROUGH SLIDING windows
  			for (wind in 1:(length(windows)-1)) {
  				#Tmp dataframe with only slinding window
    			tmp = gt_tmp[gt_tmp[,2]>windows[wind] & gt_tmp[,2] <= windows[wind+1],]
    			#Creating output dataframe. CHR = chr, POS = mid of sliding window INDIV et HET)
    			tmp_out = as.data.frame(cbind(CHR = rep(chr, length(individuals)), POS = as.numeric(rep((windows[wind] + ((windows[wind+1] - windows[wind])/2)), length(individuals))), INDIV = individuals, HET = vector(length = length(individuals))))
    			tmp_out[,2] = as.numeric(as.character(tmp_out[,2]))
    			tmp_out[,4] = as.numeric(as.character(tmp_out[,4]))
    			tmp_out[,3] = individuals

    			#Looping through individuals to calculate HET
    			for (indiv in 1:length(individuals)) {
    				tmp_out[,4][tmp_out[,3] == individuals[indiv]] = mean(tmp[,4][tmp[,3] == individuals[indiv]], na.rm = T)
    			}

    		#Add the simuID
    		tmp_out = cbind(rep("${SimuID}", nrow(tmp_out)), tmp_out)
    		colnames(tmp_out)[1] = "SimuID"

    		#Fusing data in one dataframe
    		data_wind = rbind(data_wind,tmp_out)
    		
    		}

    		print(paste0(chr, " has been processed with sliding-window!"))
		
		}

		write.table(data_wind, file = "./hetero_along_genome/TMP_${SimuID}_heterozygosity_sliding_windows_50kb.txt", quote = F, row.names = F, col.names = T)


EOF

}

export -f VCF_HetAlongGenome
parallel -j10 ::: VCF_HetAlongGenome ::: $(find ./Final_VCF -name "FINAL*.vcf.gz")


FINAL_FILES=$(find ./hetero_along_genome/ -name "TMP_*_heterozygosity_sliding_windows_50kb.txt")

#Get the header from one file
awk 'NR == 1 {print $0}' ./hetero_along_genome/TMP_Sim1668789343933_RADSUB_heterozygosity_sliding_windows_50kb.txt > ./hetero_along_genome/HETERO_ALONG_GENOME.txt

for file in $FINAL_FILES; do awk 'NR > 1 {print $0}' $file >> ./hetero_along_genome/HETERO_ALONG_GENOME.txt; done

#rm TMP files
rm ./hetero_along_genome/TMP_*_heterozygosity_sliding_windows_50kb.txt

