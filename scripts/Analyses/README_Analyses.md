##########################################################################################
######################################### README #########################################
##########################################################################################

THIS PIPELINE IS FOR ALL THE POPULATIONS; ALL SCRIPTS are made to run in the same directory and will create directories with data at each step

To change population, the only modifications needed are genome size when estimating FROH and True/False Positive/Negative Rates !

Scripts were launched in SmallPop, LargePop and CattlePop folders !

(SM) indicates that the script·s is·are for data shown only in Supplementary Material


# 2. TRUE IBD segments (from simulations)
outputs true IBD segments (plus true farction of genome within them = true F and true IBD segments distrbutions)

	2.1_Estimating_TRUE_IBD_Segments_from_TreeSeq.sh --> uses python and TreeSequences from recapitation script (1.5_launching_python_script.sh) to estimate TRUE IBD segments !  In the mansucript this was done for TRUE IBD segments from 100 and 1000 generationsa go. This scripts runs for 100 GEN ago. but can be easily modified to run for different values ! The value AND OUTPUT directory must be updated !
	
	2.2_TRUE_FIBD_Est_and_True_IBD_Segments_Distributions.sh --> Estimate fraction of genome within IBD segments + IBD segments lengths distributions.     




# 3. WGS ROHs calling PLINK and RZooRoH


## 3.1 PLINK ROHs calling, FROH quantification and ROHs distributions

	3.1_WGS_ROHscalling_PLINK_Fest_ROHsDist.sh --> Uses plink to call ROHs (min size: 100KB) from WGS VCFs, Uses R to estimate FROH (sum of ROHs/Genome length) and ROHs distributions; Also estimates FHOM (for comparison in supp material)


## 3.2 RZooRoH ROHs calling, FROH quantification and ROHs distributions

	3.2_WGS_ROHscalling_RZooRoH_Fest_ROHsDist.sh --> uses RZooRoH to call HBD segments, estimate F and HBD distributions on WGS VCFs




# 4. SNPs subsampling (reduced genomic representations), ROHs and HBD calling (PLINK and RZooRoH) (and SNP density estimation for RAD-seq)


100 subsampling replicates were performed for RAD-sequencing subsampling

## 4.1 RAD-sequencing-like subsampling (PLINK and RZooRoH separately because the fraction of genome subsampled are different) (genome.bed file containing all CHRs regions is needed) (ONLY SMALL AND LARGE POPULATIONS)

	4.1.1_RADseq_SUBSAMPLING_PLINK_ROHs_calling.sh --> Randomly select RAD fragments, subsample SNPs, call ROHs with PLINK (min 100KB), estimate F and ROHs distributions

	4.1.2_RADseq_SUBSAMPLING_RZooRoH_ROHs_calling.sh --> Randomly select RAD fragments, subsample SNPs, call HBD segments with RZooRoH, estimate F and ROHs distributions

## 4.2 SNP-arrays-like subsampling

	#small array (window size = median distance between SNPs in the Illumina BovineSNP 50 beadchip = 40KB)

	4.2.1_SNPs_selection_small_array.sh --> Select SNPs for the small array: Filter on MAF 0.05, extract SNPs, estimate MAF, loop through genome to select 1 SNP (with higher MAF) per window

	4.2.2_SNPs_Subsampling_PLINK_ROHsCalling_small_array.sh --> Call ROHs with PLINK in the small array: Uses SNPs selected in 4.2.1, subsample VCF, call ROHs with PLINK (same parameters), estimate F and ROHs distributions

	4.2.3_SNPs_Subsampling_RZooRoH_ROHsCalling_small_array.sh --> Call HBD segments with RZooRoH in the small array: Uses SNPs selected in 4.2.1, subsample VCF, call ROHs with RZooRoH, estimate F and ROHs distributions
	

	#large array (window size = median distance between SNPs in the Illumina BovineSNP 50 beadchip = 3KB)

	4.2.4_SNPs_selection_large_array.sh --> Select SNPs for the small array: Filter on MAF 0.05, extract SNPs, estimate MAF, loop through genome to select 1 SNP (with higher MAF) per window

	4.2.5_SNPs_Subsampling_PLINK_ROHsCalling_large_array.sh --> Call ROHs with PLINK in the large array: Uses SNPs selected in 3.2.4, subsample VCF, call ROHs with PLINK (same parameters), estimate F and ROHs distributions

	4.2.6_SNPs_Subsampling_RZooRoH_ROHsCalling_large_array.sh --> Call HBD segments with RZooRoH in the large array: Uses SNPs selected in 3.2.4, subsample VCF, call ROHs with RZooRoH, estimate F and ROHs distributions


## 4.3 RANDOM SNPs subsampling (SM) (ONLY SMALL AND LARGE POPULATIONS)

	S4.3.1_RANDOM_SNPs_PERCENTAGE_SUBSAMPLING_PLINK_ROHs_calling.sh --> Randomly subsample SNPs (with various percentages), call ROHs with PLINK AND HDB segments with RZooROH, estimate F and distributions


## 4.4 Estimating SNPs densities and correlation between TRUE IBD and RAD (SM)

	S4.4.1_RADseq_SNPdensity_Correlations_PLINK.sh --> Subsample SNPs, estimate mean SNP density in reduced dataset, call ROHs, estimate correlation between RAD + TRUE IBD & slope and intercept for lm(FROH RAD ~ FROH TRUE IBD) + estimates difference between F TRUE IBD and FRAD

	S4.4.2_RADseq_SNPdensity_Correlations_RZooRoH.sh --> Subsample SNPs, estimate mean SNP density in reduced dataset, call HBD segments, estimate correlation between RAD + TRUE IBD & slope and intercept for lm(FROH RAD ~ FROH TRUE IBD) + estimates difference between F TRUE IBD and FRAD




# 5. True/False Positive/Negative Rates


## 5.1 WGS

	5.1.1_TFPN_Rates_WGS_PLINK.sh --> creates bedfiles with genomic regions within ROHs per INDIVIDUAL for WGS for PLINK and compares with TRUE IBD bedfiles

	5.1.2_TFPN_Rates_WGS_RZooRoH.sh --> creates bedfiles with genomic regions within HBD segments per INDIVIDUAL for WGS for RZooRoH and compare with TRUE IBD bedfiles

## 5.2 RAD-sequencing (ONLY SMALL AND LARGE POPULATIONS)

	5.2.1_TFPN_Rates_Radseq_PLINK.sh --> creates bedfiles with genomic regions within ROHs per INDIVIDUAL for RAD sequencing and compare these bedfiles with TRUE IBD bedfiles
	
	5.2.2_TFPN_Rates_Radseq_RZooRoH.sh --> creates bedfiles with genomic regions within HBD segments per INDIVIDUAL for RAD sequencing and compare these bedfiles with TRUE IBD bedfiles

## 5.3 SNPs-Arrays

	5.3.1_TFPN_Rates_SmallArray_PLINK.sh --> creates bedfiles with genomic regions within ROHs per INDIVIDUAL for Small Array and compare these bedfiles with TRUE IBD bedfiles
	
	5.3.2_TFPN_Rates_SmallArray_RZooRoH.sh --> creates bedfiles with genomic regions within HBD segments per INDIVIDUAL for Small Array and compare these bedfiles with TRUE IBD bedfiles

	5.3.3_TFPN_Rates_LargeArray_PLINK.sh --> creates bedfiles with genomic regions within ROHs per INDIVIDUAL for Large Array and compare these bedfiles with TRUE IBD bedfiles
	
	5.3.4_TFPN_Rates_LargeArray_RZooRoH.sh --> creates bedfiles with genomic regions within HBD segments per INDIVIDUAL for Large Array and compare these bedfiles with TRUE IBD bedfiles


# 6. Heterozygosity along the genome (SM) (ONLY SMALL AND LARGE POPULATIONS)

	S6.1_heterozygosityAlongGenome.sh --> Estimates Heterozygosity along the genome for each INDV from each simulation with non-overlapping sliding-windows (50kb)


# 7. Plotting (main figures)

	7_Main_plotting.R --> Creates all the plots from main text ! Only the data, illustrations etc. were added in Illustrator: in /scripts/Plotting/ folder


# 8. Plotting (SM)

	8_Supplementary_plotting.R --> Creates all the plots from supplementary material ! Only the data, illustrations etc. were added in Illustrator: in /scripts/Plotting/ folder


