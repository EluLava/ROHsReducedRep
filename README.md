# ROHsReducedRep

This repository contains all the necessary codes for this paper: https://doi.org/10.1101/2022.08.26.505374

## ABSTRACT

Genomic measures of inbreeding based on Identical-by-Descent (IBD) segments are increasingly used to measure inbreeding and mostly estimated on SNP arrays and whole-genome-sequencing (WGS) data. However, some softwares recurrently used for their estimation assume that genomic positions which have not been genotyped are non-variant. This might be true for WGS data, but not for reduced genomic representations and can lead to spurious IBD segments estimation. Here we simulated the outputs of WGS, two SNP arrays of different sizes and RAD-sequencing for three populations with different sizes and histories. We compare the results of IBD segments estimation with two softwares: runs of homozygosity (ROHs) estimated with <i>PLINK</i> and Homozygous-by-descent (HBD) segments estimated with <i>RZooRoH</i>. We demonstrate that to obtain meaningful estimates of inbreeding, <i>RZooRoH</i> requires a SNPs density eleven times smaller compared to <i>PLINK</i>: ranks of inbreeding coefficients were conserved among individuals above 22 SNPs/Mb for <i>PLINK</i> and 2 SNPs/Mb for <i>RZooRoH</i>. We also show that in populations with simple demographic histories, ROHs and HBD segments distributions are correctly estimated with both SNP arrays and WGS. <i>PLINK</i> correctly estimated ROHs distributions with SNP densities above 22 SNPs/Mb while <i>RZooRoH</i> correctly estimated HBD segments distribution with SNPs densities above 11 SNPs/Mb. However, in a population with a more complex demographic history, <i>RZooRoH</i> resulted in better IBD segments distributions estimation compared to <i>PLINK</i> even with WGS data. Consequently, we advise researchers to use either methods relying on excess of homozygosity averaged across SNPs or model-based HBD segments calling methods for inbreeding estimations.

## INFORMATION

- The **scripts** folder contains all the scripts needed to run all the analyses (simulations + downstream analyses) presented in the paper.

- All the scripts for the simulations part are in the Simulations folder

Simulations scripts are the same for 2 populations: the "small" and "large" population (EXCEPT THE Ne PARAMETER WHICH VARIES)

Details about each script and in which order they should be run are available in README_Small_LargePops.md and README_CattlePop.md in scripts/Simulations/Small_Large_populations and scripts/Simulations/Cattle_population folder respectively.

- Downstream analyses scripts are the same for the three populations

Details about each script and in which order they should be run are available in README_Analyses.md in /scripts/Analyses/ folder.
