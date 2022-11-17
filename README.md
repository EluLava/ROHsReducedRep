# ROHsReducedRep

This repository contains all the necessary codes for this paper: https://doi.org/10.1101/2022.08.26.505374

## ABSTRACT

Genomic measures of inbreeding based on Identical-by-Descent (IBD) segments are increasingly used to measure individual inbreeding and mostly estimated on SNPs-arrays and whole-genome-sequencing (WGS) data. However, some softwares recurrently used for their estimation assume that genomic positions which have not been genotyped are non-variant. This might be true for WGS data, but not for reduced genomic representations and can lead to spurious IBD segments estimation. In this project, we simulated the outputs of WGS, two SNPs arrays of different sizes and RAD-sequencing for three populations with different sizes and histories. We compare the results of IBD segments estimation with two softwares: runs of homozygosity (ROHs) estimated with PLINK and Homozygous-by-descent (HBD) segments estimated with RZooRoH. We demonstrate that to obtain meaningful estimates of inbreeding coefficients, RZooRoH requires fraction of genomes eleven times smaller compared to PLINK: ranks of inbreeding coefficients FROH and FHBD were conserved among individuals above 22 SNPs/Mb for PLINK and 2 SNPs/Mb for RZooRoH. We also show that in population with simple demographic histories, ROHs and HBD segments distributions can be correctly estimated with both SNP arrays and WGS and both softwares. Indeed, PLINK correctly estimated ROHs distributions with SNP densities above 22NPSs/Mb while RZooRoH correctly estimated HBD segments distribution with SNPs densities above 11 SNPs/Mb. However, in a population with a more complex demographic history, RZooRoH resulted in better IBD segments distributions estimation compared to PLINK. Consequently, we advise researchers to use either SNPs-based or model-based HBD segments calling methods for inbreeding estimations.

## INFORMATION

- The **scripts** folder contains all the scripts needed to run all the analyses (simulations + downstream analyses) presented in the paper.

- All the scripts for the simulations part are in the Simulations folder

Simulations scripts are the same for 2 populations: the "small" and "large" population (EXCEPT THE Ne PARAMETER WHICH VARIES)

Details about each script and in which order they should be run are available in README_Small_LargePops.md and README_CattlePop.md in scripts/Simulations/Small_Large_populations and scripts/Simulations/Cattle_population folder respectively.

- Downstream analyses scripts are the same for the three populations

Details about each script and in which order they should be run are available in README_Analyses.md in /scripts/Analyses/ folder.
