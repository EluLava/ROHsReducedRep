##########################################################################################
######################################### README #########################################
##########################################################################################

THIS PIPELINE IS FOR SIMULATING THE CATTLE POPULATIONS; ALL SCRIPTS are made to run in the same directory and will create directories with data at each step

This folder only contains scripts for SIMULATIONs steps since the rest of the pipeline is the same for the cattle and small and large populations !

# 1. SIMULATIONS

## 1.1 BURN-IN with MSPRIME from the domestication looong ago, to the creation of the Holstein breed (200 generation ago)

	1.1_MSPRIME_BURNIN.py --> contains the msprime BURN-IN code. From domestication to creation of Holstein (200 G ago). NO MUTATIONS.
	Saves a tree file annotated for SLiM.
	1.1_launching_BURN_IN.sh --> Will launch 1.1_MSPRIME_BURNIN.py script and run cattle simulations (for each 29 CHR).

## 1.2 SLIM CREATING PEDIGREE

	1.2.1_Slim1_launching_ped.sh --> Launch first round simulations with the script: 1.2.1_Cattle_creating_ped.slim
	1.2.2_Cattle_creating_ped.slim --> mating_seed.txt and death_seeds.txt will be recorded and be used later to have same history for all CHR
	Only for 200G (since breed creation)

## 1.3 ADDING THE REAL PED GENERATIONS TO MATING_SEED.txt

	1.3.1_ADDING_GEN_REAL_PED.R --> Artificially add Generations to the real pedigree (for Simulating in SLiM)
	1.3.2_Merging_mating_Real_ped.R --> Real pedigree is appended/merged (and trimmed of useless individuals) with the mating_seed.txt file from 1.2.2_Cattle_creating_ped.slim

## 1.4 CREATING RECOMBINATION MAPS

	1.4_Rending_RecMaps_SLIMFriendly.R --> Uses real *Bos* *taurus* maps and adapts the format for SLiM

## 1.5 SLiM APPLYING PEDIGREES

	1.5.1_Slim2_launching_sim.sh --> Will launch the applying pedigree part: 1.5.1_Cattle_applying_ped.slim
	1.5.1_Cattle_applying_ped.slim --> Uses the mating_seed.txt and death_seed.txt created in Step 2 to apply a pedigree to all CHR of a simulation


## 1.6. ADDING MUTATIONS, SUBSETTING INDIVIDUALS, OUTPUTTING VCFS

	1.6.1_launching_mut_subset.py --> Will launch 1.6.1_MSPRIME_mut_subset.py --> Load SlimTree, add neutral mutations, subset individuals, output VCF
	1.6.2_MonoFilter_Changing_MetaDat.sh --> Will first filter for monomorphic SNPs (bcftools), then edit metadata
	1.6.3_merging_CHR_VCF.sh --> Merging CHRs with bcftools