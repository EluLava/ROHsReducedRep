##########################################################################################
######################################### README #########################################
##########################################################################################

THIS PIPELINE IS FOR SIMULATING THE SMALL AND LARGE POPULATIONS; ALL SCRIPTS are made to run in the same directory and will create directories with data at each step

Scripts were launched in SmallPop and LargePop folders respectively !

# 1. SIMULATIONS


10 simulations replicates were performed

## 1.1 Simulating recombination Map with Fregene
Genetic maps are produced in this step and are formated for slim

	1.1.1_geneticMap_simulation.sh --> Will launch Fregene. Input files to simulate humans-like genetic map can be obtained from Fregene website or from Fregene downloaded folder in Fregene/Example/data/

	1.1.2_Formatting_RecMap.R --> Format the recombination map to match slim input requirements


## 1.2 SLiM simulation to record a pedigree (one per simulation replicate --> 10 in my case)
Produces a death_seed.txt and mating_seed.txt per simualtion replicate

A death_seed.txt and mating_seed.txt output files are created to record mating and deaths at each generations of the simulation(s) matching our population characteristics (Ne, inbreeding proportion, etc.)

	1.2_Simulating_Inbreeding_Coeff_nonWF.slim --> Will output matings and deaths at each generation for a population with defined (within script) Ne


## 1.3 Choosing which individuals to subsample (for further analyses) from the pedigree simulated in 1.2 (one per simulation replicate --> 10 in my case)
outputs a list of individuals to use for further analyses

For each simulation (replicate) individuals covering inbreeding range with FPED from 0 to 0.5 (20 each 0.1 if possible) are selected

	1.3_Estimating_Fped_for_ind_Susbampling.R --> select individuals to use for further analyses for each FPED category


## 1.4 SLiM Simulation to apply the pedigree (all chromosomes of one simulation replicate have the same) (one per chromosome per simulation replicate --> 300 in my case).
A Tree Sequence of the whole simulation is outputed

	1.4_launchingSlim2_ApplyingPed.sh --> will launch the SLiM script 1.4_Applying_pedigree_nonWF.slim and run one simulation per chromosome, per replicate. Uses the mating_seed.txt produced 1.2_Simulating_Inbreeding_Coeff_nonWF.slim file name to extract simulation ID (i.e. seed)


## 1.5 Recapitation, Mutation Overlay, Individual subsampling, VCF Writing
Produces a VCF (and a TReeSeq) per chromosome per simulation replicate

	1.5_launching_python_script.sh --> Will launch 1.5_MSprimeRecapitation.py and for each chromosome for each simulation replicate will: recapitate and simplify tree sequence, overlay mutations (and dump a full tree just in case), then subsample individuals (from 1.3) and output a VCF + subsampled tree sequence.


## 1.6 VCF modification
makes VCF data softwares friendly

	1.6.1_change_CHR_ID_ALT_REF_IDs_contig_METADATA.sh --> Changes in VCF: CHROM ($1) from 1 to "true" CHR, ALT ($4) from 0 to A, REF ($5) from 1 to T; add in metadata all contigs IDs; changes tsk_0 to tsk_1000 (for PLINK analyses)

	1.6.2_merging_CHR_VCF.sh --> fuse all CHRs VCFs from the same simulation in one VCF per simulation = FINAL_SimuID.vcf

