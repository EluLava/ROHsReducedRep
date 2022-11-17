#! /bin/bash

#Create output directory
mkdir -p ./Final_VCFs

#Get the unique sim IDs
files=$(find ./MSprime_Recap_Mut_VCF -name "METADATA_*.vcf")
SimuIDs=$(basename -s ".vcf" $files | cut -d'_' -f2 | sort | uniq)
#Iterating through Sim IDs
for sim in $SimuIDs; do
	
#Get the file list
Input_files=$(find ./MSprime_Recap_Mut_VCF -name "METADATA_${sim}*.vcf" | sort --version-sort)
#Print them in a file .list
for file in $Input_files;
do
echo $file >> ./Final_VCFs/INPUT_${sim}.list
done

#Use bcftools to fuse files
bcftools concat -f ./Final_VCFs/INPUT_${sim}.list -o ./Final_VCFs/${sim}_ALL_CHRs.vcf

#RM list file
rm ./Final_VCFs/INPUT_${sim}.list

done
