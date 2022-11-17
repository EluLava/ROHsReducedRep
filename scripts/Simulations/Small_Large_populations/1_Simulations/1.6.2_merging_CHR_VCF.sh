#! /bin/bash

#Get the unique sim IDs (from metadata changed VCF)
files=$(find . -name "TMP_*.vcf")
SimuIDs=$(basename -s ".vcf" $files | cut -d'_' -f2 | sort | uniq)

#Create output directory if needed
mkdir -p ./Final_VCF

#Iterating through Sim IDs
for sim in $SimuIDs; do

        #Get the file list
        Input_files=$(find . -name "TMP_${sim}*.vcf" | sort --version-sort)
        #Print them in a file .list
        for file in $Input_files;
        do
          	echo $file >> ./VCFs/INPUT_${sim}.list
        done

        #Use bcftools to merge the VCFs files and rm monomorphic sites (were monomorphic in my subsampled indvs not the entire pop)
        bcftools concat -f ./VCFs/INPUT_${sim}.list | bcftools view -c 1:minor -o ./Final_VCF/TMP_${sim}.vcf

        #RM list file
        rm ./VCFs/INPUT_${sim}.list

        #Change SNPs ID
        #Grep metadata
        grep '#' ./Final_VCF/TMP_${sim}.vcf > ./Final_VCF/FINAL_${sim}.vcf
        #Grep non metadata and change SNPID
        grep -v '#' ./Final_VCF/TMP_${sim}.vcf | awk -v OFS='\t' '{$3=NR};{print $0}' >> ./Final_VCF/FINAL_${sim}.vcf

        #Rm all TMP VCF
        rm ./VCFs/TMP_${sim}_*.vcf
        rm ./Final_VCF/TMP_${sim}.vcf

        #Zip & Index
        bgzip ./Final_VCF/FINAL_${sim}.vcf
        bcftools index ./Final_VCF/FINAL_${sim}.vcf.gz

done
