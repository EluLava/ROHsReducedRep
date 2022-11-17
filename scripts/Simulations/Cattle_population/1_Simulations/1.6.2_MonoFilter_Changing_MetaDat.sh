#! /bin/bash

VCFFiles=$(find ./MSprime_Recap_Mut_VCF -name "*.vcf")

MONO_filtering_and_METADATA_editing(){

	#First argument passed to the script is file
	file=$1

	#extract CHR
	CHR_ID=$(basename -s ".vcf" ${file} | cut -d'_' -f3)
	CHR=$(echo ${CHR_ID} | cut -c 4-)
	#Extract SimuID
	SimuID=$(basename -s ".vcf" ${file} | cut -d'_' -f1)

	#filter on MAC
	bcftools view -c 1:minor -o ./MSprime_Recap_Mut_VCF/${SimuID}_CHR_${CHR}_MONO.vcf --threads 2 ${file}

	#Store CHR_Lengths
	CHR_Lengths=(159000000,137500000,121900000,121300000,121600000,120000000,113200000,114100000,106200000,104700000,107600000,91600000,84600000,85000000,85700000,82000000,75400000,66300000,64300000,72300000,71900000,61700000,52700000,62900000,43100000,51900000,45700000,46500000,51700000)

	#Getting metadata and storing them in a new file AND change individual tsk_0 name to tsk_11000 and CHR ID and add all contigs
	head -n10 ./MSprime_Recap_Mut_VCF/${SimuID}_CHR_${CHR}_MONO.vcf | awk -v OFS='\t' -v chr=${CHR} -v chr_length=${CHR_Lengths} -v indiv="tsk_11000" 'BEGIN{split(chr_length,array,",")};NR<4{print $0};NR==4{for(i=1;i<=29;i++){print "##contig=<ID="i",length="array[i]">"}};NR==5{print $0};END{$10=indiv;print}' > ./MSprime_Recap_Mut_VCF/METADATA_${SimuID}_CHR_${CHR}.vcf
	#Getting non metadata and changing CHR ID
	tail -n +11 ./MSprime_Recap_Mut_VCF/${SimuID}_CHR_${CHR}_MONO.vcf | awk -v OFS='\t' -v chr=$CHR '{$1=chr};{$3=NR};{$4="A"};{$5="T"};{print $0}' >> ./MSprime_Recap_Mut_VCF/METADATA_${SimuID}_CHR_${CHR}.vcf


}

export -f MONO_filtering_and_METADATA_editing

parallel -j29 ::: MONO_filtering_and_METADATA_editing ::: $(for sim in ${VCFFiles}; do echo $sim; done)