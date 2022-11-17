#! /bin/bash

#This script runs on my computer

#List VCFs produced in 1.5
simulations=$(find ./VCFs -name "*.vcf")

Setting_CHR_ID(){

        Simu=$1
        #Getting the Simulation ID but string contains SimID & CHR name
        SimuID=$(basename -s ".vcf" $Simu)
        #GET the CHR
        CHR_ID=$(echo $SimuID | cut -d'_' -f3)
        CHR=$(echo $CHR_ID | cut -c 4-)

        #Getting metadata and storing them in a new file AND change individual tsk_0 name to tsk_1000 and CHR ID and add all contigs
        grep '#' $Simu | awk -v OFS='\t' -v chr=$CHR -v indiv="tsk_1000" 'NR<4{print $0};NR==4{for(i=1;i<=30;i++){print "##contig=<ID="i",length=100000000>"}};NR==5{print $0};END{$10=indiv;print}' > ./VCFs/TMP_${SimuID}.vcf

        #Getting non metadata and changing CHR ID
        grep -v '#' $Simu | awk -v OFS='\t' -v chr=$CHR '{$1=chr};{$4="A"};{$5="T"};{print $0}' >> ./VCFs/TMP_${SimuID}.vcf

        #rm original NO metadata VCF
        rm ${Simu}

}

export -f Setting_CHR_ID
parallel -j 5 ::: Setting_CHR_ID ::: $(for sim in $simulations; do echo $sim; done)
