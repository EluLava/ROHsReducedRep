#######################################################################

# RATES

#Create the function
INDIV_Rates_TRUE_IBD() {

#First argument passed to the function = TRUE IBD bedfile !
BEDFILE=$1

#Extract Simulation
SimuID=$(basename -s ".bed" ${BEDFILE} | cut -d'_' -f3)

INDIV=$(basename -s ".bed" ${BEDFILE} | cut -d'_' -f5)

#Get corresponding WGS PLINK file
WGS_FILE=./TFPN_Rates_PLINK/WGS/bedfiles/WGS_BED_${SimuID}_Ind_${INDIV}.bed

#Estimate the TP rate
TP=$(bedtools intersect -a ${BEDFILE} -b ${WGS_FILE} -g ./genome.bed | awk -v TP=0 -v genlength=3000000000 '{TP = TP + ($3 - $2)};END{print TP/genlength}')

#TN rate with complement function from bedtools, but only works with one single file --> I need to "merge" WGS and RAd file first, then complement the result
TN=$(cat ${BEDFILE} ${WGS_FILE} | sort -k1,1n -k2,2n | bedtools complement -i stdin -g ./genome.bed | awk -v TN=0 -v genlength=3000000000 '{TN = TN + ($3 - $2)};END{print TN/genlength}')

#Estimate the FP rate
FN=$(bedtools subtract -a ${BEDFILE} -b ${WGS_FILE} -g ./genome.bed | awk -v FN=0 -v genlength=3000000000 '{FN = FN + ($3 - $2)};END{print FN/genlength}')

#Estimate the TP rate
FP=$(bedtools subtract -a ${WGS_FILE} -b ${BEDFILE} -g ./genome.bed | awk -v FP=0 -v genlength=3000000000 '{FP = FP + ($3 - $2)};END{print FP/genlength}')

#Complete the TFPN Rates FINAL file
echo -e "${SimuID}\tPLINK\t${INDIV}\t${TP}\t${TN}\t${FP}\t${FN}" >> ./Analyses/WGS_PLINK_TFPN_Rates_GEN1000.txt

}

export -f INDIV_Rates_TRUE_IBD

#Create the header for the file where we'll record TP, FP, TN, FN rates per indiv
echo -e "SimuID\tNB_RAD_FRAG\tINDIV\tTP\tTN\tFP\tFN" > ./Analyses/WGS_PLINK_TFPN_Rates_GEN1000.txt
#Apply function to all BEDFILES
find ./TFPN_Rates_PLINK/TRUE_IBD/1000GEN/bedfiles/ -name "TRUEIBD_BED_*_Ind_*.bed" | xargs -I {} bash -c 'INDIV_Rates_TRUE_IBD {}'



