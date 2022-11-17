#List the VCF
VCFs=$(find ./Final_VCF/ -name "*.vcf.gz")

#Create output directory
mkdir -p ./True_IBD/100GEN
mkdir -p ./True_IBD/1000GEN

#Loop through VCFs
for file in ${VCFs};
do

	#Extarct SimID
	SimID=$(basename -s ".vcf.gz" ${file} | cut -d'_' -f2)

	#Loop through Chromosomes
	for CHR in {1..30}; do

		#Python
		python3 - << EOF

import msprime, tskit, pyslim, sys, pandas as pd

SimuID = "${SimID}"
CHR = "${CHR}"

#TMRCAs in GEN HERE 100, just change to 1000 if needed !
MAXibd = 100

myPATH='.'

#Load the corresonding tree
trees = tskit.load('{}/MSprime_Recap_Mut_VCF/Subset_TreeSeq_{}_CHR{}.trees'.format(myPATH, SimuID, CHR))

#Read Subset individuals
#Create an emtpy vector first
lines = []
#Read the individuals file we want to subset
with open('{}/subsampled_individuals/subsampled_Individuals_{}'.format(myPATH, SimuID), "r") as file:
	for line in file:
		lines.append(int(line.strip()))

#Extract SUBSET individuals
subs_inds = []
#Loop through individuals  ALIVE (faster) in the dataset
for i in pyslim.individuals_alive_at(trees,0):
	ind = trees.individual(i)
	if(ind.metadata['pedigree_id'] in lines):
		subs_inds.append(i)

#/!\ IF INDIV NAME = 0, we need to change it to 1000 in the final df

#Create the names for
df_names = ['SimuID', 'ID', 'CHR', 'POS1', 'POS2']

#Create an empty panda.df with SimuID, ID, CHR, POS1, POS2
ibd_segments = pd.DataFrame(columns = df_names)

#Loop through individuals
for i in subs_inds:
	#Extract individual
	ind = trees.individual(i)
	#Extract nodes for this individual
	node1 = ind.nodes[0]
	node2 = ind.nodes[1]
	#IF INDIVIDUAL = 0, change its name to 1000
	if(i == 0):
		ID = 1000
	else:
		ID = subs_inds.index(i)
	#Loop through trees
	for ts in trees.trees():
		#Estimate TMRCA between both nodes
		time = ts.tmrca(node1,node2)
		#IF TMRCA > our treshold, NOT IBD, if < = IBD
		if(time < MAXibd):
			#Extract start ans stop position of the genomic segment
			POS1 = round(ts.get_interval()[0])
			POS2 = round(ts.get_interval()[1])
			#Check if the ibd df already contains on line with this individual (becasue if not ERROR, SUBSCRIPT OUT OF BOND); IF YES, check if last position etc.; if not, just add the line
			if(ID in ibd_segments.ID.values):
				#if the last segment of this individual stops where the new one starts, just REPLACE POS2 by new POS2; else just add a line
				if(ibd_segments.POS2.values[-1] == POS1):
					#JUST REPLACE POS2 (5th position)
					ibd_segments.iloc[-1,ibd_segments.columns.get_loc('POS2')] = POS2
				else:
					#Append the SEGMENT ROW to our final DF
					ibd_segments = ibd_segments.append({'SimuID' : SimuID, 'ID' : ID, 'CHR' : CHR, 'POS1' : POS1, 'POS2' : POS2}, ignore_index = True)
			else:
				#Append the SEGMENT ROW to our final DF
				ibd_segments = ibd_segments.append({'SimuID' : SimuID, 'ID' : ID, 'CHR' : CHR, 'POS1' : POS1, 'POS2' : POS2}, ignore_index = True)

#Write the FINAL IIBD SEGMENTS DF
ibd_segments.to_csv('{}/TRUE_IBD/{}GEN/TRUE_IBD_SEGMENTS_{}_{}.csv'.format(myPATH, MAXibd, SimuID, CHR), index = False, sep = "\t")

EOF

##Merge final DFs CHANGE HERE TO if 1000 GEN !
SimuIDs=$(find ./Final_VCFs/ -name "FINAL_*" | xargs -I {} bash -c 'basename -s ".vcf.gz" {} | cut -d'_' -f2' | sort | uniq)
echo -e "SimuID\tID\tCHR\tPOS1\tPOS2" > ./TRUE_IBD/100GEN/stats/TRUE_IBD_SEGMENTS_GEN100.txt
for sim in $SimuIDs; do FILES=$(find ./TRUE_IBD/100GEN/stats/ -name "TRUE_IBD_SEGMENTS_${sim}_*.csv" | sort --version-sort); for f in $FILES; do awk 'NR > 1 {print $0}' $f >> ./TRUE_IBD/100GEN/stats/TRUE_IBD_SEGMENTS_GEN100.txt; done; done
