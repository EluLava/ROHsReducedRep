################################ Python msprime script #############################

import msprime, pyslim, random, sys

#First argument given to the script is the tree file
tree = sys.argv[1]
#Second argument is the random SimID
simID = sys.argv[2]
#Third argument is CHR
CHR = sys.argv[3]

#### Load the tree ####
tree_slim = pyslim.load(tree).simplify() #with simplification

#### Overlay mutations ####
ts = pyslim.SlimTreeSequence(msprime.mutate(tree_slim, rate = 2.5e-8))

ts.dump("./MSprime_Recap_Mut_VCF/Full_TreeSeq_{}_{}.trees".format(simID,CHR))

#### Subset individuals ####

#Extract metadata
tables = ts.tables

subs_inds = []
for i in ts.individuals_alive_at(0):
	ind = ts.individual(i)
	md = ind.metadata
	if(md["age"] == 0):
		subs_inds.append(i)

with open("./MSprime_Recap_Mut_VCF/{}_CHR_{}.vcf".format(simID, CHR), "w") as simuvcf:
	ts.write_vcf(simuvcf, individuals = subs_inds, position_transform = "legacy")

#To have the reduced tree

#Create an empty list which will contain nodes associated with subsampled individuals
sampled_nodes = []
#Looping through individuals to fill the nodes list
for ind in subs_inds:
	sampled_nodes.extend(ts.individual(ind).nodes)

#Simplify to remove the rest of the tree sequence
tree_simp = ts.simplify(sampled_nodes)

tree_simp.dump("./MSprime_Recap_Mut_VCF/Subset_TreeSeq_{}_{}.trees".format(simID, CHR))