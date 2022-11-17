################################ Python msprime script #############################

import msprime, pyslim, random, sys

#First argument given to the script is the tree file
tree = sys.argv[1]
#Second argument is the random SimID
simID = sys.argv[2]
#Third argument is CHR
CHR = sys.argv[3]

#Defining Ne
Ne = 1000

#### Load the tree ####
tree_slim = pyslim.load(tree) #No simplify for recapitation

#### Recapitate ####
recap = tree_slim.recapitate(recombination_rate = 1e-8, Ne = Ne).simplify()

#### Overlay mutations ####
ts = pyslim.SlimTreeSequence(msprime.mutate(recap, rate = 2.5e-8))

ts.dump("./TreeSequences/Full_TreeSeq_{}_{}.trees".format(simID,CHR))

#### Subset individuals  ####
#Create an emtpy vector first
lines = []
#Read the individuals file we want to subset
with open("./subsampled_Individuals/subsampled_Individuals_{}".format(simID), "r") as file:
        for line in file:
                lines.append(int(line.strip()))

subs_inds = []
for i in ts.individuals_alive_at(0):
        ind = ts.individual(i)
        if(ind.metadata["pedigree_id"] in lines):
                subs_inds.append(i)

with open("./VCFs/{}_Ne{}_{}.vcf".format(simID, int(Ne), CHR), "w") as simuvcf:
        ts.write_vcf(simuvcf, individuals = subs_inds, position_transform = "legacy")

#To have the reduced tree

#Create an empty list which will contain nodes associated with subsampled individuals
sampled_nodes = []
#Looping through individuals to fill the nodes list
for ind in subs_inds:
        sampled_nodes.extend(ts.individual(ind).nodes)

#Simplify to remove the rest of the tree sequence
tree_simp = ts.simplify(sampled_nodes)

tree_simp.dump("./TreeSequences/Subset_TreeSeq_{}_{}.trees".format(simID, CHR))
