import msprime, pyslim, random, sys, os

mypath='.'

simID = int(sys.argv[1])

CHR = int(sys.argv[2])

LEN = int(sys.argv[3])

#Effective population size WILD
Ne_Wild = 50000
#Effective population size
Ne_Domestic = 1500

#Time Wild and Domestic populations split (forward) or merge (backward) IN GENERATIONS
T1 = 10000

#Migration rate between Wild and domestic Populations
M1 = 3e-5

#The number of individuals which will be sampled --> The creation of the breed bottleneck
N_FIN=200

#Populations Configurations --> Ne and sample size in the end
population_configurations = [msprime.PopulationConfiguration(sample_size = 0, initial_size = Ne_Wild),
								msprime.PopulationConfiguration(sample_size = N_FIN, initial_size = Ne_Domestic)]

#Create the migration matrix
migration_matrix = [[0, M1], [M1, 0]]

#DEMOGRAPHY
demographic_events = [msprime.MassMigration(time = T1, source = 1, destination = 0, proportion = 1.0), msprime.MigrationRateChange(time = T1, rate = 0)]

#SIMULATING THE TREE
ts = msprime.simulate(population_configurations = population_configurations, migration_matrix = migration_matrix, demographic_events = demographic_events,
	recombination_rate = 1e-8, length = LEN, random_seed = simID)

#ANNOTATING THE TREE with PYSLIM
tables = ts.dump_tables()
pyslim.annotate_defaults_tables(tables, model_type = "nonWF", slim_generation=1)

#Read the sex table

#INDIVIDUAL METADATA
individual_metadata=list(pyslim.extract_individual_metadata(tables))
for j in range(len(individual_metadata)):
	#Assign F or M sex 1 indiv out of 2
	individual_metadata[j].sex = [pyslim.INDIVIDUAL_TYPE_MALE, pyslim.INDIVIDUAL_TYPE_FEMALE][j%2]
	individual_metadata[j].age = random.choice([0,1,2,3])

pyslim.annotate_individual_metadata(tables, individual_metadata)
slim_tree = pyslim.load_tables(tables)

#CLUSTER
#slim_tree.dump("/scratch/wally/FAC/FBM/DEE/jgoudet/default/elavanch/Which_data/Data/Cattle_Simulations/BURN_IN/BI_Full_TreeSeq_Cattle_{}_CHR_{}".format(simID, CHR))

#MyComupter
slim_tree.dump(os.path.join(mypath, 'BURN_IN/BI_Full_TreeSeq_Cattle_{}_CHR_{}.trees'.format(simID, CHR)))

