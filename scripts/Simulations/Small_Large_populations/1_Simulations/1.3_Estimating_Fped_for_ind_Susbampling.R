#This script runs on my computer

#Install packages if not installed already
my_packages <- c("pedigree", "foreach", "doParallel")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)

#Load packages
library(pedigree)
library(foreach)
library(doParallel)

#Creating output directory if needed
system("mkdir -p ./subsampled_Individuals")

#Listing the mating files in dir
filesmating = list.files("./Simulations", pattern = "mating")
filesdeaths = list.files("./Simulations", pattern = "death")
#Defining inbreeding classes
inbreeding_classes = c(0,0.1,0.2,0.3,0.4,0.5)

#Choosing the nb of cores
registerDoParallel(cores=5)
#Looping through the files
foreach(file=1:length(filesmating)) %dopar% {

  #Extracting simu IDs
  SimID = strsplit(strsplit(filesmating[file], "_", fixed = TRUE)[[1]][3], ".", fixed = TRUE)[[1]][1]

  #Reading mating files
  dta_mating = read.table(paste0("./Simulations/",filesmating[file]))
  colnames(dta_mating) = c("Generation", "Mate1", "Mate2", "Child")

  #Reading deaths files
  dta_death = read.table(paste0("./Simulations/",filesdeaths[file]))
  colnames(dta_death) = c("Generation", "Individuals")

  #Subsetting the individuals alive at last generation
  ind_alive = dta_mating$Child[!dta_mating$Child %in% dta_death$Individuals]

  #Subsetting only last 15 generations
  dta_mating = dta_mating[dta_mating$Generation > 985,]

  #Creating pedigree format df
  ped = as.data.frame(cbind(dta_mating$Child, dta_mating$Mate1, dta_mating$Mate2))
  colnames(ped) = c("Child", "Mate1", "Mate2")

  #Calculating inbreeding coefficient
  inbreeding = calcInbreeding(ped)
  #Adding the individuals tag as names
  names(inbreeding) = ped$Child

  #Extracting only the alive individuals
  inbreeding_alive = inbreeding[names(inbreeding) %in% ind_alive]

  #Creating empty list to store individuals I'll subsample
  sub_ind = vector(length = 0)
  #Loop through calsses of inbreeding to select INDVs
  for (class in 2:length(inbreeding_classes)) {
    #INDVs within the class
    inds = inbreeding_alive[inbreeding_alive > inbreeding_classes[class - 1] & inbreeding_alive < inbreeding_classes[class]]
    #If less or extactly 20 individuals in this class, take them all
    if(length(inds) <= 20){
      sub_ind = append(sub_ind, inds)
    #else randomly subsample
    } else {
      sub_ind = append(sub_ind, sample(inds, 20, replace = FALSE))
    }
  }

  #save sub inds in a object including the simID in its name
  assign(paste0("sub_", SimID), names(sub_ind))

  #saving the list of INDVs
  write.table(names(sub_ind), paste0("./subsampled_Individuals/subsampled_Individuals_", SimID), row.names = FALSE, quote = FALSE, col.names = FALSE)
}