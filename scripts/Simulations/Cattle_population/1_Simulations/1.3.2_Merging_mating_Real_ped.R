library(pedigree)

#Read ped table
ped = read.table("./Real_pedigree_recmap/ped_with_gen", header = T)

#Change gen number as gen+125 --> Last = gen 200
ped$gen = ped$gen + 125

#Change columns order
pedordered = as.data.frame(cbind(as.integer(ped$gen), as.character(ped$Dam), as.character(ped$Sire), as.character(ped$ID)), stringsAsFactors = F)
colnames(pedordered) = c("Gen", "Dam", "Sire", "ID")
pedordered$Gen = as.integer(pedordered$Gen)

#Read mating table
mating = read.table("./Slim_creating_inbreeding/mating_Sim17390.txt")
#Add colnames
colnames(mating) = c("Gen", "Dam", "Sire", "ID")
mating$Dam = as.character(mating$Dam)
mating$Sire = as.character(mating$Sire)
mating$ID = as.character(mating$ID)

#Read death table
deaths = read.table("./Slim_creating_inbreeding/death_Sim17390.txt")
#Add colnames
colnames(deaths) = c("Gen", "ID")

#Merge both PED + MATING
dtafull = rbind(pedordered, mating)
#Order by generation
dtafull = dtafull[order(dtafull$Gen),]

#Loop through generations to assign random individuals to NA from real ped
for (gen in 125:200) {
  #Subsample the generation pedigree
  dta_sub = dtafull[dtafull$Gen == gen,]
  
  #Loop through NA DAM
  for (indiv in dta_sub$ID[is.na(dta_sub$Dam)]) {
    #Select the available DAM from mating file (not real ped)
    dams = mating$Dam[mating$Gen == gen]
    #Randomly assign a DAM
    dtafull$Dam[dtafull$ID == indiv] = sample(dams, 1, replace = F)
  }
  
  #Loop through NA SIRE
  for (indiv in dta_sub$ID[is.na(dta_sub$Sire)]) {
    #Select the available SIRES from mating file (not real ped)
    sires = mating$Sire[mating$Gen == gen]
    #Randomly assign a SIRE
    dtafull$Sire[dtafull$ID == indiv] = sample(sires, 1, replace = F)
  }
}

#Now all indiv have parents --> We are only interested in the indiv from last gen in real ped
#We can trim the pedigree and remove all indiv not ancestors to these ones

#Extract individuals from real ped and last gen
interesting_INDIVs = pedordered$ID[pedordered$Gen == 200]

#LOGICAL vector for all indiv
TFdata = dtafull$ID %in% interesting_INDIVs

#Extract individuals
indiv_trim = trimPed(ped = as.data.frame(cbind(dtafull$ID, dtafull$Dam, dtafull$Sire)), data = TFdata)

#Trim the pedigree
ped_trim = dtafull[indiv_trim,]

#Now we need to pass the weird GR-NAMES to INTEGERS

#extract weird GR names
weird_names = ped_trim$ID[ped_trim$ID %in% unique(ped$ID)]
#Strsplit GRR and numbers
str_weird_names = unlist(strsplit(weird_names, "-"))
#rm GRR & add 90000 to all indiv to ensure they are not in the SLIM ID tag
str_weird_names = as.numeric(as.character(str_weird_names[str_weird_names != "GRR"])) + 90000
#Create a df with weird name as col1 and new int tag as col2 (starts at 90000 to be safe)
name_change = as.data.frame(cbind(str_weird_names, weird_names))
name_change$str_weird_names = as.character(name_change$str_weird_names)
name_change$weird_names = as.character(name_change$weird_names)

#change ped_trim to chr
ped_trim$Dam = as.character(ped_trim$Dam)
ped_trim$Sire = as.character(ped_trim$Sire)
ped_trim$ID = as.character(ped_trim$ID)

#Loop through pedigree and change weird names to new int names
for(indiv in 1:nrow(name_change)){
  ped_trim$ID[ped_trim$ID == (name_change[indiv,2])] = name_change[indiv,1]
  ped_trim$Dam[ped_trim$Dam == (name_change[indiv,2])] = name_change[indiv,1]
  ped_trim$Sire[ped_trim$Sire == (name_change[indiv,2])] = name_change[indiv,1]
}

#Now we need to re-create the matings & death files with new individiduals
#ped_trim is the UPDATED mating.txt file

#Now we need to add one final column: SEX
ped_trim = cbind(ped_trim, Sex=vector(mode = "character", length = nrow(ped_trim)))
ped_trim$Sex = factor(ped_trim$Sex, levels = c(0, 1))
#We loop through individuals in the ped and see if they are in DAM or Sire to add "F" or "M" in sex column
for(ind in unique(ped_trim$ID)){
  #If ind in both ERROR SANITY CHECK
  if(ind %in% unique(ped_trim$Dam) & ind %in% unique(ped_trim$Sire)){
    print(paste0("ERROR: Individual ", ind," in both Dam & Sire"))
  } else if(ind %in% unique(ped_trim$Dam)){
    ped_trim$Sex[ped_trim$ID == ind] = 0
  } else if(ind %in% unique(ped_trim$Sire)){
    ped_trim$Sex[ped_trim$ID == ind] = 1
  }
}

#We add random sex for individuals from GEN 200
for(ind in unique(ped_trim$ID[ped_trim$Gen == 200])){
  ped_trim$Sex[ped_trim$ID == ind] = sample(c(0,1), 1)
}

#RM DAM & SIRE AS FACTOR
ped_trim$Dam = as.character(ped_trim$Dam)
ped_trim$Dam = unlist(ped_trim$Dam)
ped_trim$Sire = as.character(ped_trim$Sire)
ped_trim$Sire = unlist(ped_trim$Sire)

#Write UPDATED MATING FILE
write.table(ped_trim, "./Slim_creating_inbreeding/mating_complete.txt", quote = F, col.names = F, row.names = F)

#Now to create the deaths file we need to record for each ind of the pedigree what is the last gen is appears in

#Create empty df
deathsfull = as.data.frame(matrix(nrow = 0, ncol = 2))

#Loop through all indivdividuals in the trim ped
for (ind in unique(ped_trim$ID)) {
  lastg = max(ped_trim$Gen[ped_trim$Dam == ind | ped_trim$Sire == ind])
  line = cbind(lastg, ind)
  deathsfull = rbind(deathsfull, line)
}

#rm the -Inf values in GEN
deathsfull = deathsfull[deathsfull$lastg != -Inf,]
deathsfull$lastg = as.numeric(as.character(deathsfull$lastg))
deathsfull$ind = as.character(deathsfull$ind)

#Sort the death df by generation
deathsfullordered = deathsfull[order(deathsfull$lastg),]

write.table(deathsfullordered, "./Slim_creating_inbreeding/deaths_complete.txt", quote = F, col.names = F, row.names = F)




