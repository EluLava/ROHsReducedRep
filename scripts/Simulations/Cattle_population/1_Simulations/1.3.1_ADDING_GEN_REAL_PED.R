ped = read.table("Druet_Rec_pedigree.txt")
colnames(ped) = c("ID", "Dam", "Sire")

#Change 0 to NA for no parents
ped[ped == 0] = NA

#identify founders
founders = ped[is.na(ped$Dam) & is.na(ped$Sire),]

#Add a generation column to the ped
pedgen = cbind(gen=vector(mode = "numeric", length = nrow(ped)), ped)

#Create an empty vector we'll fill with every "adult" individuals
#first, only the founders
adults=founders$ID

#Create empty vector we'll fill with juveniles
juveniles=vector(mode = "character", length = 0)

#create the generation count start at ZERO for founders
gen=0

#Loop through the rows of the gen pedigree
for(row in 1:nrow(pedgen)){
  #If both parents in adults OR NA, we keep same gen
  if((pedgen$Dam[row] %in% adults & pedgen$Sire[row] %in% adults) | (is.na(pedgen$Dam[row]) & is.na(pedgen$Sire[row]))){
    pedgen$gen[row] = gen
    #Add ID to juvenile list
    juveniles = append(as.character(pedgen$ID[row]), juveniles)
  #else if one parent in juveniles, change generation
  } else if(pedgen$Dam[row] %in% juveniles | pedgen$Sire[row] %in% juveniles) {
    #Change generation
    gen = gen + 1
    pedgen$gen[row] = gen
    #Pass all juveniles to adults
    adults = append(juveniles, adults)
    #Reset juveniles to 0 becasue all became adults BUT new indiv = NEW JUVENILE
    juveniles = as.vector(pedgen$ID[row],mode = "character")
    
  } else {
    #Keep same generation
    pedgen$gen[row] = gen
    #Add kid in juveniles
    juveniles = append(as.character(pedgen$ID[row]), juveniles)
  }
}

system("mkdir -p ./Real_pedigree_recmap/")
write.table(pedgen, "./Real_pedigree_recmap/ped_with_gen", quote = F, col.names = T, row.names = F)
