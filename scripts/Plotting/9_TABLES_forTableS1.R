###############################################################################################
#### TABLE 1: RAD-SEQUENCING PLINK100KB MEAN Rsquared, slope & Intercept among simulations ####
###############################################################################################

#Read the data FROH ne 1'000
dtaFROH1000 = read.table("./SmallPop/Analyses/RADSEQ_PLINK_FROH.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH1000 = dtaFROH1000[dtaFROH1000$WINDOWS.NB %in% c(180000, 240000, 600000),]

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_1000 = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(dtaTable1_1000) = c("PopulationSize", "SimuID", "PercGenSeq", "Replicate", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH1000$SimID)){
  #Loop trhough RAD percentages
  for(rad in sort(unique(dtaFROH1000$WINDOWS.NB))){
    #Loop through replicates
    for(repl in 1:100){
      #Subsample the data
      dta_sub = dtaFROH1000[dtaFROH1000$SimID == sim & dtaFROH1000$WINDOWS.NB == rad & dtaFROH1000$REP == repl,]
      #Estimate Rsquared
      Rsquared = cor(dta_sub$FROH_sub, dta_sub$Froh_TRUE_IBD_100GEN)
      #Build linear model
      mod = lm(dta_sub$FROH_sub ~ dta_sub$Froh_TRUE_IBD_100GEN)
      #Extract Intercept
      inter = mod$coefficients[1]
      #Extract slope
      slo = mod$coefficients[2]
      #Fuse every info in one row table
      dtatablesub = data.frame("1'000", sim, round(((rad)*500)/30000000, digits = 2), repl, Rsquared, inter, slo)
      names(dtatablesub) = colnames(dtaTable1_1000)
      
      #Fill the data
      dtaTable1_1000 = rbind(dtaTable1_1000, dtatablesub)
    }
  }
}

#rm row-names
row.names(dtaTable1_1000) = NULL
#Change 4.92 to 5 (hum)
dtaTable1_1000$PercGenSeq[dtaTable1_1000$PercGenSeq == 4.92] = 5.00
#This is the full table small NE

#Read the data FROH ne 10'000
dtaFROH10000 = read.table("./LargePop/Analyses/RADSEQ_PLINK_FROH.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH10000 = dtaFROH10000[dtaFROH10000$WINDOWS.NB %in% c(15000, 25000, 60000),]

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_10000 = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(dtaTable1_10000) = c("PopulationSize", "SimuID", "PercGenSeq", "Replicate", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH10000$SimID)){
  #Loop trhough RAD percentages
  for(rad in sort(unique(dtaFROH10000$WINDOWS.NB))){
    #Loop through replicates
    for(repl in 1:100){
      #Subsample the data
      dta_sub = dtaFROH10000[dtaFROH10000$SimID == sim & dtaFROH10000$WINDOWS.NB == rad & dtaFROH10000$REP == repl,]
      #Estimate Rsquared
      Rsquared = cor(dta_sub$FROH_sub, dta_sub$Froh_TRUE_IBD_100GEN)
      #Build linear model
      mod = lm(dta_sub$FROH_sub ~ dta_sub$Froh_TRUE_IBD_100GEN)
      #Extract Intercept
      inter = mod$coefficients[1]
      #Extract slope
      slo = mod$coefficients[2]
      #Fuse every info in one row table
      dtatablesub = data.frame("10'000", sim, round(((rad)*500)/30000000, digits = 2), repl, Rsquared, inter, slo)
      names(dtatablesub) = colnames(dtaTable1_10000)
      
      #Fill the data
      dtaTable1_10000 = rbind(dtaTable1_10000, dtatablesub)
    }
  }
}

#rm row-names
row.names(dtaTable1_10000) = NULL
#This is the full table large NE

#system("mkdir -p ./TABLES")

#Merge both small and large POP for full
dtaFullS1 = rbind(dtaTable1_1000, dtaTable1_10000)
#Write this as SUPPLEMENTARY TABLE
#write.csv(dtaFullS1, "./TABLES/FROH_RAD_SEQUENCING_PLINK100KB_bothPops_SUPPLEMENTARY_fullTABLE.csv",
          row.names = F, quote = F)

#Create the TABLE 1 (mean Rsquared, slope & Intercept per PercSeqGen)
dtaTable1MEAN = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$PercGenSeq), FUN = mean)
dtaTable1MEAN = dtaTable1MEAN[order(dtaTable1MEAN$Group.1, dtaTable1MEAN$Group.2),]
dtaTable1SD = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$PercGenSeq), FUN = sd)
dtaTable1SD = dtaTable1SD[order(dtaTable1SD$Group.1, dtaTable1SD$Group.2),]

#Create emtpy dataframe we'll fill with mean ± sd
dtaTable1 = data.frame(matrix(nrow = 0, ncol = 5))
colnames(dtaTable1) = c("PopulationSize", "PercGenSeq", "Rsquared", "Intercept", "Slope")

#Loop through rows of mean & sd df
for (row in 1:nrow(dtaTable1MEAN)) {
  #Extract PopSize
  Popsize = dtaTable1MEAN$Group.1[row]
  #Extract PerSeq
  PercSeq = dtaTable1MEAN$Group.2[row]
  #Extract mean RSquared ± sd
  Rsquared = paste(round(dtaTable1MEAN$Rsquared[row], digits = 3), round(dtaTable1SD$Rsquared[row], digits = 3), sep = " ± ")
  #Extract mean Intercept ± sd
  Intercept = paste(round(dtaTable1MEAN$Intercept[row], digits = 3), round(dtaTable1SD$Intercept[row], digits = 3), sep = " ± ")
  #Extract mean Slope ± sd
  Slope = paste(round(dtaTable1MEAN$Slope[row], digits = 3), round(dtaTable1SD$Slope[row], digits = 3), sep = " ± ")
  #Fuse in one row
  rowsub = data.frame(Popsize, PercSeq, Rsquared, Intercept, Slope)
  #Set colnames for rbind
  colnames(rowsub) = colnames(dtaTable1)
  #Rbind
  dtaTable1 = rbind(dtaTable1, rowsub)
}


#This is the complete TABLE 1
write.csv(dtaTable1, "./TABLES/FROH_RAD_SEQUENCING_PLINK100KB_bothPops_TABLE1.csv",
          row.names = F, quote = F)

######################################################################################################
#### TABLE 2: RAD-SEQUENCING PLINK SEQ FHOM MEAN Rsquared, slope & Intercept among simulations #######
######################################################################################################

#Read the data FROH ne 1'000
dtaFROH1000 = read.table("./SmallPop/Analyses/RADSEQ_PLINK_FHOM.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH1000 = dtaFROH1000[dtaFROH1000$WINDOWS.NB %in% c(180000, 240000, 600000),]

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_1000 = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(dtaTable1_1000) = c("PopulationSize", "SimuID", "PercGenSeq", "Replicate", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH1000$SimID)){
  #Loop trhough RAD percentages
  for(rad in sort(unique(dtaFROH1000$WINDOWS.NB))){
    #Loop through replicates
      #Subsample the data
      dta_sub = dtaFROH1000[dtaFROH1000$SimID == sim & dtaFROH1000$WINDOWS.NB == rad,]
      #Estimate Rsquared
      Rsquared = cor(dta_sub$FHOM_sub, dta_sub$Froh_TRUE_IBD_100GEN)
      #Build linear model
      mod = lm(dta_sub$FHOM_sub ~ dta_sub$Froh_TRUE_IBD_100GEN)
      #Extract Intercept
      inter = mod$coefficients[1]
      #Extract slope
      slo = mod$coefficients[2]
      #Fuse every info in one row table
      dtatablesub = data.frame("1'000", sim, round(((rad)*500)/30000000, digits = 2), repl, Rsquared, inter, slo)
      names(dtatablesub) = colnames(dtaTable1_1000)
      
      #Fill the data
      dtaTable1_1000 = rbind(dtaTable1_1000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_1000) = NULL
#Change 4.92 to 5 (hum)
dtaTable1_1000$PercGenSeq[dtaTable1_1000$PercGenSeq == 4.92] = 5.00
#This is the full table small NE

#Read the data FROH ne 10'000
dtaFROH10000 = read.table("./LargePop/Analyses/RADSEQ_PLINK_FHOM.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH10000 = dtaFROH10000[dtaFROH10000$WINDOWS.NB %in% c(15000, 25000, 60000),]

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_10000 = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(dtaTable1_10000) = c("PopulationSize", "SimuID", "PercGenSeq", "Replicate", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH10000$SimID)){
  #Loop trhough RAD percentages
  for(rad in sort(unique(dtaFROH10000$WINDOWS.NB))){
      #Subsample the data
      dta_sub = dtaFROH10000[dtaFROH10000$SimID == sim & dtaFROH10000$WINDOWS.NB == rad,]
      #Estimate Rsquared
      Rsquared = cor(dta_sub$FHOM_sub, dta_sub$Froh_TRUE_IBD_100GEN)
      #Build linear model
      mod = lm(dta_sub$FHOM_sub ~ dta_sub$Froh_TRUE_IBD_100GEN)
      #Extract Intercept
      inter = mod$coefficients[1]
      #Extract slope
      slo = mod$coefficients[2]
      #Fuse every info in one row table
      dtatablesub = data.frame("10'000", sim, round(((rad)*500)/30000000, digits = 2), repl, Rsquared, inter, slo)
      names(dtatablesub) = colnames(dtaTable1_10000)
      
      #Fill the data
      dtaTable1_10000 = rbind(dtaTable1_10000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_10000) = NULL
#This is the full table large NE

#system("mkdir -p ./TABLES")

#Merge both small and large POP for full
dtaFullS1 = rbind(dtaTable1_1000, dtaTable1_10000)
#Write this as SUPPLEMENTARY TABLE
#write.csv(dtaFullS1, "./TABLES/FROH_RAD_SEQUENCING_PLINK100KB_bothPops_SUPPLEMENTARY_fullTABLE.csv",
#row.names = F, quote = F)

#Create the TABLE 1 (mean Rsquared, slope & Intercept per PercSeqGen)
dtaTable1MEAN = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$PercGenSeq), FUN = mean)
dtaTable1MEAN = dtaTable1MEAN[order(dtaTable1MEAN$Group.1, dtaTable1MEAN$Group.2),]
dtaTable1SD = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$PercGenSeq), FUN = sd)
dtaTable1SD = dtaTable1SD[order(dtaTable1SD$Group.1, dtaTable1SD$Group.2),]

#Create emtpy dataframe we'll fill with mean ± sd
dtaTable1 = data.frame(matrix(nrow = 0, ncol = 5))
colnames(dtaTable1) = c("PopulationSize", "PercGenSeq", "Rsquared", "Intercept", "Slope")

#Loop through rows of mean & sd df
for (row in 1:nrow(dtaTable1MEAN)) {
  #Extract PopSize
  Popsize = dtaTable1MEAN$Group.1[row]
  #Extract PerSeq
  PercSeq = dtaTable1MEAN$Group.2[row]
  #Extract mean RSquared ± sd
  Rsquared = paste(round(dtaTable1MEAN$Rsquared[row], digits = 3), round(dtaTable1SD$Rsquared[row], digits = 3), sep = " ± ")
  #Extract mean Intercept ± sd
  Intercept = paste(round(dtaTable1MEAN$Intercept[row], digits = 3), round(dtaTable1SD$Intercept[row], digits = 3), sep = " ± ")
  #Extract mean Slope ± sd
  Slope = paste(round(dtaTable1MEAN$Slope[row], digits = 3), round(dtaTable1SD$Slope[row], digits = 3), sep = " ± ")
  #Fuse in one row
  rowsub = data.frame(Popsize, PercSeq, Rsquared, Intercept, Slope)
  #Set colnames for rbind
  colnames(rowsub) = colnames(dtaTable1)
  #Rbind
  dtaTable1 = rbind(dtaTable1, rowsub)
}


#This is the complete TABLE 1
write.csv(dtaTable1, "./TABLES/FHOM_RAD_SEQUENCING_PLINK100KB_bothPops_TABLE1.csv",
          row.names = F, quote = F)

###############################################################################################
#### TABLE 3: RAD-SEQUENCING RZooROH MEAN Rsquared, slope & Intercept among simulations #######
###############################################################################################

#Read the data FROH ne 1'000
dtaFROH1000 = read.table("./SmallPop/Analyses/RADSEQ_RZooRoH_FROH.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH1000 = dtaFROH1000[dtaFROH1000$WINDOWS.NB %in% c(3000, 9000, 60000),]

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_1000 = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(dtaTable1_1000) = c("PopulationSize", "SimuID", "PercGenSeq", "Replicate", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH1000$SimID)){
  #Loop trhough RAD percentages
  for(rad in sort(unique(dtaFROH1000$WINDOWS.NB))){
    #Loop through replicates
    for(repl in 1:100){
      #Subsample the data
      dta_sub = dtaFROH1000[dtaFROH1000$SimID == sim & dtaFROH1000$WINDOWS.NB == rad & dtaFROH1000$REP == repl,]
      #Estimate Rsquared
      Rsquared = cor(dta_sub$FROH_sub, dta_sub$Froh_TRUE_IBD_100GEN)
      #Build linear model
      mod = lm(dta_sub$FROH_sub ~ dta_sub$Froh_TRUE_IBD_100GEN)
      #Extract Intercept
      inter = mod$coefficients[1]
      #Extract slope
      slo = mod$coefficients[2]
      #Fuse every info in one row table
      dtatablesub = data.frame("1'000", sim, round(((rad)*500)/30000000, digits = 2), repl, Rsquared, inter, slo)
      names(dtatablesub) = colnames(dtaTable1_1000)
      
      #Fill the data
      dtaTable1_1000 = rbind(dtaTable1_1000, dtatablesub)
    }
  }
}

#rm row-names
row.names(dtaTable1_1000) = NULL
#Change 4.92 to 5 (hum)
dtaTable1_1000$PercGenSeq[dtaTable1_1000$PercGenSeq == 4.92] = 5.00
#This is the full table small NE

#Read the data FROH ne 10'000
dtaFROH10000 = read.table("./LargePop/Analyses/RADSEQ_RZooRoH_FROH.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH10000 = dtaFROH10000[dtaFROH10000$WINDOWS.NB %in% c(250, 500, 7500),]

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_10000 = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(dtaTable1_10000) = c("PopulationSize", "SimuID", "PercGenSeq", "Replicate", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH10000$SimID)){
  #Loop trhough RAD percentages
  for(rad in sort(unique(dtaFROH10000$WINDOWS.NB))){
    #Loop through replicates
    for(repl in 1:100){
      #Subsample the data
      dta_sub = dtaFROH10000[dtaFROH10000$SimID == sim & dtaFROH10000$WINDOWS.NB == rad & dtaFROH10000$REP == repl,]
      #Estimate Rsquared
      Rsquared = cor(dta_sub$FROH_sub, dta_sub$Froh_TRUE_IBD_100GEN)
      #Build linear model
      mod = lm(dta_sub$FROH_sub ~ dta_sub$Froh_TRUE_IBD_100GEN)
      #Extract Intercept
      inter = mod$coefficients[1]
      #Extract slope
      slo = mod$coefficients[2]
      #Fuse every info in one row table
      dtatablesub = data.frame("10'000", sim, round(((rad)*500)/30000000, digits = 2), repl, Rsquared, inter, slo)
      names(dtatablesub) = colnames(dtaTable1_10000)
      
      #Fill the data
      dtaTable1_10000 = rbind(dtaTable1_10000, dtatablesub)
    }
  }
}

#rm row-names
row.names(dtaTable1_10000) = NULL
#This is the full table large NE

#system("mkdir -p ./TABLES")

#Merge both small and large POP for full
dtaFullS1 = rbind(dtaTable1_1000, dtaTable1_10000)
#Write this as SUPPLEMENTARY TABLE
#write.csv(dtaFullS1, "./TABLES/FROH_RAD_SEQUENCING_PLINK100KB_bothPops_SUPPLEMENTARY_fullTABLE.csv",
#row.names = F, quote = F)

#Create the TABLE 1 (mean Rsquared, slope & Intercept per PercSeqGen)
dtaTable1MEAN = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$PercGenSeq), FUN = mean)
dtaTable1MEAN = dtaTable1MEAN[order(dtaTable1MEAN$Group.1, dtaTable1MEAN$Group.2),]
dtaTable1SD = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$PercGenSeq), FUN = sd)
dtaTable1SD = dtaTable1SD[order(dtaTable1SD$Group.1, dtaTable1SD$Group.2),]

#Create emtpy dataframe we'll fill with mean ± sd
dtaTable1 = data.frame(matrix(nrow = 0, ncol = 5))
colnames(dtaTable1) = c("PopulationSize", "PercGenSeq", "Rsquared", "Intercept", "Slope")

#Loop through rows of mean & sd df
for (row in 1:nrow(dtaTable1MEAN)) {
  #Extract PopSize
  Popsize = dtaTable1MEAN$Group.1[row]
  #Extract PerSeq
  PercSeq = dtaTable1MEAN$Group.2[row]
  #Extract mean RSquared ± sd
  Rsquared = paste(round(dtaTable1MEAN$Rsquared[row], digits = 3), round(dtaTable1SD$Rsquared[row], digits = 3), sep = " ± ")
  #Extract mean Intercept ± sd
  Intercept = paste(round(dtaTable1MEAN$Intercept[row], digits = 3), round(dtaTable1SD$Intercept[row], digits = 3), sep = " ± ")
  #Extract mean Slope ± sd
  Slope = paste(round(dtaTable1MEAN$Slope[row], digits = 3), round(dtaTable1SD$Slope[row], digits = 3), sep = " ± ")
  #Fuse in one row
  rowsub = data.frame(Popsize, PercSeq, Rsquared, Intercept, Slope)
  #Set colnames for rbind
  colnames(rowsub) = colnames(dtaTable1)
  #Rbind
  dtaTable1 = rbind(dtaTable1, rowsub)
}


#This is the complete TABLE 1
write.csv(dtaTable1, "./TABLES/FROH_RAD_SEQUENCING_RZooROH_bothPops_TABLE1.csv",
          row.names = F, quote = F)

######################################################################################################
#### TABLE 4: RAD-SEQUENCING RZooROH SEQ FHOM MEAN Rsquared, slope & Intercept among simulations #######
######################################################################################################

#Read the data FROH ne 1'000
dtaFROH1000 = read.table("./SmallPop/Analyses/RADSEQ_RZooRoH_FHOM.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH1000 = dtaFROH1000[dtaFROH1000$WINDOWS.NB %in% c(3000, 9000, 60000),]

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_1000 = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(dtaTable1_1000) = c("PopulationSize", "SimuID", "PercGenSeq", "Replicate", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH1000$SimID)){
  #Loop trhough RAD percentages
  for(rad in sort(unique(dtaFROH1000$WINDOWS.NB))){
    #Loop through replicates
    #Subsample the data
    dta_sub = dtaFROH1000[dtaFROH1000$SimID == sim & dtaFROH1000$WINDOWS.NB == rad,]
    #Estimate Rsquared
    Rsquared = cor(dta_sub$FHOM_sub, dta_sub$Froh_TRUE_IBD_100GEN)
    #Build linear model
    mod = lm(dta_sub$FHOM_sub ~ dta_sub$Froh_TRUE_IBD_100GEN)
    #Extract Intercept
    inter = mod$coefficients[1]
    #Extract slope
    slo = mod$coefficients[2]
    #Fuse every info in one row table
    dtatablesub = data.frame("1'000", sim, round(((rad)*500)/30000000, digits = 2), repl, Rsquared, inter, slo)
    names(dtatablesub) = colnames(dtaTable1_1000)
    
    #Fill the data
    dtaTable1_1000 = rbind(dtaTable1_1000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_1000) = NULL
#Change 4.92 to 5 (hum)
dtaTable1_1000$PercGenSeq[dtaTable1_1000$PercGenSeq == 4.92] = 5.00
#This is the full table small NE

#Read the data FROH ne 10'000
dtaFROH10000 = read.table("./LargePop/Analyses/RADSEQ_RZooRoH_FHOM.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH10000 = dtaFROH10000[dtaFROH10000$WINDOWS.NB %in% c(250, 500, 7500),]

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_10000 = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(dtaTable1_10000) = c("PopulationSize", "SimuID", "PercGenSeq", "Replicate", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH10000$SimID)){
  #Loop trhough RAD percentages
  for(rad in sort(unique(dtaFROH10000$WINDOWS.NB))){
    #Subsample the data
    dta_sub = dtaFROH10000[dtaFROH10000$SimID == sim & dtaFROH10000$WINDOWS.NB == rad,]
    #Estimate Rsquared
    Rsquared = cor(dta_sub$FHOM_sub, dta_sub$Froh_TRUE_IBD_100GEN)
    #Build linear model
    mod = lm(dta_sub$FHOM_sub ~ dta_sub$Froh_TRUE_IBD_100GEN)
    #Extract Intercept
    inter = mod$coefficients[1]
    #Extract slope
    slo = mod$coefficients[2]
    #Fuse every info in one row table
    dtatablesub = data.frame("10'000", sim, round(((rad)*500)/30000000, digits = 2), repl, Rsquared, inter, slo)
    names(dtatablesub) = colnames(dtaTable1_10000)
    
    #Fill the data
    dtaTable1_10000 = rbind(dtaTable1_10000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_10000) = NULL
#This is the full table large NE

#system("mkdir -p ./TABLES")

#Merge both small and large POP for full
dtaFullS1 = rbind(dtaTable1_1000, dtaTable1_10000)
#Write this as SUPPLEMENTARY TABLE
#write.csv(dtaFullS1, "./TABLES/FROH_RAD_SEQUENCING_PLINK100KB_bothPops_SUPPLEMENTARY_fullTABLE.csv",
#row.names = F, quote = F)

#Create the TABLE 1 (mean Rsquared, slope & Intercept per PercSeqGen)
dtaTable1MEAN = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$PercGenSeq), FUN = mean)
dtaTable1MEAN = dtaTable1MEAN[order(dtaTable1MEAN$Group.1, dtaTable1MEAN$Group.2),]
dtaTable1SD = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$PercGenSeq), FUN = sd)
dtaTable1SD = dtaTable1SD[order(dtaTable1SD$Group.1, dtaTable1SD$Group.2),]

#Create emtpy dataframe we'll fill with mean ± sd
dtaTable1 = data.frame(matrix(nrow = 0, ncol = 5))
colnames(dtaTable1) = c("PopulationSize", "PercGenSeq", "Rsquared", "Intercept", "Slope")

#Loop through rows of mean & sd df
for (row in 1:nrow(dtaTable1MEAN)) {
  #Extract PopSize
  Popsize = dtaTable1MEAN$Group.1[row]
  #Extract PerSeq
  PercSeq = dtaTable1MEAN$Group.2[row]
  #Extract mean RSquared ± sd
  Rsquared = paste(round(dtaTable1MEAN$Rsquared[row], digits = 3), round(dtaTable1SD$Rsquared[row], digits = 3), sep = " ± ")
  #Extract mean Intercept ± sd
  Intercept = paste(round(dtaTable1MEAN$Intercept[row], digits = 3), round(dtaTable1SD$Intercept[row], digits = 3), sep = " ± ")
  #Extract mean Slope ± sd
  Slope = paste(round(dtaTable1MEAN$Slope[row], digits = 3), round(dtaTable1SD$Slope[row], digits = 3), sep = " ± ")
  #Fuse in one row
  rowsub = data.frame(Popsize, PercSeq, Rsquared, Intercept, Slope)
  #Set colnames for rbind
  colnames(rowsub) = colnames(dtaTable1)
  #Rbind
  dtaTable1 = rbind(dtaTable1, rowsub)
}


#This is the complete TABLE 1
write.csv(dtaTable1, "./TABLES/FHOM_RAD_SEQUENCING_RZooROH_bothPops_TABLE1.csv",
          row.names = F, quote = F)

###########################################################################################
#### TABLE 5: SNP-ARRAYs PLINK100KB MEAN Rsquared, slope & Intercept among simulations ####
###########################################################################################

#Read the data Ne 1'000 50K
dtaFROH100050k = read.table("./SmallPop/Analyses/SmallArray_PLINK_FROH.txt", header = T)[,c(1:3,8)]
dtaFROH100050k = cbind(dtaFROH100050k[,1], SUB_meth=rep("small", nrow(dtaFROH100050k)), dtaFROH100050k[,2:4])
colnames(dtaFROH100050k) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")
#Read the data Ne 1'000 700K
dtaFROH1000700k = read.table("./SmallPop/Analyses/LargeArray_PLINK_FROH.txt", header = T)[,c(1:3,8)]
dtaFROH1000700k = cbind(dtaFROH1000700k[,1], SUB_meth=rep("large", nrow(dtaFROH1000700k)), dtaFROH1000700k[,2:4])
colnames(dtaFROH1000700k) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")
#Merge both
dtaFROH1000 = rbind(dtaFROH100050k, dtaFROH1000700k)

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_1000 = as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(dtaTable1_1000) = c("PopulationSize", "SimuID", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH1000$Simu_ID)){
  #Loop trhough array sizes
  for(rad in sort(unique(dtaFROH1000$SUB_meth))){
    #Subsample the data
    dta_sub = dtaFROH1000[dtaFROH1000$Simu_ID == sim & dtaFROH1000$SUB_meth == rad,]
    #Estimate Rsquared
    Rsquared = cor(dta_sub$Froh_array, dta_sub$Froh_TRUE_IBD_100GEN)
    #Build linear model
    mod = lm(dta_sub$Froh_array ~ dta_sub$Froh_TRUE_IBD_100GEN)
    #Extract Intercept
    inter = mod$coefficients[1]
    #Extract slope
    slo = mod$coefficients[2]
    #Fuse every info in one row table
    dtatablesub = data.frame("1'000", sim, rad, Rsquared, inter, slo)
    names(dtatablesub) = colnames(dtaTable1_1000)
    
    #Fill the data
    dtaTable1_1000 = rbind(dtaTable1_1000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_1000) = NULL
#This is the full table small NE

#Read the data Ne 1'000 50K
dtaFROH1000050k = read.table("./LargePop/Analyses/SmallArray_PLINK_FROH.txt", header = T)[,c(1:3,8)]
dtaFROH1000050k = cbind(dtaFROH1000050k[,1], SUB_meth=rep("small", nrow(dtaFROH1000050k)), dtaFROH1000050k[,2:4])
colnames(dtaFROH1000050k) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")
#Read the data Ne 1'000 700K
dtaFROH10000700k = read.table("./LargePop/Analyses/LargeArray_PLINK_FROH.txt", header = T)[,c(1:3,8)]
dtaFROH10000700k = cbind(dtaFROH10000700k[,1], SUB_meth=rep("large", nrow(dtaFROH10000700k)), dtaFROH10000700k[,2:4])
colnames(dtaFROH10000700k) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")
#Merge both
dtaFROH10000 = rbind(dtaFROH1000050k, dtaFROH10000700k)

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_10000 = as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(dtaTable1_10000) = c("PopulationSize", "SimuID", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH10000$Simu_ID)){
  #Loop trhough array sizes
  for(rad in sort(unique(dtaFROH10000$SUB_meth))){
    #Subsample the data
    dta_sub = dtaFROH10000[dtaFROH10000$Simu_ID == sim & dtaFROH10000$SUB_meth == rad,]
    #Estimate Rsquared
    Rsquared = cor(dta_sub$Froh_array, dta_sub$Froh_TRUE_IBD_100GEN)
    #Build linear model
    mod = lm(dta_sub$Froh_array ~ dta_sub$Froh_TRUE_IBD_100GEN)
    #Extract Intercept
    inter = mod$coefficients[1]
    #Extract slope
    slo = mod$coefficients[2]
    #Fuse every info in one row table
    dtatablesub = data.frame("10'000", sim, rad, Rsquared, inter, slo)
    names(dtatablesub) = colnames(dtaTable1_10000)
    
    #Fill the data
    dtaTable1_10000 = rbind(dtaTable1_10000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_10000) = NULL
#This is the full table small NE

#Merge both small and large POP for full
dtaFullS1 = rbind(dtaTable1_1000, dtaTable1_10000)
#Write this as SUPPLEMENTARY TABLE
#write.csv(dtaFullS1, "./TABLES/FROH_SNPsArrays_PLINK100KB_bothPops_SUPPLEMENTARY_fullTABLE.csv",
#          row.names = F, quote = F)

#Create the TABLE 2 (mean Rsquared, slope & Intercept per PercSeqGen)
dtaTable1MEAN = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$ArraySize), FUN = mean)
dtaTable1MEAN = dtaTable1MEAN[order(dtaTable1MEAN$Group.1, dtaTable1MEAN$Group.2),]
dtaTable1SD = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$ArraySize), FUN = sd)
dtaTable1SD = dtaTable1SD[order(dtaTable1SD$Group.1, dtaTable1SD$Group.2),]

#Create emtpy dataframe we'll fill with mean ± sd
dtaTable1 = data.frame(matrix(nrow = 0, ncol = 5))
colnames(dtaTable1) = c("PopulationSize", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop through rows of mean & sd df
for (row in 1:nrow(dtaTable1MEAN)) {
  #Extract PopSize
  Popsize = dtaTable1MEAN$Group.1[row]
  #Extract PerSeq
  ArraySize = dtaTable1MEAN$Group.2[row]
  #Extract mean RSquared ± sd
  Rsquared = paste(round(dtaTable1MEAN$Rsquared[row], digits = 3), round(dtaTable1SD$Rsquared[row], digits = 3), sep = " ± ")
  #Extract mean Intercept ± sd
  Intercept = paste(round(dtaTable1MEAN$Intercept[row], digits = 3), round(dtaTable1SD$Intercept[row], digits = 3), sep = " ± ")
  #Extract mean Slope ± sd
  Slope = paste(round(dtaTable1MEAN$Slope[row], digits = 3), round(dtaTable1SD$Slope[row], digits = 3), sep = " ± ")
  #Fuse in one row
  rowsub = data.frame(Popsize, ArraySize, Rsquared, Intercept, Slope)
  #Set colnames for rbind
  colnames(rowsub) = colnames(dtaTable1)
  #Rbind
  dtaTable1 = rbind(dtaTable1, rowsub)
}


#This is the complete TABLE 1
write.csv(dtaTable1, "./TABLES/FROH_SNPsArrays_PLINK100KB_bothPops_TABLE5.csv",
          row.names = F, quote = F)

########################################################################################
#### TABLE 6: SNP-ARRAYs RZooROH MEAN Rsquared, slope & Intercept among simulations ####
########################################################################################

#Read the data Ne 1'000 50K
dtaFROH100050k = read.table("./SmallPop/Analyses/SmallArray_RZooRoH_FROH.txt", header = T)[,c(1:2,4,6)]
dtaFROH100050k = cbind(dtaFROH100050k[,1], SUB_meth=rep("small", nrow(dtaFROH100050k)), dtaFROH100050k[,2:4])
colnames(dtaFROH100050k) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")
#Read the data Ne 1'000 700K
dtaFROH1000700k = read.table("./SmallPop/Analyses/LargeArray_RZooRoH_FROH.txt", header = T)[,c(1:2,4,6)]
dtaFROH1000700k = cbind(dtaFROH1000700k[,1], SUB_meth=rep("large", nrow(dtaFROH1000700k)), dtaFROH1000700k[,2:4])
colnames(dtaFROH1000700k) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")
#Merge both
dtaFROH1000 = rbind(dtaFROH100050k, dtaFROH1000700k)

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_1000 = as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(dtaTable1_1000) = c("PopulationSize", "SimuID", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH1000$Simu_ID)){
  #Loop trhough array sizes
  for(rad in sort(unique(dtaFROH1000$SUB_meth))){
    #Subsample the data
    dta_sub = dtaFROH1000[dtaFROH1000$Simu_ID == sim & dtaFROH1000$SUB_meth == rad,]
    #Estimate Rsquared
    Rsquared = cor(dta_sub$Froh_array, dta_sub$Froh_TRUE_IBD_100GEN)
    #Build linear model
    mod = lm(dta_sub$Froh_array ~ dta_sub$Froh_TRUE_IBD_100GEN)
    #Extract Intercept
    inter = mod$coefficients[1]
    #Extract slope
    slo = mod$coefficients[2]
    #Fuse every info in one row table
    dtatablesub = data.frame("1'000", sim, rad, Rsquared, inter, slo)
    names(dtatablesub) = colnames(dtaTable1_1000)
    
    #Fill the data
    dtaTable1_1000 = rbind(dtaTable1_1000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_1000) = NULL
#This is the full table small NE

#Read the data Ne 1'000 50K
dtaFROH1000050k = read.table("./LargePop/Analyses/SmallArray_RZooRoH_FROH.txt", header = T)[,c(1:2,4,6)]
dtaFROH1000050k = cbind(dtaFROH1000050k[,1], SUB_meth=rep("small", nrow(dtaFROH1000050k)), dtaFROH1000050k[,2:4])
colnames(dtaFROH1000050k) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")
#Read the data Ne 1'000 700K
dtaFROH10000700k = read.table("./LargePop/Analyses/LargeArray_RZooRoH_FROH.txt", header = T)[,c(1:2,4,6)]
dtaFROH10000700k = cbind(dtaFROH10000700k[,1], SUB_meth=rep("large", nrow(dtaFROH10000700k)), dtaFROH10000700k[,2:4])
colnames(dtaFROH10000700k) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")
#Merge both
dtaFROH10000 = rbind(dtaFROH1000050k, dtaFROH10000700k)

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_10000 = as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(dtaTable1_10000) = c("PopulationSize", "SimuID", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH10000$Simu_ID)){
  #Loop trhough array sizes
  for(rad in sort(unique(dtaFROH10000$SUB_meth))){
    #Subsample the data
    dta_sub = dtaFROH10000[dtaFROH10000$Simu_ID == sim & dtaFROH10000$SUB_meth == rad,]
    #Estimate Rsquared
    Rsquared = cor(dta_sub$Froh_array, dta_sub$Froh_TRUE_IBD_100GEN)
    #Build linear model
    mod = lm(dta_sub$Froh_array ~ dta_sub$Froh_TRUE_IBD_100GEN)
    #Extract Intercept
    inter = mod$coefficients[1]
    #Extract slope
    slo = mod$coefficients[2]
    #Fuse every info in one row table
    dtatablesub = data.frame("10'000", sim, rad, Rsquared, inter, slo)
    names(dtatablesub) = colnames(dtaTable1_10000)
    
    #Fill the data
    dtaTable1_10000 = rbind(dtaTable1_10000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_10000) = NULL
#This is the full table small NE

#Merge both small and large POP for full
dtaFullS1 = rbind(dtaTable1_1000, dtaTable1_10000)
#Write this as SUPPLEMENTARY TABLE
#write.csv(dtaFullS1, "./TABLES/FROH_SNPsArrays_PLINK100KB_bothPops_SUPPLEMENTARY_fullTABLE.csv",
#          row.names = F, quote = F)

#Create the TABLE 2 (mean Rsquared, slope & Intercept per PercSeqGen)
dtaTable1MEAN = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$ArraySize), FUN = mean)
dtaTable1MEAN = dtaTable1MEAN[order(dtaTable1MEAN$Group.1, dtaTable1MEAN$Group.2),]
dtaTable1SD = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$ArraySize), FUN = sd)
dtaTable1SD = dtaTable1SD[order(dtaTable1SD$Group.1, dtaTable1SD$Group.2),]

#Create emtpy dataframe we'll fill with mean ± sd
dtaTable1 = data.frame(matrix(nrow = 0, ncol = 5))
colnames(dtaTable1) = c("PopulationSize", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop through rows of mean & sd df
for (row in 1:nrow(dtaTable1MEAN)) {
  #Extract PopSize
  Popsize = dtaTable1MEAN$Group.1[row]
  #Extract PerSeq
  ArraySize = dtaTable1MEAN$Group.2[row]
  #Extract mean RSquared ± sd
  Rsquared = paste(round(dtaTable1MEAN$Rsquared[row], digits = 3), round(dtaTable1SD$Rsquared[row], digits = 3), sep = " ± ")
  #Extract mean Intercept ± sd
  Intercept = paste(round(dtaTable1MEAN$Intercept[row], digits = 3), round(dtaTable1SD$Intercept[row], digits = 3), sep = " ± ")
  #Extract mean Slope ± sd
  Slope = paste(round(dtaTable1MEAN$Slope[row], digits = 3), round(dtaTable1SD$Slope[row], digits = 3), sep = " ± ")
  #Fuse in one row
  rowsub = data.frame(Popsize, ArraySize, Rsquared, Intercept, Slope)
  #Set colnames for rbind
  colnames(rowsub) = colnames(dtaTable1)
  #Rbind
  dtaTable1 = rbind(dtaTable1, rowsub)
}


#This is the complete TABLE 1
write.csv(dtaTable1, "./TABLES/FROH_SNPsArrays_RZooROH_bothPops_TABLE6.csv",
          row.names = F, quote = F)

###############################################################################
#### TABLE 7: WGS PLINK MEAN Rsquared, slope & Intercept among simulations ####
###############################################################################

#Read the data Ne 1'000 50K
dtaFROH1000 = read.table("./SmallPop/Analyses/WGS_PLINK_FROH.txt", header = T)[,c(1:3,6)]
dtaFROH1000 = cbind(dtaFROH1000[,1], SUB_meth=rep("wgs", nrow(dtaFROH1000)), dtaFROH1000[,2:4])
colnames(dtaFROH1000) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_1000 = as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(dtaTable1_1000) = c("PopulationSize", "SimuID", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH1000$Simu_ID)){
  #Loop trhough array sizes
  for(rad in sort(unique(dtaFROH1000$SUB_meth))){
    #Subsample the data
    dta_sub = dtaFROH1000[dtaFROH1000$Simu_ID == sim & dtaFROH1000$SUB_meth == rad,]
    #Estimate Rsquared
    Rsquared = cor(dta_sub$Froh_array, dta_sub$Froh_TRUE_IBD_100GEN)
    #Build linear model
    mod = lm(dta_sub$Froh_array ~ dta_sub$Froh_TRUE_IBD_100GEN)
    #Extract Intercept
    inter = mod$coefficients[1]
    #Extract slope
    slo = mod$coefficients[2]
    #Fuse every info in one row table
    dtatablesub = data.frame("1'000", sim, rad, Rsquared, inter, slo)
    names(dtatablesub) = colnames(dtaTable1_1000)
    
    #Fill the data
    dtaTable1_1000 = rbind(dtaTable1_1000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_1000) = NULL
#This is the full table small NE

#Read the data Ne 1'000 50K
dtaFROH10000 = read.table("./LargePop/Analyses/WGS_PLINK_FROH.txt", header = T)[,c(1:3,6)]
dtaFROH10000 = cbind(dtaFROH10000[,1], SUB_meth=rep("wgs", nrow(dtaFROH10000)), dtaFROH10000[,2:4])
colnames(dtaFROH10000) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_10000 = as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(dtaTable1_10000) = c("PopulationSize", "SimuID", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH10000$Simu_ID)){
  #Loop trhough array sizes
  for(rad in sort(unique(dtaFROH10000$SUB_meth))){
    #Subsample the data
    dta_sub = dtaFROH10000[dtaFROH10000$Simu_ID == sim & dtaFROH10000$SUB_meth == rad,]
    #Estimate Rsquared
    Rsquared = cor(dta_sub$Froh_array, dta_sub$Froh_TRUE_IBD_100GEN)
    #Build linear model
    mod = lm(dta_sub$Froh_array ~ dta_sub$Froh_TRUE_IBD_100GEN)
    #Extract Intercept
    inter = mod$coefficients[1]
    #Extract slope
    slo = mod$coefficients[2]
    #Fuse every info in one row table
    dtatablesub = data.frame("10'000", sim, rad, Rsquared, inter, slo)
    names(dtatablesub) = colnames(dtaTable1_10000)
    
    #Fill the data
    dtaTable1_10000 = rbind(dtaTable1_10000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_10000) = NULL
#This is the full table small NE

#Merge both small and large POP for full
dtaFullS1 = rbind(dtaTable1_1000, dtaTable1_10000)
#Write this as SUPPLEMENTARY TABLE
#write.csv(dtaFullS1, "./TABLES/FROH_SNPsArrays_PLINK100KB_bothPops_SUPPLEMENTARY_fullTABLE.csv",
#          row.names = F, quote = F)

#Create the TABLE 2 (mean Rsquared, slope & Intercept per PercSeqGen)
dtaTable1MEAN = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$ArraySize), FUN = mean)
dtaTable1MEAN = dtaTable1MEAN[order(dtaTable1MEAN$Group.1, dtaTable1MEAN$Group.2),]
dtaTable1SD = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$ArraySize), FUN = sd)
dtaTable1SD = dtaTable1SD[order(dtaTable1SD$Group.1, dtaTable1SD$Group.2),]

#Create emtpy dataframe we'll fill with mean ± sd
dtaTable1 = data.frame(matrix(nrow = 0, ncol = 5))
colnames(dtaTable1) = c("PopulationSize", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop through rows of mean & sd df
for (row in 1:nrow(dtaTable1MEAN)) {
  #Extract PopSize
  Popsize = dtaTable1MEAN$Group.1[row]
  #Extract PerSeq
  ArraySize = dtaTable1MEAN$Group.2[row]
  #Extract mean RSquared ± sd
  Rsquared = paste(round(dtaTable1MEAN$Rsquared[row], digits = 3), round(dtaTable1SD$Rsquared[row], digits = 3), sep = " ± ")
  #Extract mean Intercept ± sd
  Intercept = paste(round(dtaTable1MEAN$Intercept[row], digits = 3), round(dtaTable1SD$Intercept[row], digits = 3), sep = " ± ")
  #Extract mean Slope ± sd
  Slope = paste(round(dtaTable1MEAN$Slope[row], digits = 3), round(dtaTable1SD$Slope[row], digits = 3), sep = " ± ")
  #Fuse in one row
  rowsub = data.frame(Popsize, ArraySize, Rsquared, Intercept, Slope)
  #Set colnames for rbind
  colnames(rowsub) = colnames(dtaTable1)
  #Rbind
  dtaTable1 = rbind(dtaTable1, rowsub)
}


#This is the complete TABLE 1
write.csv(dtaTable1, "./TABLES/FROH_WGS_PLINK_bothPops_TABLE7.csv",
          row.names = F, quote = F)
#################################################################################
#### TABLE 8: WGS RZooROH MEAN Rsquared, slope & Intercept among simulations ####
#################################################################################

#Read the data Ne 1'000 50K
dtaFROH1000 = read.table("./SmallPop/Analyses/WGS_RZooRoH_FROH.txt", header = T)[,c(1:2,4,7)]
dtaFROH1000 = cbind(dtaFROH1000[,1], SUB_meth=rep("wgs", nrow(dtaFROH1000)), dtaFROH1000[,2:4])
colnames(dtaFROH1000) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_1000 = as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(dtaTable1_1000) = c("PopulationSize", "SimuID", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH1000$Simu_ID)){
  #Loop trhough array sizes
  for(rad in sort(unique(dtaFROH1000$SUB_meth))){
    #Subsample the data
    dta_sub = dtaFROH1000[dtaFROH1000$Simu_ID == sim & dtaFROH1000$SUB_meth == rad,]
    #Estimate Rsquared
    Rsquared = cor(dta_sub$Froh_array, dta_sub$Froh_TRUE_IBD_100GEN)
    #Build linear model
    mod = lm(dta_sub$Froh_array ~ dta_sub$Froh_TRUE_IBD_100GEN)
    #Extract Intercept
    inter = mod$coefficients[1]
    #Extract slope
    slo = mod$coefficients[2]
    #Fuse every info in one row table
    dtatablesub = data.frame("1'000", sim, rad, Rsquared, inter, slo)
    names(dtatablesub) = colnames(dtaTable1_1000)
    
    #Fill the data
    dtaTable1_1000 = rbind(dtaTable1_1000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_1000) = NULL
#This is the full table small NE

#Read the data Ne 1'000 50K
dtaFROH10000 = read.table("./LargePop/Analyses/WGS_RZooRoH_FROH.txt", header = T)[,c(1:2,4,7)]
dtaFROH10000 = cbind(dtaFROH10000[,1], SUB_meth=rep("wgs", nrow(dtaFROH10000)), dtaFROH10000[,2:4])
colnames(dtaFROH10000) = c("Simu_ID", "SUB_meth", "individuals", "Froh_array", "Froh_TRUE_IBD_100GEN")

#Create new dataframe that we will fill with RSquared, slope & Intercept
dtaTable1_10000 = as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(dtaTable1_10000) = c("PopulationSize", "SimuID", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop trhough SimuIDs
for(sim in unique(dtaFROH10000$Simu_ID)){
  #Loop trhough array sizes
  for(rad in sort(unique(dtaFROH10000$SUB_meth))){
    #Subsample the data
    dta_sub = dtaFROH10000[dtaFROH10000$Simu_ID == sim & dtaFROH10000$SUB_meth == rad,]
    #Estimate Rsquared
    Rsquared = cor(dta_sub$Froh_array, dta_sub$Froh_TRUE_IBD_100GEN)
    #Build linear model
    mod = lm(dta_sub$Froh_array ~ dta_sub$Froh_TRUE_IBD_100GEN)
    #Extract Intercept
    inter = mod$coefficients[1]
    #Extract slope
    slo = mod$coefficients[2]
    #Fuse every info in one row table
    dtatablesub = data.frame("10'000", sim, rad, Rsquared, inter, slo)
    names(dtatablesub) = colnames(dtaTable1_10000)
    
    #Fill the data
    dtaTable1_10000 = rbind(dtaTable1_10000, dtatablesub)
  }
}

#rm row-names
row.names(dtaTable1_10000) = NULL
#This is the full table small NE

#Merge both small and large POP for full
dtaFullS1 = rbind(dtaTable1_1000, dtaTable1_10000)
#Write this as SUPPLEMENTARY TABLE
#write.csv(dtaFullS1, "./TABLES/FROH_SNPsArrays_PLINK100KB_bothPops_SUPPLEMENTARY_fullTABLE.csv",
#          row.names = F, quote = F)

#Create the TABLE 2 (mean Rsquared, slope & Intercept per PercSeqGen)
dtaTable1MEAN = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$ArraySize), FUN = mean)
dtaTable1MEAN = dtaTable1MEAN[order(dtaTable1MEAN$Group.1, dtaTable1MEAN$Group.2),]
dtaTable1SD = aggregate(dtaFullS1[,c("Rsquared", "Intercept", "Slope")], by = list(dtaFullS1$PopulationSize, dtaFullS1$ArraySize), FUN = sd)
dtaTable1SD = dtaTable1SD[order(dtaTable1SD$Group.1, dtaTable1SD$Group.2),]

#Create emtpy dataframe we'll fill with mean ± sd
dtaTable1 = data.frame(matrix(nrow = 0, ncol = 5))
colnames(dtaTable1) = c("PopulationSize", "ArraySize", "Rsquared", "Intercept", "Slope")

#Loop through rows of mean & sd df
for (row in 1:nrow(dtaTable1MEAN)) {
  #Extract PopSize
  Popsize = dtaTable1MEAN$Group.1[row]
  #Extract PerSeq
  ArraySize = dtaTable1MEAN$Group.2[row]
  #Extract mean RSquared ± sd
  Rsquared = paste(round(dtaTable1MEAN$Rsquared[row], digits = 3), round(dtaTable1SD$Rsquared[row], digits = 3), sep = " ± ")
  #Extract mean Intercept ± sd
  Intercept = paste(round(dtaTable1MEAN$Intercept[row], digits = 3), round(dtaTable1SD$Intercept[row], digits = 3), sep = " ± ")
  #Extract mean Slope ± sd
  Slope = paste(round(dtaTable1MEAN$Slope[row], digits = 3), round(dtaTable1SD$Slope[row], digits = 3), sep = " ± ")
  #Fuse in one row
  rowsub = data.frame(Popsize, ArraySize, Rsquared, Intercept, Slope)
  #Set colnames for rbind
  colnames(rowsub) = colnames(dtaTable1)
  #Rbind
  dtaTable1 = rbind(dtaTable1, rowsub)
}


#This is the complete TABLE 1
write.csv(dtaTable1, "./TABLES/FROH_WGS_RZooROH_bothPops_TABLE8.csv",
          row.names = F, quote = F)


