##############################################################################
################################# FUNCTIONS  #################################
##############################################################################

setwd("~/PhD/ROH/Which_data/Data/REVIEWed/")

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


##################################################################################
################ FIGURE S1 FSUB VS FTRUE IBD 3 PERC RAD + ARRAYs UR ##############
##################################################################################

#Read the data FROH ne 1'000
dtaFROH.1000.100KB.RAD = read.table("./SmallPop/Analyses/RADSEQ_PLINK_FROH.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH.1000.100KB.RAD = dtaFROH.1000.100KB.RAD[dtaFROH.1000.100KB.RAD$WINDOWS.NB %in% c(180000, 240000, 600000),]

#Read the data FROH ne 10'000
dtaFROH.10000.100KB.RAD = read.table("./LargePop/Analyses/RADSEQ_PLINK_FROH.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH.10000.100KB.RAD = dtaFROH.10000.100KB.RAD[dtaFROH.10000.100KB.RAD$WINDOWS.NB %in% c(15000, 25000, 120000),]

#Read the data FROH ne 1'000
dtaFROH.1000.Rzoo.RAD = read.table("./SmallPop/Analyses/RADSEQ_RZooRoH_FROH.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH.1000.Rzoo.RAD = dtaFROH.1000.Rzoo.RAD[dtaFROH.1000.Rzoo.RAD$WINDOWS.NB %in% c(3000, 9000, 60000),]

#Read the data FROH ne 10'000
dtaFROH.10000.Rzoo.RAD = read.table("./LargePop/Analyses/RADSEQ_RZooRoH_FROH.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH.10000.Rzoo.RAD = dtaFROH.10000.Rzoo.RAD[dtaFROH.10000.Rzoo.RAD$WINDOWS.NB %in% c(250, 500, 7500),]

##Read the data Ne 1'000 50K
dtaFROH.1000.100KB.50K = read.table("./SmallPop/Analyses/SmallArray_PLINK_FROH.txt", header = T)
dtaFROH.1000.100KB.50K = cbind(dtaFROH.1000.100KB.50K, ARRAY=rep("SMALL", nrow(dtaFROH.1000.100KB.50K)))
#Read the data Ne 1'000 700K
dtaFROH.1000.100KB.700K = read.table("./SmallPop/Analyses/LargeArray_PLINK_FROH.txt", header = T)
dtaFROH.1000.100KB.700K = cbind(dtaFROH.1000.100KB.700K, ARRAY=rep("LARGE", nrow(dtaFROH.1000.100KB.700K)))
dtaFROH.1000.100KB.ARRAY = rbind(dtaFROH.1000.100KB.50K,dtaFROH.1000.100KB.700K)
rm(dtaFROH.1000.100KB.50K);rm(dtaFROH.1000.100KB.700K)
dtaFROH.1000.100KB.ARRAY$ARRAY = factor(dtaFROH.1000.100KB.ARRAY$ARRAY, levels = c("SMALL","LARGE"))

##Read the data Ne 1'000 50K
dtaFROH.10000.100KB.50K = read.table("./LargePop/Analyses/SmallArray_PLINK_FROH.txt", header = T)
dtaFROH.10000.100KB.50K = cbind(dtaFROH.10000.100KB.50K, ARRAY=rep("SMALL", nrow(dtaFROH.10000.100KB.50K)))
#Read the data Ne 1'000 700K
dtaFROH.10000.100KB.700K = read.table("./LargePop/Analyses/LargeArray_PLINK_FROH.txt", header = T)
dtaFROH.10000.100KB.700K = cbind(dtaFROH.10000.100KB.700K, ARRAY=rep("LARGE", nrow(dtaFROH.10000.100KB.700K)))
dtaFROH.10000.100KB.ARRAY = rbind(dtaFROH.10000.100KB.50K,dtaFROH.10000.100KB.700K)
rm(dtaFROH.10000.100KB.50K);rm(dtaFROH.10000.100KB.700K)
dtaFROH.10000.100KB.ARRAY$ARRAY = factor(dtaFROH.10000.100KB.ARRAY$ARRAY, levels = c("SMALL","LARGE"))

#Read the data Ne 1'000 50K
dtaFROH.1000.Rzoo.50K = read.table("./SmallPop/Analyses/SmallArray_RZooRoH_FROH.txt", header = T)
colnames(dtaFROH.1000.Rzoo.50K) = c("Simu_ID","individuals", "ARRAY", "Froh_array", "Replicate", "Froh_trueIBD100Gen", "Froh_trueIBD1000Gen")
#change model/array into "SMALL"
dtaFROH.1000.Rzoo.50K$ARRAY = "SMALL"
#Read the data Ne 1'000 700K
dtaFROH.1000.Rzoo.700K = read.table("./SmallPop/Analyses/LargeArray_RZooRoH_FROH.txt", header = T)
colnames(dtaFROH.1000.Rzoo.700K) =  c("Simu_ID","individuals", "ARRAY", "Froh_array", "Replicate", "Froh_trueIBD100Gen", "Froh_trueIBD1000Gen")
#change model/array into "SMALL"
dtaFROH.1000.Rzoo.700K$ARRAY = "LARGE"
dtaFROH.1000.Rzoo.ARRAY = rbind(dtaFROH.1000.Rzoo.50K,dtaFROH.1000.Rzoo.700K)
rm(dtaFROH.1000.Rzoo.50K);rm(dtaFROH.1000.Rzoo.700K)
dtaFROH.1000.Rzoo.ARRAY$ARRAY = factor(dtaFROH.1000.Rzoo.ARRAY$ARRAY, levels = c("SMALL", "LARGE"))

#Read the data Ne 10'000 50K
dtaFROH.10000.Rzoo.50K = read.table("./LargePop/Analyses/SmallArray_RZooRoH_FROH.txt", header = T)
colnames(dtaFROH.10000.Rzoo.50K) = c("Simu_ID","individuals", "ARRAY", "Froh_array", "Replicate", "Froh_trueIBD100Gen", "Froh_trueIBD1000Gen")
#change model/array into "SMALL"
dtaFROH.10000.Rzoo.50K$ARRAY = "SMALL"
#Read the data Ne 10'000 700K
dtaFROH.10000.Rzoo.700K = read.table("./LargePop/Analyses/LargeArray_RZooRoH_FROH.txt", header = T)
colnames(dtaFROH.10000.Rzoo.700K) =  c("Simu_ID","individuals", "ARRAY", "Froh_array", "Replicate", "Froh_trueIBD100Gen", "Froh_trueIBD1000Gen")
#change model/array into "SMALL"
dtaFROH.10000.Rzoo.700K$ARRAY = "LARGE"
dtaFROH.10000.Rzoo.ARRAY = rbind(dtaFROH.10000.Rzoo.50K,dtaFROH.10000.Rzoo.700K)
rm(dtaFROH.10000.Rzoo.50K);rm(dtaFROH.10000.Rzoo.700K)
dtaFROH.10000.Rzoo.ARRAY$ARRAY = factor(dtaFROH.10000.Rzoo.ARRAY$ARRAY, levels = c("SMALL", "LARGE"))

jpeg("~/PhD/ROH/Which_data/Manuscript/SUBMITTED_VERSION/Review_round_I/FIGURES/SUPPLEMENTARY_FIGURES/FIGURES1.jpg", width = 5000, height = 3500)

#Layout for having all plots
layout(matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5)), 5, 5, byrow = T), heights = c(.5,1,.7,1,.3), widths = c(.5,1,1.5,1,1))

#Par for margin and outer margings
par(mar = c(3.5,3.5,3.5,3.5), oma = c(.2,.2,.2,.2), xpd = T)

#### PANEL A: FROH SUB ~ FROH TRUE IBD 100 RAD PLINK ####

#Plot PLINK RADs SMALL POP EMPTY PLOT
plot(dtaFROH.1000.100KB.RAD$FROH_sub ~ dtaFROH.1000.100KB.RAD$Froh_TRUE_IBD_100GEN, pch = 19,
     ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(0,0.9), ylim = c(0,0.9), cex.axis = 4, type = "n", frame = F)

#Colors RAD SMALL NE
RADcolours = c("#82EEFD", "#0492C2", "#016064")
names(RADcolours) = sort(unique(dtaFROH.1000.100KB.RAD$WINDOWS.NB))
#Shape points SMALL NE
POINTSshapes = c(15,16,17)
names(POINTSshapes) = sort(unique(dtaFROH.1000.100KB.RAD$WINDOWS.NB))

#Loop through the RAD subset SMALL NE
for (sub in sort(unique(dtaFROH.1000.100KB.RAD$WINDOWS.NB))) {
  
  #Add the points SMALL NE
  points(dtaFROH.1000.100KB.RAD$FROH_sub[dtaFROH.1000.100KB.RAD$WINDOWS.NB == sub] ~
           dtaFROH.1000.100KB.RAD$Froh_TRUE_IBD_100GEN[dtaFROH.1000.100KB.RAD$WINDOWS.NB == sub],
         col = add.alpha(RADcolours[names(RADcolours) == sub],.5), pch = POINTSshapes[names(POINTSshapes) == sub], cex = 3)
  
}

#Colors RAD LARGE NE
RADcolours = c("#FFBF00", "#FD6A02", "#FDA172")
names(RADcolours) = sort(unique(dtaFROH.10000.100KB.RAD$WINDOWS.NB))
#Shape points LARGE NE
POINTSshapes = c(1,4,8)
names(POINTSshapes) = sort(unique(dtaFROH.10000.100KB.RAD$WINDOWS.NB))

#Loop through the RAD subset LARGE NE
for (sub in sort(unique(dtaFROH.10000.100KB.RAD$WINDOWS.NB))) {
  
  #Add the points LARGE NE
  points(dtaFROH.10000.100KB.RAD$FROH_sub[dtaFROH.10000.100KB.RAD$WINDOWS.NB == sub] ~
           dtaFROH.10000.100KB.RAD$Froh_TRUE_IBD_100GEN[dtaFROH.10000.100KB.RAD$WINDOWS.NB == sub],
         col = add.alpha(RADcolours[names(RADcolours) == sub],.5), pch = POINTSshapes[names(POINTSshapes) == sub], cex = 3)
  
}

## ADD AXES
axis(1,xlab = NULL, at=c(0:9)*0.1, padj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, lwd = 3)
axis(2,xlab = NULL, at=c(0:9)*0.1, hadj = 1.5, cex.axis= 10, lwd.ticks = 3, las = 1, tck = -0.03, lwd = 3)

#Add axes titles
mtext(expression(F["IBD"]), side = 1, outer = FALSE, line = 30, cex = 8)
mtext(expression(F["ROH RAD"]), side = 2, outer = FALSE, line = 30, cex = 8)

#Add perfect line (must be below everything else)
abline(0,1, lwd = 4)

#Add panel legend
mtext("A", side = 3, outer = FALSE, line = 20, cex = 10, at = -.3)
#Add panel lab
mtext("PLINK", side = 3, outer = FALSE, line = 20, cex = 8, at = .35)

#### PANEL B: FROH SUB ~ FROH TRUE IBD 100 RAD RZooRoH ####

#Plot PLINK RADs SMALL POP EMPTY PLOT
plot(dtaFROH.1000.Rzoo.RAD$FROH_sub ~ dtaFROH.1000.Rzoo.RAD$Froh_TRUE_IBD_100GEN, pch = 19,
     ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(0,0.75), ylim = c(0,0.75), cex.axis = 4, type = "n", frame = F)

#Colors RAD SMALL NE
RADcolours = c("#82EEFD", "#0492C2", "#016064")
names(RADcolours) = sort(unique(dtaFROH.1000.Rzoo.RAD$WINDOWS.NB))
#Shape points SMALL NE
POINTSshapes = c(15,16,17)
names(POINTSshapes) = sort(unique(dtaFROH.1000.Rzoo.RAD$WINDOWS.NB))

#Loop through the RAD subset SMALL NE
for (sub in sort(unique(dtaFROH.1000.Rzoo.RAD$WINDOWS.NB))) {
  
  #Add the points SMALL NE
  points(dtaFROH.1000.Rzoo.RAD$FROH_sub[dtaFROH.1000.Rzoo.RAD$WINDOWS.NB == sub] ~
           dtaFROH.1000.Rzoo.RAD$Froh_TRUE_IBD_100GEN[dtaFROH.1000.Rzoo.RAD$WINDOWS.NB == sub],
         col = add.alpha(RADcolours[names(RADcolours) == sub],.5), pch = POINTSshapes[names(POINTSshapes) == sub], cex = 3)
  
}

#Colors RAD LARGE NE
RADcolours = c("#FFBF00", "#FD6A02", "#FDA172")
names(RADcolours) = sort(unique(dtaFROH.10000.Rzoo.RAD$WINDOWS.NB))
#Shape points LARGE NE
POINTSshapes = c(1,4,8)
names(POINTSshapes) = sort(unique(dtaFROH.10000.Rzoo.RAD$WINDOWS.NB))

#Loop through the RAD subset LARGE NE
for (sub in sort(unique(dtaFROH.10000.Rzoo.RAD$WINDOWS.NB))) {
  
  #Add the points LARGE NE
  points(dtaFROH.10000.Rzoo.RAD$FROH_sub[dtaFROH.10000.Rzoo.RAD$WINDOWS.NB == sub] ~
           dtaFROH.10000.Rzoo.RAD$Froh_TRUE_IBD_100GEN[dtaFROH.10000.Rzoo.RAD$WINDOWS.NB == sub],
         col = add.alpha(RADcolours[names(RADcolours) == sub],.5), pch = POINTSshapes[names(POINTSshapes) == sub], cex = 3)
  
}

## ADD AXES
axis(1,xlab = NULL, at=c(0:8)*0.1, padj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, lwd = 3)
axis(2,xlab = NULL, at=c(0:8)*0.1, hadj = 1.5, cex.axis= 10, lwd.ticks = 3, las = 1, tck = -0.03, lwd = 3)

#Add axes titles
mtext(expression(F["IBD"]), side = 1, outer = FALSE, line = 30, cex = 8)
mtext(expression(F["HBD RAD"]), side = 2, outer = FALSE, line = 30, cex = 8)

#Add perfect line (must be below everything else)
abline(0,1, lwd = 4)

#Add panel lab
mtext("B", side = 3, outer = FALSE, line = 20, cex = 10, at = -.3)
#Add panel lab
mtext("RZooRoH", side = 3, outer = FALSE, line = 20, cex = 8, at = .35)

#### PANNEL C: ARRAYs FROH PLINK ####

#Plot PLINK ARRAYs SMALL POP EMPTY PLOT
plot(dtaFROH.1000.100KB.ARRAY$FROH_sub ~ dtaFROH.1000.100KB.ARRAY$Froh_TRUE_IBD_100GEN, pch = 19,
     ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(0,0.7), ylim = c(0,0.7), cex.axis = 4, type = "n", frame = F)

#Add the points SMALL NE
points(dtaFROH.1000.100KB.ARRAY$FROH_sub ~ dtaFROH.1000.100KB.ARRAY$Froh_TRUE_IBD_100GEN, col = add.alpha(c("#2832C2", "#0A1172")[dtaFROH.1000.100KB.ARRAY$ARRAY],1),
       pch = c(15,16)[dtaFROH.1000.100KB.ARRAY$ARRAY], cex = 3)

#Add the points LARGE NE
points(dtaFROH.10000.100KB.ARRAY$FROH_sub ~ dtaFROH.10000.100KB.ARRAY$Froh_TRUE_IBD_100GEN, col = add.alpha(c("#CC7722", "#793802")[dtaFROH.10000.100KB.ARRAY$ARRAY],1),
       pch = c(17,18)[dtaFROH.10000.100KB.ARRAY$ARRAY], cex = 3)

## ADD AXES
axis(1,xlab = NULL, at=c(0:7)*0.1, padj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, lwd = 3)
axis(2,xlab = NULL, at=c(0:7)*0.1, hadj = 1.5, cex.axis= 10, lwd.ticks = 3, las = 1, tck = -0.03, lwd = 3)

#Add axes titles
mtext(expression(F["IBD"]), side = 1, outer = FALSE, line = 26, cex = 8)
mtext(expression(F["ROH ARRAYs"]), side = 2, outer = FALSE, line = 26, cex = 8)

#Add perfect line (must be after everything else)
abline(0,1, lwd = 4)

#Add panel legend
mtext("C", side = 3, outer = FALSE, line = 20, cex = 10, at = -.3)
#Add panel lab
mtext("PLINK", side = 3, outer = FALSE, line = 20, cex = 8, at = .35)

#### PANNEL D: ARRAYs FROH RZooRoH ####

#Plot RZooRoH ARRAYs SMALL POP EMPTY PLOT
plot(dtaFROH.1000.Rzoo.ARRAY$Froh_array ~ dtaFROH.1000.Rzoo.ARRAY$Froh_trueIBD100Gen, pch = 19,
     ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(0,0.75), ylim = c(0,0.75), cex.axis = 4, type = "n", frame = F)

#Add the points SMALL NE
points(dtaFROH.1000.Rzoo.ARRAY$Froh_array ~ dtaFROH.1000.Rzoo.ARRAY$Froh_trueIBD100Gen, col = add.alpha(c("#2832C2", "#0A1172")[dtaFROH.1000.Rzoo.ARRAY$ARRAY],1),
       pch = c(15,16)[dtaFROH.1000.Rzoo.ARRAY$ARRAY], cex = 3)

#Add the points LARGE NE
points(dtaFROH.10000.Rzoo.ARRAY$Froh_array ~ dtaFROH.10000.Rzoo.ARRAY$Froh_trueIBD100Gen, col = add.alpha(c("#CC7722", "#793802")[dtaFROH.10000.Rzoo.ARRAY$ARRAY],1),
       pch = c(17,18)[dtaFROH.1000.Rzoo.ARRAY$ARRAY], cex = 3)

## ADD AXES
axis(1,xlab = NULL, at=c(0:8)*0.1, padj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, lwd = 3)
axis(2,xlab = NULL, at=c(0:8)*0.1, hadj = 1.5, cex.axis= 10, lwd.ticks = 3, las = 1, tck = -0.03, lwd = 3)
#Add axes titles
mtext(expression(F["IBD"]), side = 1, outer = FALSE, line = 26, cex = 8)
mtext(expression(F["HBD ARRAYs"]), side = 2, outer = FALSE, line = 26, cex = 8)

#Add perfect line (must be after everything else)
abline(0,1, lwd = 4)

#Add panel legend
mtext("D", side = 3, outer = FALSE, line = 20, cex = 10, at = -.3)
#Add panel lab
mtext("RZooRoH", side = 3, outer = FALSE, line = 20, cex = 8, at = .35)

dev.off()

#####################################################################
################ FIGURE S2 FHOM VS FTRUE IBD BOTH Pops ##############
#####################################################################

#Read the data
dtaNe1000PLINK = read.table("./SmallPop/Analyses/RADSEQ_PLINK_FHOM.txt", header = T)
dtaNe1000Rzoo = read.table("./SmallPop/Analyses/RADSEQ_RZooRoH_FHOM.txt", header = T)
dtaNe10000PLINK = read.table("./LargePop/Analyses/RADSEQ_PLINK_FHOM.txt", header = T)
dtaNe10000Rzoo = read.table("./LargePop/Analyses/RADSEQ_RZooRoH_FHOM.txt", header = T)

#Open PDF
jpeg("~/PhD/ROH/Which_data/Manuscript/SUBMITTED_VERSION/Review_round_I/FIGURES/SUPPLEMENTARY_FIGURES/FIGURES2.jpg", width = 5500, height = 3500)

layout(matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5)), 5, 5, byrow = T), heights = c(.5,1,.7,1,.3), widths = c(.5,1,1.5,1,1.5))

#### Pannel A FHOM with SUBsample from PLINK SMALL POP #####

#Create empty graph
plot(dtaNe1000PLINK$FHOM_sub ~ dtaNe1000PLINK$Froh_TRUE_IBD_100GEN, type = 'n', xlim = c(-0.1,0.9), ylim = c(-0.1,0.9), ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', frame = F)
RADcolours = c("#82EEFD", "#0492C2", "#016064")
names(RADcolours) = c(180000,240000,600000)
POINTSshapes = c(15,16,17)
names(POINTSshapes) = c(180000,240000,600000)

#Loop through sub from PLINK SMALL NE
for (sub in c(180000,240000,600000)) {
  
  #Subsample the data
  dta_sub = dtaNe1000PLINK[dtaNe1000PLINK$WINDOWS.NB == sub,]
  
  #Fill the empty plot
  points(dta_sub$FHOM_sub ~ dta_sub$Froh_TRUE_IBD_100GEN, pch = POINTSshapes[names(POINTSshapes) == sub], col = add.alpha(RADcolours[names(RADcolours) == sub],0.6))
  
}

#Add the one-to-one line
abline(0,1, lwd = 4)

#Add x axis
axis(1,xlab = NULL, at=c(-1:9)*0.1, padj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, lwd = 3)
#Add y axis
axis(2,xlab = NULL, at=c(-1:9)*0.1, hadj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, las = 1, lwd = 3)

#Add the axes names
mtext(expression(F[HOM]), side = 2, outer = FALSE, line = 25, cex = 8, at = .35)
mtext(expression(F["IBD"]), side = 1, outer = FALSE, line = 25, cex = 8, at = .35)
mtext("A",side = 3, outer = FALSE, line = 20, cex = 10, at = -.6)

#### Pannel B FHOM with SUBsample from RZOOROH SMALL POP #####

#Create empty graph
plot(dtaNe1000Rzoo$FHOM_sub ~ dtaNe1000Rzoo$Froh_TRUE_IBD_100GEN, type = 'n', xlim = c(-.1,0.9), ylim = c(-.3,0.9), ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', frame = F)
RADcolours = c("#82EEFD", "#0492C2", "#016064")
names(RADcolours) = c(3000, 9000, 120000)
POINTSshapes = c(15,16,17)
names(POINTSshapes) = c(3000, 9000, 120000)

#Loop through sub from PLINK SMALL NE
for (sub in c(3000, 9000, 60000)) {
  
  #Subsample the data
  dta_sub = dtaNe1000Rzoo[dtaNe1000Rzoo$WINDOWS.NB == sub,]
  
  #Fill the empty plot
  points(dta_sub$FHOM_sub ~ dta_sub$Froh_TRUE_IBD_100GEN, pch = POINTSshapes[names(POINTSshapes) == sub], col = add.alpha(RADcolours[names(RADcolours) == sub],0.6))
  
}

#Add the one-to-one line
abline(0,1, lwd = 4)

#Add x axis
axis(1,xlab = NULL, at=c(-1:9)*0.1, padj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, lwd = 3)
#Add y axis
axis(2,xlab = NULL, at=c(-4:9)*0.1, hadj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, las = 1, lwd = 3)

#Add the axes names
mtext(expression(F[HOM]), side = 2, outer = FALSE, line = 25, cex = 8, at = .3)
mtext(expression(F["IBD"]), side = 1, outer = FALSE, line = 25, cex = 8, at = .35)
mtext("B",side = 3, outer = FALSE, line = 20, cex = 10, at = -.6)

#### Pannel C FHOM with SUBsample from PLINK LARGE POP #####

#Create empty graph
plot(dtaNe10000PLINK$FHOM_sub ~ dtaNe10000PLINK$Froh_TRUE_IBD_100GEN, type = 'n', xlim = c(-.1,0.9), ylim = c(-.1,0.9), ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', frame = F)
RADcolours = c("#FFBF00", "#FD6A02", "#FDA172")
names(RADcolours) = c(15000, 25000, 60000)
POINTSshapes = c(1,4,8)
names(POINTSshapes) = c(15000, 25000, 60000)

#Loop through sub from PLINK SMALL NE
for (sub in c(15000, 25000, 60000)) {
  
  #Subsample the data
  dta_sub = dtaNe10000PLINK[dtaNe10000PLINK$WINDOWS.NB == sub,]
  
  #Fill the empty plot
  points(dta_sub$FHOM_sub ~ dta_sub$Froh_TRUE_IBD_100GEN, pch = POINTSshapes[names(POINTSshapes) == sub], col = add.alpha(RADcolours[names(RADcolours) == sub],0.6))
  
}

#Add the one-to-one line
abline(0,1, lwd = 4)

#Add x axis
axis(1,xlab = NULL, at=c(-4:9)*0.1, padj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, lwd = 3)
#Add y axis
axis(2,xlab = NULL, at=c(-4:9)*0.1, hadj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, las = 1, lwd = 3)

#Add the axes names
mtext(expression(F[HOM]), side = 2, outer = FALSE, line = 25, cex = 8, at = .35)
mtext(expression(F["IBD"]), side = 1, outer = FALSE, line = 25, cex = 8, at = .35)
mtext("C",side = 3, outer = FALSE, line = 20, cex = 10, at = -.6)

#### Pannel D FHOM with SUBsample from RZOOROH LARGE POP #####

#Create empty graph
plot(dtaNe10000Rzoo$FHOM_sub ~ dtaNe10000Rzoo$Froh_TRUE_IBD_100GEN, type = 'n', xlim = c(-.1,0.9), ylim = c(-.4,0.9), ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', frame = F)
RADcolours = c("#FFBF00", "#FD6A02", "#FDA172")
names(RADcolours) = c(250, 500, 7500)
POINTSshapes = c(1,4,8)
names(POINTSshapes) = c(250, 500, 7500)

#Loop through sub from PLINK SMALL NE
for (sub in c(250, 500, 7500)) {
  
  #Subsample the data
  dta_sub = dtaNe10000Rzoo[dtaNe10000Rzoo$WINDOWS.NB == sub,]
  
  #Fill the empty plot
  points(dta_sub$FHOM_sub ~ dta_sub$Froh_TRUE_IBD_100GEN, pch = POINTSshapes[names(POINTSshapes) == sub], col = add.alpha(RADcolours[names(RADcolours) == sub],0.6))
  
}

#Add the one-to-one line
abline(0,1, lwd = 4)

#Add x axis
axis(1,xlab = NULL, at=c(-1:9)*0.1, padj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, lwd = 3)
#Add y axis
axis(2,xlab = NULL, at=c(-4:9)*0.1, hadj = 1.5, cex.axis= 10, lwd.ticks = 3, srt = 45, tck = -0.03, las = 1, lwd = 3)

#Add the axes names
mtext(expression(F[HOM]), side = 2, outer = FALSE, line = 25, cex = 8, at = .3)
mtext(expression(F["IBD"]), side = 1, outer = FALSE, line = 25, cex = 8, at = .35)
mtext("D",side = 3, outer = FALSE, line = 20, cex = 10, at = -.6)

dev.off()

####################################################################
############# FIGURE S3 ROHs DIST. 1000 GEN AGO TRUE IBD ###########
####################################################################

#### Read the data ROHsDist PLINK100KB Ne 1'000 ####

dtaROHSDIST.1000.PLINK.RAD = read.table("./SmallPop/Analyses/RADSEQ_PLINK_ROHsDistribution.txt", header = T)
#Take only the perc we are interested in
dtaROHSDIST.1000.PLINK.RAD = dtaROHSDIST.1000.PLINK.RAD[dtaROHSDIST.1000.PLINK.RAD$WINDOWS.NB %in% c("180000", "240000", "600000"),]
colnames(dtaROHSDIST.1000.PLINK.RAD) = c("SimuID","SUB_method","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.RAD = aggregate(dtaROHSDIST.1000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.RAD$SUB_method), FUN = mean)
#Calculate sd
ROHsDist.SD.1000.PLINK.RAD = as.data.frame(as.matrix(aggregate(dtaROHSDIST.1000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.RAD$SUB_method), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.RAD = merge(dtaROHSDIST.plot.1000.PLINK.RAD, ROHsDist.SD.1000.PLINK.RAD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.RAD) = c("CLASS", "PERC", "x", "sd")
#Pass sd to numeric
dtaROHSDIST.plot.1000.PLINK.RAD$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.RAD$sd))
#Order the df
dtaROHSDIST.plot.1000.PLINK.RAD = dtaROHSDIST.plot.1000.PLINK.RAD[order(dtaROHSDIST.plot.1000.PLINK.RAD$CLASS, dtaROHSDIST.plot.1000.PLINK.RAD$PERC),]
dtaROHSDIST.plot.1000.PLINK.RAD$PERC = factor(dtaROHSDIST.plot.1000.PLINK.RAD$PERC)
#rm RAD and SD
rm(dtaROHSDIST.1000.PLINK.RAD); rm(ROHsDist.SD.1000.PLINK.RAD)

#Read the ne 1'000 50k array PLINK
dtaROHSDIST.1000.PLINK.50k = read.table("./SmallPop/Analyses/SmallArray_PLINK_ROHsDistributions.txt", header = T)
#Read the ne 1'000 700k array PLINK
dtaROHSDIST.1000.PLINK.700k = read.table("./SmallPop/Analyses/LargeArray_PLINK_ROHsDistributions.txt", header = T)
dtaROHSDIST.1000.PLINK.ARRAY = unique(rbind(dtaROHSDIST.1000.PLINK.50k,dtaROHSDIST.1000.PLINK.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.ARRAY = aggregate(dtaROHSDIST.1000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.ARRAY$SUB_TECH), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.ARRAY = aggregate(dtaROHSDIST.1000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.ARRAY$SUB_TECH, SimuIDs=dtaROHSDIST.1000.PLINK.ARRAY$SimuID, REP=dtaROHSDIST.1000.PLINK.ARRAY$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.ARRAY = merge(dtaROHSDIST.plot.1000.PLINK.ARRAY, ROHsDist.SD.1000.PLINK.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.ARRAY = dtaROHSDIST.plot.1000.PLINK.ARRAY[order(dtaROHSDIST.plot.1000.PLINK.ARRAY$CLASS, dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC),]
dtaROHSDIST.plot.1000.PLINK.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.50k);rm(dtaROHSDIST.1000.PLINK.700k);rm(dtaROHSDIST.1000.PLINK.ARRAY);rm(ROHsDist.MeanforCI.1000.PLINK.ARRAY);rm(ROHsDist.SD.1000.PLINK.ARRAY)

#Read the ne 1'000 WGS PLINK
dtaROHSDIST.1000.PLINK.WGS = read.table("./SmallPop/Analyses//WGS_PLINK_ROHsDistribution.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.WGS = aggregate(dtaROHSDIST.1000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.WGS$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.WGS$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.WGS$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.WGS$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.WGS = aggregate(dtaROHSDIST.1000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.WGS$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.PLINK.WGS$SimID, REP=dtaROHSDIST.1000.PLINK.WGS$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.WGS = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.WGS$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.WGS$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.WGS$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.WGS = merge(dtaROHSDIST.plot.1000.PLINK.WGS, ROHsDist.SD.1000.PLINK.WGS, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.WGS) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.WGS = dtaROHSDIST.plot.1000.PLINK.WGS[order(dtaROHSDIST.plot.1000.PLINK.WGS$CLASS, dtaROHSDIST.plot.1000.PLINK.WGS$PERC),]
dtaROHSDIST.plot.1000.PLINK.WGS$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.WGS$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.WGS);rm(ROHsDist.MeanforCI.1000.PLINK.WGS);rm(ROHsDist.SD.1000.PLINK.WGS)

#Read the ne 1'000 TRUE IBD 100GEN PLINK
dtaROHSDIST.1000.PLINK.TRUEIBD100 = read.table("./SmallPop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD100$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD100$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.PLINK.TRUEIBD100$SimID, REP=dtaROHSDIST.1000.PLINK.TRUEIBD100$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.TRUEIBD100 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100 = merge(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100, ROHsDist.SD.1000.PLINK.TRUEIBD100, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100 = dtaROHSDIST.plot.1000.PLINK.TRUEIBD100[order(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$CLASS, dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC),]
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.TRUEIBD100);rm(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100);rm(ROHsDist.SD.1000.PLINK.TRUEIBD100)

#Read the ne 1'000 TRUE IBD 1000GEN PLINK
dtaROHSDIST.1000.PLINK.TRUEIBD1000 = read.table("./SmallPop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD1000$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD1000$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.PLINK.TRUEIBD1000$SimID, REP=dtaROHSDIST.1000.PLINK.TRUEIBD1000$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.TRUEIBD1000 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000 = merge(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000, ROHsDist.SD.1000.PLINK.TRUEIBD1000, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000 = dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000[order(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$CLASS, dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC),]
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.TRUEIBD1000);rm(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000);rm(ROHsDist.SD.1000.PLINK.TRUEIBD1000)

#Merge ALL and remove duplicated WGS lines
dtaROHSDIST.1000.PLINK = unique(rbind(dtaROHSDIST.plot.1000.PLINK.RAD, dtaROHSDIST.plot.1000.PLINK.ARRAY, dtaROHSDIST.plot.1000.PLINK.WGS, dtaROHSDIST.plot.1000.PLINK.TRUEIBD100, dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000))
#Change factor levels
dtaROHSDIST.1000.PLINK$PERC = factor(dtaROHSDIST.1000.PLINK$PERC, levels = c("180000", "240000", "600000","SMALL_ARRAY","LARGE_ARRAY","WGS","TRUE_IBD_100", "TRUE_IBD_1000"))
#Order the df
dtaROHSDIST.1000.PLINK = dtaROHSDIST.1000.PLINK[order(dtaROHSDIST.1000.PLINK$CLASS, dtaROHSDIST.1000.PLINK$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.1000.PLINK.RAD); rm(dtaROHSDIST.plot.1000.PLINK.ARRAY,dtaROHSDIST.plot.1000.PLINK.WGS,dtaROHSDIST.plot.1000.PLINK.TRUEIBD100,dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000)
#Divide TRUE IBD by 1000 (because in basepair while PLINK in kb)
dtaROHSDIST.1000.PLINK$x[dtaROHSDIST.1000.PLINK$PERC %in% c("TRUE_IBD_100", "TRUE_IBD_1000")] = dtaROHSDIST.1000.PLINK$x[dtaROHSDIST.1000.PLINK$PERC %in% c("TRUE_IBD_100", "TRUE_IBD_1000")]/1000

#### Read the data ROHsDist PLINK100KB Ne 1'000 ####

dtaROHSDIST.10000.PLINK.RAD = read.table("./LargePop/Analyses/RADSEQ_PLINK_ROHsDistribution.txt", header = T)
#Take only the perc we are interested in
dtaROHSDIST.10000.PLINK.RAD = dtaROHSDIST.10000.PLINK.RAD[dtaROHSDIST.10000.PLINK.RAD$WINDOWS.NB %in% c("15000", "25000", "60000"),]
colnames(dtaROHSDIST.10000.PLINK.RAD) = c("SimuID","SUB_method","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.RAD = aggregate(dtaROHSDIST.10000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.RAD$SUB_method), FUN = mean)
#Calculate sd
ROHsDist.SD.10000.PLINK.RAD = as.data.frame(as.matrix(aggregate(dtaROHSDIST.10000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.RAD$SUB_method), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.RAD = merge(dtaROHSDIST.plot.10000.PLINK.RAD, ROHsDist.SD.10000.PLINK.RAD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.RAD) = c("CLASS", "PERC", "x", "sd")
#Pass sd to numeric
dtaROHSDIST.plot.10000.PLINK.RAD$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.RAD$sd))
#Order the df
dtaROHSDIST.plot.10000.PLINK.RAD = dtaROHSDIST.plot.10000.PLINK.RAD[order(dtaROHSDIST.plot.10000.PLINK.RAD$CLASS, dtaROHSDIST.plot.10000.PLINK.RAD$PERC),]
dtaROHSDIST.plot.10000.PLINK.RAD$PERC = factor(dtaROHSDIST.plot.10000.PLINK.RAD$PERC)
#rm RAD and SD
rm(dtaROHSDIST.10000.PLINK.RAD); rm(ROHsDist.SD.10000.PLINK.RAD)

#Read the ne 1'000 50k array PLINK
dtaROHSDIST.10000.PLINK.50k = read.table("./LargePop/Analyses/SmallArray_PLINK_ROHsDistributions.txt", header = T)
#Read the ne 1'000 700k array PLINK
dtaROHSDIST.10000.PLINK.700k = read.table("./LargePop/Analyses/LargeArray_PLINK_ROHsDistributions.txt", header = T)
dtaROHSDIST.10000.PLINK.ARRAY = unique(rbind(dtaROHSDIST.10000.PLINK.50k,dtaROHSDIST.10000.PLINK.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.ARRAY = aggregate(dtaROHSDIST.10000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.ARRAY$SUB_TECH), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.ARRAY = aggregate(dtaROHSDIST.10000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.ARRAY$SUB_TECH, SimuIDs=dtaROHSDIST.10000.PLINK.ARRAY$SimuID, REP=dtaROHSDIST.10000.PLINK.ARRAY$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.ARRAY = merge(dtaROHSDIST.plot.10000.PLINK.ARRAY, ROHsDist.SD.10000.PLINK.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.ARRAY = dtaROHSDIST.plot.10000.PLINK.ARRAY[order(dtaROHSDIST.plot.10000.PLINK.ARRAY$CLASS, dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC),]
dtaROHSDIST.plot.10000.PLINK.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.50k);rm(dtaROHSDIST.10000.PLINK.700k);rm(dtaROHSDIST.10000.PLINK.ARRAY);rm(ROHsDist.MeanforCI.10000.PLINK.ARRAY);rm(ROHsDist.SD.10000.PLINK.ARRAY)

#Read the ne 1'000 WGS PLINK
dtaROHSDIST.10000.PLINK.WGS = read.table("./LargePop/Analyses/WGS_PLINK_ROHsDistribution.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.WGS = aggregate(dtaROHSDIST.10000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.WGS$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.WGS$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.WGS$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.WGS$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.WGS = aggregate(dtaROHSDIST.10000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.WGS$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.PLINK.WGS$SimID, REP=dtaROHSDIST.10000.PLINK.WGS$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.WGS = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.WGS$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.WGS$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.WGS$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.WGS = merge(dtaROHSDIST.plot.10000.PLINK.WGS, ROHsDist.SD.10000.PLINK.WGS, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.WGS) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.WGS = dtaROHSDIST.plot.10000.PLINK.WGS[order(dtaROHSDIST.plot.10000.PLINK.WGS$CLASS, dtaROHSDIST.plot.10000.PLINK.WGS$PERC),]
dtaROHSDIST.plot.10000.PLINK.WGS$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.WGS$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.WGS);rm(ROHsDist.MeanforCI.10000.PLINK.WGS);rm(ROHsDist.SD.10000.PLINK.WGS)

#Read the ne 1'000 TRUE IBD 100GEN PLINK
dtaROHSDIST.10000.PLINK.TRUEIBD100 = read.table("./LargePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD100$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD100$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.PLINK.TRUEIBD100$SimID, REP=dtaROHSDIST.10000.PLINK.TRUEIBD100$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.TRUEIBD100 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100 = merge(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100, ROHsDist.SD.10000.PLINK.TRUEIBD100, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100 = dtaROHSDIST.plot.10000.PLINK.TRUEIBD100[order(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$CLASS, dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC),]
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.TRUEIBD100);rm(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100);rm(ROHsDist.SD.10000.PLINK.TRUEIBD100)

#Read the ne 1'000 TRUE IBD 10000GEN PLINK
dtaROHSDIST.10000.PLINK.TRUEIBD1000 = read.table("./LargePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD1000$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.TRUEIBD1000 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD1000$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.PLINK.TRUEIBD1000$SimID, REP=dtaROHSDIST.10000.PLINK.TRUEIBD1000$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.TRUEIBD1000 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD1000$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD1000$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD1000$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000 = merge(dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000, ROHsDist.SD.10000.PLINK.TRUEIBD1000, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000 = dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000[order(dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000$CLASS, dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000$PERC),]
dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.TRUEIBD1000);rm(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD1000);rm(ROHsDist.SD.10000.PLINK.TRUEIBD1000)

#Merge ALL and remove duplicated WGS lines
dtaROHSDIST.10000.PLINK = unique(rbind(dtaROHSDIST.plot.10000.PLINK.RAD, dtaROHSDIST.plot.10000.PLINK.ARRAY, dtaROHSDIST.plot.10000.PLINK.WGS, dtaROHSDIST.plot.10000.PLINK.TRUEIBD100, dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000))
#Change factor levels
dtaROHSDIST.10000.PLINK$PERC = factor(dtaROHSDIST.10000.PLINK$PERC, levels = c("15000", "25000", "60000","SMALL_ARRAY","LARGE_ARRAY","WGS","TRUE_IBD_100", "TRUE_IBD_1000"))
#Order the df
dtaROHSDIST.10000.PLINK = dtaROHSDIST.10000.PLINK[order(dtaROHSDIST.10000.PLINK$CLASS, dtaROHSDIST.10000.PLINK$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.10000.PLINK.RAD); rm(dtaROHSDIST.plot.10000.PLINK.ARRAY,dtaROHSDIST.plot.10000.PLINK.WGS,dtaROHSDIST.plot.10000.PLINK.TRUEIBD100,dtaROHSDIST.plot.10000.PLINK.TRUEIBD1000)
#Divide TRUE IBD by 1000 (because in basepair while PLINK in kb)
dtaROHSDIST.10000.PLINK$x[dtaROHSDIST.10000.PLINK$PERC %in% c("TRUE_IBD_100", "TRUE_IBD_1000")] = dtaROHSDIST.10000.PLINK$x[dtaROHSDIST.10000.PLINK$PERC %in% c("TRUE_IBD_100", "TRUE_IBD_1000")]/1000

#### Read the data ROHsDist RZooRoH Ne 1'000 ####

dtaROHSDIST.1000.Rzoo.RAD = read.table("./SmallPop/Analyses/Radseq_RZooRoH_ROHsDistribution.txt", header = T)
#Take only the perc we are interested in
dtaROHSDIST.1000.Rzoo.RAD = dtaROHSDIST.1000.Rzoo.RAD[dtaROHSDIST.1000.Rzoo.RAD$WINDOWS.NB %in% c("3000", "9000", "120000"),]
colnames(dtaROHSDIST.1000.Rzoo.RAD) = c("SimuID","SUB_method","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.1000.Rzoo.RAD = aggregate(dtaROHSDIST.1000.Rzoo.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.Rzoo.RAD$ROHs_CLASS, PERC=dtaROHSDIST.1000.Rzoo.RAD$SUB_method), FUN = mean)
#Calculate sd
ROHsDist.SD.1000.Rzoo.RAD = as.data.frame(as.matrix(aggregate(dtaROHSDIST.1000.Rzoo.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.Rzoo.RAD$ROHs_CLASS, PERC=dtaROHSDIST.1000.Rzoo.RAD$SUB_method), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.Rzoo.RAD = merge(dtaROHSDIST.plot.1000.Rzoo.RAD, ROHsDist.SD.1000.Rzoo.RAD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.Rzoo.RAD) = c("CLASS", "PERC", "x", "sd")
#Pass sd to numeric
dtaROHSDIST.plot.1000.Rzoo.RAD$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.Rzoo.RAD$sd))
#Order the df
dtaROHSDIST.plot.1000.Rzoo.RAD = dtaROHSDIST.plot.1000.Rzoo.RAD[order(dtaROHSDIST.plot.1000.Rzoo.RAD$CLASS, dtaROHSDIST.plot.1000.Rzoo.RAD$PERC),]
dtaROHSDIST.plot.1000.Rzoo.RAD$PERC = factor(dtaROHSDIST.plot.1000.Rzoo.RAD$PERC)
#rm RAD and SD
rm(dtaROHSDIST.1000.Rzoo.RAD); rm(ROHsDist.SD.1000.Rzoo.RAD)

#Read the ne 1'000 50k array Rzoo
dtaROHSDIST.1000.Rzoo.50k = read.table("./SmallPop/Analyses/SmallArray_RZooRoH_ROHsDistributions.txt", header = T)
#Read the ne 1'000 700k array Rzoo
dtaROHSDIST.1000.Rzoo.700k = read.table("./SmallPop/Analyses/LargeArray_RZooRoH_ROHsDistributions.txt", header = T)
dtaROHSDIST.1000.Rzoo.ARRAY = unique(rbind(dtaROHSDIST.1000.Rzoo.50k,dtaROHSDIST.1000.Rzoo.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.1000.Rzoo.ARRAY = aggregate(dtaROHSDIST.1000.Rzoo.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.Rzoo.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.1000.Rzoo.ARRAY$SUB_TECH), FUN = mean)
dtaROHSDIST.plot.1000.Rzoo.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.1000.Rzoo.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.1000.Rzoo.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.Rzoo.ARRAY = aggregate(dtaROHSDIST.1000.Rzoo.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.Rzoo.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.1000.Rzoo.ARRAY$SUB_TECH, SimuIDs=dtaROHSDIST.1000.Rzoo.ARRAY$SimuID, REP=dtaROHSDIST.1000.Rzoo.ARRAY$REP), FUN = mean)
ROHsDist.SD.1000.Rzoo.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.Rzoo.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.1000.Rzoo.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.1000.Rzoo.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.Rzoo.ARRAY = merge(dtaROHSDIST.plot.1000.Rzoo.ARRAY, ROHsDist.SD.1000.Rzoo.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.Rzoo.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.Rzoo.ARRAY = dtaROHSDIST.plot.1000.Rzoo.ARRAY[order(dtaROHSDIST.plot.1000.Rzoo.ARRAY$CLASS, dtaROHSDIST.plot.1000.Rzoo.ARRAY$PERC),]
dtaROHSDIST.plot.1000.Rzoo.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.Rzoo.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.Rzoo.50k);rm(dtaROHSDIST.1000.Rzoo.700k);rm(dtaROHSDIST.1000.Rzoo.ARRAY);rm(ROHsDist.MeanforCI.1000.Rzoo.ARRAY);rm(ROHsDist.SD.1000.Rzoo.ARRAY)

#Read the ne 1'000 WGS Rzoo
dtaROHSDIST.1000.Rzoo.WGS = read.table("./SmallPop/Analyses//WGS_RZooRoH_ROHsDistribution.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.Rzoo.WGS = aggregate(dtaROHSDIST.1000.Rzoo.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.Rzoo.WGS$ROHs_CLASS, PERC=dtaROHSDIST.1000.Rzoo.WGS$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.Rzoo.WGS$PERC = factor(as.character(dtaROHSDIST.plot.1000.Rzoo.WGS$PERC), levels = unique(dtaROHSDIST.plot.1000.Rzoo.WGS$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.Rzoo.WGS = aggregate(dtaROHSDIST.1000.Rzoo.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.Rzoo.WGS$ROHs_CLASS, PERC=dtaROHSDIST.1000.Rzoo.WGS$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.Rzoo.WGS$SimID, REP=dtaROHSDIST.1000.Rzoo.WGS$REP), FUN = mean)
ROHsDist.SD.1000.Rzoo.WGS = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.Rzoo.WGS$x, by = list(CLASS=ROHsDist.MeanforCI.1000.Rzoo.WGS$CLASS, PERC=ROHsDist.MeanforCI.1000.Rzoo.WGS$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.Rzoo.WGS = merge(dtaROHSDIST.plot.1000.Rzoo.WGS, ROHsDist.SD.1000.Rzoo.WGS, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.Rzoo.WGS) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.Rzoo.WGS = dtaROHSDIST.plot.1000.Rzoo.WGS[order(dtaROHSDIST.plot.1000.Rzoo.WGS$CLASS, dtaROHSDIST.plot.1000.Rzoo.WGS$PERC),]
dtaROHSDIST.plot.1000.Rzoo.WGS$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.Rzoo.WGS$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.Rzoo.WGS);rm(ROHsDist.MeanforCI.1000.Rzoo.WGS);rm(ROHsDist.SD.1000.Rzoo.WGS)

#Read the ne 1'000 TRUE IBD 100GEN Rzoo
dtaROHSDIST.1000.Rzoo.TRUEIBD100 = read.table("./SmallPop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100 = aggregate(dtaROHSDIST.1000.Rzoo.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.Rzoo.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.1000.Rzoo.TRUEIBD100$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100$PERC = factor(as.character(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100$PERC), levels = unique(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.Rzoo.TRUEIBD100 = aggregate(dtaROHSDIST.1000.Rzoo.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.Rzoo.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.1000.Rzoo.TRUEIBD100$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.Rzoo.TRUEIBD100$SimID, REP=dtaROHSDIST.1000.Rzoo.TRUEIBD100$REP), FUN = mean)
ROHsDist.SD.1000.Rzoo.TRUEIBD100 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.Rzoo.TRUEIBD100$x, by = list(CLASS=ROHsDist.MeanforCI.1000.Rzoo.TRUEIBD100$CLASS, PERC=ROHsDist.MeanforCI.1000.Rzoo.TRUEIBD100$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100 = merge(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100, ROHsDist.SD.1000.Rzoo.TRUEIBD100, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100 = dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100[order(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100$CLASS, dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100$PERC),]
dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.Rzoo.TRUEIBD100);rm(ROHsDist.MeanforCI.1000.Rzoo.TRUEIBD100);rm(ROHsDist.SD.1000.Rzoo.TRUEIBD100)

#Read the ne 1'000 TRUE IBD 1000GEN Rzoo
dtaROHSDIST.1000.Rzoo.TRUEIBD1000 = read.table("./SmallPop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000 = aggregate(dtaROHSDIST.1000.Rzoo.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.Rzoo.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.1000.Rzoo.TRUEIBD1000$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000$PERC = factor(as.character(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000$PERC), levels = unique(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.Rzoo.TRUEIBD1000 = aggregate(dtaROHSDIST.1000.Rzoo.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.Rzoo.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.1000.Rzoo.TRUEIBD1000$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.Rzoo.TRUEIBD1000$SimID, REP=dtaROHSDIST.1000.Rzoo.TRUEIBD1000$REP), FUN = mean)
ROHsDist.SD.1000.Rzoo.TRUEIBD1000 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.Rzoo.TRUEIBD1000$x, by = list(CLASS=ROHsDist.MeanforCI.1000.Rzoo.TRUEIBD1000$CLASS, PERC=ROHsDist.MeanforCI.1000.Rzoo.TRUEIBD1000$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000 = merge(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000, ROHsDist.SD.1000.Rzoo.TRUEIBD1000, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000 = dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000[order(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000$CLASS, dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000$PERC),]
dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.Rzoo.TRUEIBD1000);rm(ROHsDist.MeanforCI.1000.Rzoo.TRUEIBD1000);rm(ROHsDist.SD.1000.Rzoo.TRUEIBD1000)

#Merge ALL and remove duplicated WGS lines
dtaROHSDIST.1000.Rzoo = unique(rbind(dtaROHSDIST.plot.1000.Rzoo.RAD, dtaROHSDIST.plot.1000.Rzoo.ARRAY, dtaROHSDIST.plot.1000.Rzoo.WGS, dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100, dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000))
#Change factor levels
dtaROHSDIST.1000.Rzoo$PERC = factor(dtaROHSDIST.1000.Rzoo$PERC, levels = c("3000", "9000", "120000","SMALL_ARRAY","LARGE_ARRAY","WGS","TRUE_IBD_100", "TRUE_IBD_1000"))
#Order the df
dtaROHSDIST.1000.Rzoo = dtaROHSDIST.1000.Rzoo[order(dtaROHSDIST.1000.Rzoo$CLASS, dtaROHSDIST.1000.Rzoo$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.1000.Rzoo.RAD); rm(dtaROHSDIST.plot.1000.Rzoo.ARRAY,dtaROHSDIST.plot.1000.Rzoo.WGS,dtaROHSDIST.plot.1000.Rzoo.TRUEIBD100,dtaROHSDIST.plot.1000.Rzoo.TRUEIBD1000)


#### Read the data ROHsDist RZooRoH Ne 10'000 ####

dtaROHSDIST.10000.Rzoo.RAD = read.table("./LargePop/Analyses/Radseq_RZooRoH_ROHsDistribution.txt", header = T)
#Take only the perc we are interested in
dtaROHSDIST.10000.Rzoo.RAD = dtaROHSDIST.10000.Rzoo.RAD[dtaROHSDIST.10000.Rzoo.RAD$WINDOWS.NB %in% c("250", "500", "7500"),]
colnames(dtaROHSDIST.10000.Rzoo.RAD) = c("SimuID","SUB_method","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.10000.Rzoo.RAD = aggregate(dtaROHSDIST.10000.Rzoo.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.Rzoo.RAD$ROHs_CLASS, PERC=dtaROHSDIST.10000.Rzoo.RAD$SUB_method), FUN = mean)
#Calculate sd
ROHsDist.SD.10000.Rzoo.RAD = as.data.frame(as.matrix(aggregate(dtaROHSDIST.10000.Rzoo.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.Rzoo.RAD$ROHs_CLASS, PERC=dtaROHSDIST.10000.Rzoo.RAD$SUB_method), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.Rzoo.RAD = merge(dtaROHSDIST.plot.10000.Rzoo.RAD, ROHsDist.SD.10000.Rzoo.RAD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.Rzoo.RAD) = c("CLASS", "PERC", "x", "sd")
#Pass sd to numeric
dtaROHSDIST.plot.10000.Rzoo.RAD$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.Rzoo.RAD$sd))
#Order the df
dtaROHSDIST.plot.10000.Rzoo.RAD = dtaROHSDIST.plot.10000.Rzoo.RAD[order(dtaROHSDIST.plot.10000.Rzoo.RAD$CLASS, dtaROHSDIST.plot.10000.Rzoo.RAD$PERC),]
dtaROHSDIST.plot.10000.Rzoo.RAD$PERC = factor(dtaROHSDIST.plot.10000.Rzoo.RAD$PERC)
#rm RAD and SD
rm(dtaROHSDIST.10000.Rzoo.RAD); rm(ROHsDist.SD.10000.Rzoo.RAD)

#Read the ne 10'000 50k array Rzoo
dtaROHSDIST.10000.Rzoo.50k = read.table("./LargePop//Analyses/SmallArray_RZooRoH_ROHsDistributions.txt", header = T)
#Read the ne 1'000 700k array Rzoo
dtaROHSDIST.10000.Rzoo.700k = read.table("./LargePop//Analyses/LargeArray_RZooRoH_ROHsDistributions.txt", header = T)
dtaROHSDIST.10000.Rzoo.ARRAY = unique(rbind(dtaROHSDIST.10000.Rzoo.50k,dtaROHSDIST.10000.Rzoo.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.10000.Rzoo.ARRAY = aggregate(dtaROHSDIST.10000.Rzoo.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.Rzoo.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.10000.Rzoo.ARRAY$SUB_TECH), FUN = mean)
dtaROHSDIST.plot.10000.Rzoo.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.10000.Rzoo.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.10000.Rzoo.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.Rzoo.ARRAY = aggregate(dtaROHSDIST.10000.Rzoo.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.Rzoo.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.10000.Rzoo.ARRAY$SUB_TECH, SimuIDs=dtaROHSDIST.10000.Rzoo.ARRAY$SimuID, REP=dtaROHSDIST.10000.Rzoo.ARRAY$REP), FUN = mean)
ROHsDist.SD.10000.Rzoo.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.Rzoo.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.10000.Rzoo.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.10000.Rzoo.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.Rzoo.ARRAY = merge(dtaROHSDIST.plot.10000.Rzoo.ARRAY, ROHsDist.SD.10000.Rzoo.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.Rzoo.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.Rzoo.ARRAY = dtaROHSDIST.plot.10000.Rzoo.ARRAY[order(dtaROHSDIST.plot.10000.Rzoo.ARRAY$CLASS, dtaROHSDIST.plot.10000.Rzoo.ARRAY$PERC),]
dtaROHSDIST.plot.10000.Rzoo.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.Rzoo.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.Rzoo.50k);rm(dtaROHSDIST.10000.Rzoo.700k);rm(dtaROHSDIST.10000.Rzoo.ARRAY);rm(ROHsDist.MeanforCI.10000.Rzoo.ARRAY);rm(ROHsDist.SD.10000.Rzoo.ARRAY)

#Read the ne 10'000 WGS Rzoo
dtaROHSDIST.10000.Rzoo.WGS = read.table("./LargePop//Analyses//WGS_RZooRoH_ROHsDistribution.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.Rzoo.WGS = aggregate(dtaROHSDIST.10000.Rzoo.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.Rzoo.WGS$ROHs_CLASS, PERC=dtaROHSDIST.10000.Rzoo.WGS$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.Rzoo.WGS$PERC = factor(as.character(dtaROHSDIST.plot.10000.Rzoo.WGS$PERC), levels = unique(dtaROHSDIST.plot.10000.Rzoo.WGS$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.Rzoo.WGS = aggregate(dtaROHSDIST.10000.Rzoo.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.Rzoo.WGS$ROHs_CLASS, PERC=dtaROHSDIST.10000.Rzoo.WGS$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.Rzoo.WGS$SimID, REP=dtaROHSDIST.10000.Rzoo.WGS$REP), FUN = mean)
ROHsDist.SD.10000.Rzoo.WGS = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.Rzoo.WGS$x, by = list(CLASS=ROHsDist.MeanforCI.10000.Rzoo.WGS$CLASS, PERC=ROHsDist.MeanforCI.10000.Rzoo.WGS$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.Rzoo.WGS = merge(dtaROHSDIST.plot.10000.Rzoo.WGS, ROHsDist.SD.10000.Rzoo.WGS, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.Rzoo.WGS) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.Rzoo.WGS = dtaROHSDIST.plot.10000.Rzoo.WGS[order(dtaROHSDIST.plot.10000.Rzoo.WGS$CLASS, dtaROHSDIST.plot.10000.Rzoo.WGS$PERC),]
dtaROHSDIST.plot.10000.Rzoo.WGS$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.Rzoo.WGS$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.Rzoo.WGS);rm(ROHsDist.MeanforCI.10000.Rzoo.WGS);rm(ROHsDist.SD.10000.Rzoo.WGS)

#Read the ne 10'000 TRUE IBD 100GEN Rzoo
dtaROHSDIST.10000.Rzoo.TRUEIBD100 = read.table("./LargePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100 = aggregate(dtaROHSDIST.10000.Rzoo.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.Rzoo.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.10000.Rzoo.TRUEIBD100$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100$PERC = factor(as.character(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100$PERC), levels = unique(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.Rzoo.TRUEIBD100 = aggregate(dtaROHSDIST.10000.Rzoo.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.Rzoo.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.10000.Rzoo.TRUEIBD100$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.Rzoo.TRUEIBD100$SimID, REP=dtaROHSDIST.10000.Rzoo.TRUEIBD100$REP), FUN = mean)
ROHsDist.SD.10000.Rzoo.TRUEIBD100 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.Rzoo.TRUEIBD100$x, by = list(CLASS=ROHsDist.MeanforCI.10000.Rzoo.TRUEIBD100$CLASS, PERC=ROHsDist.MeanforCI.10000.Rzoo.TRUEIBD100$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100 = merge(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100, ROHsDist.SD.10000.Rzoo.TRUEIBD100, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100 = dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100[order(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100$CLASS, dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100$PERC),]
dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.Rzoo.TRUEIBD100);rm(ROHsDist.MeanforCI.10000.Rzoo.TRUEIBD100);rm(ROHsDist.SD.10000.Rzoo.TRUEIBD100)

#Read the ne 10'000 TRUE IBD 1000GEN Rzoo
dtaROHSDIST.10000.Rzoo.TRUEIBD1000 = read.table("./LargePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000 = aggregate(dtaROHSDIST.10000.Rzoo.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.Rzoo.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.10000.Rzoo.TRUEIBD1000$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000$PERC = factor(as.character(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000$PERC), levels = unique(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.Rzoo.TRUEIBD1000 = aggregate(dtaROHSDIST.10000.Rzoo.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.Rzoo.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.10000.Rzoo.TRUEIBD1000$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.Rzoo.TRUEIBD1000$SimID, REP=dtaROHSDIST.10000.Rzoo.TRUEIBD1000$REP), FUN = mean)
ROHsDist.SD.10000.Rzoo.TRUEIBD1000 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.Rzoo.TRUEIBD1000$x, by = list(CLASS=ROHsDist.MeanforCI.10000.Rzoo.TRUEIBD1000$CLASS, PERC=ROHsDist.MeanforCI.10000.Rzoo.TRUEIBD1000$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000 = merge(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000, ROHsDist.SD.10000.Rzoo.TRUEIBD1000, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000 = dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000[order(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000$CLASS, dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000$PERC),]
dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.Rzoo.TRUEIBD1000);rm(ROHsDist.MeanforCI.10000.Rzoo.TRUEIBD1000);rm(ROHsDist.SD.10000.Rzoo.TRUEIBD1000)

#Merge ALL and remove duplicated WGS lines
dtaROHSDIST.10000.Rzoo = unique(rbind(dtaROHSDIST.plot.10000.Rzoo.RAD, dtaROHSDIST.plot.10000.Rzoo.ARRAY, dtaROHSDIST.plot.10000.Rzoo.WGS, dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100, dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000))
#Change factor levels
dtaROHSDIST.10000.Rzoo$PERC = factor(dtaROHSDIST.10000.Rzoo$PERC, levels = c("250", "500", "7500","SMALL_ARRAY","LARGE_ARRAY","WGS","TRUE_IBD_100", "TRUE_IBD_1000"))
#Order the df
dtaROHSDIST.10000.Rzoo = dtaROHSDIST.10000.Rzoo[order(dtaROHSDIST.10000.Rzoo$CLASS, dtaROHSDIST.10000.Rzoo$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.10000.Rzoo.RAD); rm(dtaROHSDIST.plot.10000.Rzoo.ARRAY,dtaROHSDIST.plot.10000.Rzoo.WGS,dtaROHSDIST.plot.10000.Rzoo.TRUEIBD100,dtaROHSDIST.plot.10000.Rzoo.TRUEIBD1000)

#Define ROHs classes
classesROH = c("< 2", "2 - 4","4 - 6","6 - 10","10 - 16", "> 16")

pdf("~/PhD/ROH/Which_data/Manuscript/SUBMITTED_VERSION/Review_round_I/FIGURES/SUPPLEMENTARY_FIGURES/FIGURES3.pdf", width = 55, height = 28)

#Layout for having all plots
layout(matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5)), 5, 5, byrow = T), heights = c(.5,1,.8,1,.5), widths = c(.5,2,1,2,1))

#Par for margin and outer margings
par(mar = c(3.5,3.5,3.5,3.5), oma = c(.2,10,.2,.2))

#### PANNEL A: ROHsDist Ne 1'000 PLINK 100KB RAD + ARRAYS ####

#Define colors for subsampled dataset
subCOLORS = c("#82EEFD", "#0492C2", "#016064", "#2832C2", "#0A1172", "grey18")

#Create empty plot
plot(dtaROHSDIST.1000.PLINK$CLASS, dtaROHSDIST.1000.PLINK$x/1000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0,700))

#Loop through classes to draw boxplots
for(class in sort(unique(dtaROHSDIST.1000.PLINK$CLASS))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.36, class-0.24, class-0.12, class, class+0.12, class+0.24,class+0.36)
  
  #Loop through reduced representation
  for(reduced in 1:6){
    
    #subset the data
    dta_sub = dtaROHSDIST.1000.PLINK[dtaROHSDIST.1000.PLINK$CLASS == class & dtaROHSDIST.1000.PLINK$PERC == levels(dtaROHSDIST.1000.PLINK$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$CLASS == class & dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_1000"],
         ytop = dta_sub$x/1000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.06), x1 = (polyXcoord[reduced] + 0.06), y0 = ((dta_sub$x - dta_sub$sd)/1000), y1 = ((dta_sub$x + dta_sub$sd)/1000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_1000" & dtaROHSDIST.1000.PLINK$CLASS == class],
             y1 = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_1000" & dtaROHSDIST.1000.PLINK$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add the y axis
axis(2, cex.axis = 4, hadj = 1.5, las = 1, cex.axis = 5, lwd.ticks = 5, tck = -0.03)
#Add y axis title
mtext(text = "Mean Mb ± SD (among individuals)", side = 2, cex = 3.5, line = 20)
#Add x axis title
mtext(text = "ROHs Length Classes [Mb]", side = 1, cex = 4, line = 15)
#Add panel lab
mtext(text = "A", side = 3, cex = 7, line = 15, at = -2)
mtext(text = "PLINK", side = 3, cex = 5, line = 15, at = 3)

#### PANNEL B: ROHsDist Ne 1'000 RZooRoH RAD + ARRAYS ####

#Define colors for subsampled dataset
subCOLORS = c("#82EEFD", "#0492C2", "#016064", "#2832C2", "#0A1172", "grey18")

#Create empty plot
plot(dtaROHSDIST.1000.Rzoo$CLASS, dtaROHSDIST.1000.Rzoo$x/1000000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0,700))

#Loop through classes to draw boxplots
for(class in sort(unique(dtaROHSDIST.1000.Rzoo$CLASS))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.36, class-0.24, class-0.12, class, class+0.12, class+0.24,class+0.36)
  
  #Loop through reduced representation
  for(reduced in 1:6){
    
    #subset the data
    dta_sub = dtaROHSDIST.1000.Rzoo[dtaROHSDIST.1000.Rzoo$CLASS == class & dtaROHSDIST.1000.Rzoo$PERC == levels(dtaROHSDIST.1000.Rzoo$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.1000.Rzoo$x/1000000)[dtaROHSDIST.1000.Rzoo$CLASS == class & dtaROHSDIST.1000.Rzoo$PERC == "TRUE_IBD_1000"],
         ytop = dta_sub$x/1000000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.06), x1 = (polyXcoord[reduced] + 0.06), y0 = ((dta_sub$x - dta_sub$sd)/1000000), y1 = ((dta_sub$x + dta_sub$sd)/1000000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.1000.Rzoo$x/1000000)[dtaROHSDIST.1000.Rzoo$PERC == "TRUE_IBD_1000" & dtaROHSDIST.1000.Rzoo$CLASS == class],
             y1 = (dtaROHSDIST.1000.Rzoo$x/1000000)[dtaROHSDIST.1000.Rzoo$PERC == "TRUE_IBD_1000" & dtaROHSDIST.1000.Rzoo$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add the y axis
axis(2, cex.axis = 4, hadj = 1.5, las = 1, cex.axis = 5, lwd.ticks = 5, tck = -0.03)
#Add y axis title
mtext(text = "Mean Mb ± SD (among individuals)", side = 2, cex = 3.5, line = 20)
#Add x axis title
mtext(text = "HBD Length Classes [Mb]", side = 1, cex = 4, line = 15)
#Add panel lab
mtext(text = "B", side = 3, cex = 7, line = 15, at = -2)
mtext(text = "RZooRoH", side = 3, cex = 5, line = 15, at = 3)

#### PANNEL C: ROHsDist Ne 10'000 PLINK 100KB RAD + ARRAYS ####

subCOLORS = c("#FFBF00", "#FD6A02", "#FDA172","#CC7722", "#793802", "grey18")

#Create empty plot
plot(dtaROHSDIST.10000.PLINK$CLASS, dtaROHSDIST.10000.PLINK$x/1000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0,600))

#Loop through classes to draw boxplots
for(class in sort(unique(dtaROHSDIST.10000.PLINK$CLASS))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.36, class-0.24, class-0.12, class, class+0.12, class+0.24,class+0.36)
  
  #Loop through reduced representation
  for(reduced in 1:6){
    
    #subset the data
    dta_sub = dtaROHSDIST.10000.PLINK[dtaROHSDIST.10000.PLINK$CLASS == class & dtaROHSDIST.10000.PLINK$PERC == levels(dtaROHSDIST.10000.PLINK$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$CLASS == class & dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_1000"],
         ytop = dta_sub$x/1000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.06), x1 = (polyXcoord[reduced] + 0.06), y0 = ((dta_sub$x - dta_sub$sd)/1000), y1 = ((dta_sub$x + dta_sub$sd)/1000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_1000" & dtaROHSDIST.10000.PLINK$CLASS == class],
             y1 = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_1000" & dtaROHSDIST.10000.PLINK$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add the y axis
axis(2, cex.axis = 4, hadj = 1.5, las = 1, cex.axis = 5, lwd.ticks = 5, tck = -0.03)
#Add y axis title
mtext(text = "Mean Mb ± SD (among individuals)", side = 2, cex = 3.5, line = 20)
#Add x axis title
mtext(text = "ROHs Length Classes [Mb]", side = 1, cex = 4, line = 15)
#Add panel lab
mtext(text = "C", side = 3, cex = 7, line = 15, at = -2)
mtext(text = "PLINK", side = 3, cex = 5, line = 15, at = 3)


#### PANNEL D: ROHsDist Ne 10'000 RZooRoH RAD + ARRAYS ####

subCOLORS = c("#FFBF00", "#FD6A02", "#FDA172","#CC7722", "#793802", "grey18")

#Create empty plot
plot(dtaROHSDIST.10000.Rzoo$CLASS, dtaROHSDIST.10000.Rzoo$x/1000000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0,600))

#Loop through classes to draw boxplots
for(class in sort(unique(dtaROHSDIST.10000.Rzoo$CLASS))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.36, class-0.24, class-0.12, class, class+0.12, class+0.24,class+0.36)
  
  #Loop through reduced representation
  for(reduced in 1:6){
    
    #subset the data
    dta_sub = dtaROHSDIST.10000.Rzoo[dtaROHSDIST.10000.Rzoo$CLASS == class & dtaROHSDIST.10000.Rzoo$PERC == levels(dtaROHSDIST.10000.Rzoo$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.10000.Rzoo$x/1000000)[dtaROHSDIST.10000.Rzoo$CLASS == class & dtaROHSDIST.10000.Rzoo$PERC == "TRUE_IBD_1000"],
         ytop = dta_sub$x/1000000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.06), x1 = (polyXcoord[reduced] + 0.06), y0 = ((dta_sub$x - dta_sub$sd)/1000000), y1 = ((dta_sub$x + dta_sub$sd)/1000000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.10000.Rzoo$x/1000000)[dtaROHSDIST.10000.Rzoo$PERC == "TRUE_IBD_1000" & dtaROHSDIST.10000.Rzoo$CLASS == class],
             y1 = (dtaROHSDIST.10000.Rzoo$x/1000000)[dtaROHSDIST.10000.Rzoo$PERC == "TRUE_IBD_1000" & dtaROHSDIST.10000.Rzoo$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add the y axis
axis(2, cex.axis = 4, hadj = 1.5, las = 1, cex.axis = 5, lwd.ticks = 5, tck = -0.03)
#Add y axis title
mtext(text = "Mean Mb ± SD (among individuals)", side = 2, cex = 3.5, line = 20)
#Add x axis title
mtext(text = "HBD Length Classes [Mb]", side = 1, cex = 4, line = 15)
#Add panel lab
mtext(text = "D", side = 3, cex = 7, line = 15, at = -2)
mtext(text = "RZooRoH", side = 3, cex = 5, line = 15, at = 3)

#Close the plot
dev.off()

######################################################################################################
################## FIGURE S4 PLINK DEFAULT PARAMETERS CORRELATIONS FROH ~SNP DENSITY #################
######################################################################################################

#Read the data ne 1'000 100KB
dtaCOR.1000.PLINK = read.table("./SmallPop/Analyses/CORRELATION_according_to_DENSITY_PLINK_DefaultPARAM.txt", header = T)
#Read the data ne 10'000 100KB
dtaCOR.10000.PLINK = read.table("./LargePop/Analyses/CORRELATION_according_to_DENSITY_PLINK_DefaultPARAM.txt", header = T)

##### JPEG ####

jpeg("~/PhD/ROH/Which_data/Manuscript/SUBMITTED_VERSION/Review_round_I/FIGURES/SUPPLEMENTARY_FIGURES/FIGURES4.jpg", width = 4500, height = 2100)

#Layout for having all plots
layout(matrix(c(rep(0,5),0,1,0,2,0,rep(0,5)), 3, 5, byrow = T), heights = c(.2,1,.2), widths = c(.5,1,1,1,.5))

#Par for margin and outer margings
par(mar = c(3.5,3.5,3.5,3.5), oma = c(.2,.2,.2,.2))

#### PANNEL A: CORR ~ DENSITY PLINK SmallNe ####

## CORRELATION
plot(dtaCOR.1000.PLINK$CorrelarionFsubFwgs ~ dtaCOR.1000.PLINK$NSNPsPerMb,  pch = 20, yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA, col = add.alpha("#82EEFD",0.8), bg = add.alpha("#82EEFD",0.1), lwd = 0, ylim = c(-1,1))
#Add thre x axis
axis(1, cex.axis = 8, padj = 1.5, lwd.ticks = 3)
#Add thre y axis
axis(2, cex.axis = 8, hadj = 1.5, las = 1, lwd.ticks = 3)

#Add x and Y axes names
#mtext(text = "Mean SNPs Density [SNPs/Mb]", side = 1, cex = 4, line = 15)
mtext(text = expression(r^{2} ~ " (" ~ F["IBD"] ~ ", " ~ F["ROH RAD"] ~ ")"), side = 2, cex = 6, line = 18)
mtext(text = "A", side = 3, cex = 10, line = 18, at = 1)
mtext(text = "PLINK", side = 3, cex = 8, line = 18, at = 15)
mtext(text = "Mean SNPs Density [SNPs/Mb]", side = 1, cex = 4, line = 20)

#### PANNEL B: CORR ~ DENSITY PLINK ####

## CORRELATION
plot(dtaCOR.10000.PLINK$CorrelarionFsubFwgs ~ dtaCOR.10000.PLINK$NSNPsPerMb,  pch = 20, yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA, col = add.alpha("#FCAE1E",0.8), bg = add.alpha("#FCAE1E",0.1), lwd = 0, ylim = c(-1,1))
#Add thre x axis
axis(1, cex.axis = 8, padj = 1.5, lwd.ticks = 3)
#Add thre y axis
axis(2, cex.axis = 8, hadj = 1.5, las = 1, lwd.ticks = 3)

#Add x and Y axes names
#mtext(text = "Mean SNPs Density [SNPs/Mb]", side = 1, cex = 4, line = 15)
mtext(text = expression(r^{2} ~ " (" ~ F["IBD"] ~ ", " ~ F["ROH RAD"] ~ ")"), side = 2, cex = 6, line = 18)
mtext(text = "B", side = 3, cex = 10, line = 18, at = 1)
mtext(text = "PLINK", side = 3, cex = 8, line = 18, at = 15)
mtext(text = "Mean SNPs Density [SNPs/Mb]", side = 1, cex = 4, line = 20)

dev.off()

########################################################################################
################## FIGURE S5 PLINK DEFAULT PARAMETERS FSUB SUB VS FIBD #################
########################################################################################

#Read the data FROH ne 1'000
dtaFROH.1000.100KB.RAD = read.table("./SmallPop/Analyses/RADSEQ_PLINK_defaultPARAM_FROH.txt", header = T)
#Subsample the 4 PERCENTAGES
dtaFROH.1000.100KB.RAD = dtaFROH.1000.100KB.RAD[dtaFROH.1000.100KB.RAD$RAD_WINDOW_NB %in% c(180000,240000,600000),]
#read the WGS data to merge wwith TRUE IBD
dtaWGS = read.table("./SmallPop/Analyses/WGS_PLINK_FROH.txt", header = T)[,c(1,2,6)]
#change indv to match
dtaWGS$individuals = paste0("tsk_", dtaWGS$individuals)
colnames(dtaWGS)[1] = "Simu_ID"
#merge with sub data default param
dtaFROH.1000.100KB.RAD = merge(dtaFROH.1000.100KB.RAD,dtaWGS, by = c("Simu_ID","individuals"))

#Read the data FROH ne 10'000
dtaFROH.10000.100KB.RAD = read.table("./LargePop/Analyses/RADSEQ_PLINK_defaultPARAM_FROH.txt", header = T)
#Subsample the 4 PERCENTAGEs
dtaFROH.10000.100KB.RAD = dtaFROH.10000.100KB.RAD[dtaFROH.10000.100KB.RAD$RAD_WINDOW_NB %in% c(15000, 25000, 60000),]
#read the WGS data to merge wwith TRUE IBD
dtaWGS = read.table("./LargePop/Analyses/WGS_PLINK_FROH.txt", header = T)[,c(1,2,6)]
#change indv to match
dtaWGS$individuals = paste0("tsk_", dtaWGS$individuals)
colnames(dtaWGS)[1] = "Simu_ID"
#merge with sub data default param
dtaFROH.10000.100KB.RAD = merge(dtaFROH.10000.100KB.RAD,dtaWGS, by = c("Simu_ID","individuals"))

##Read the data Ne 1'000 50K
dtaFROH.1000.100KB.50K = read.table("./SmallPop/Analyses/SmallArray_PLINK_defaultPARAM_FROH.txt", header = T)
dtaFROH.1000.100KB.50K = cbind(dtaFROH.1000.100KB.50K, ARRAY=rep("SMALL", nrow(dtaFROH.1000.100KB.50K)))
colnames(dtaFROH.1000.100KB.50K) = c("Simu_ID", "individuals", "Froh_array", "NROH_array", "SROH_array", "Replicate", "Froh_WGS", "NROH_WGS", "SROH_WGS", "ARRAY")
#Read the data Ne 1'000 700K
dtaFROH.1000.100KB.700K = read.table("./SmallPop/Analyses/LargeArray_PLINK_defaultPARAM_FROH.txt", header = T)
dtaFROH.1000.100KB.700K = cbind(dtaFROH.1000.100KB.700K, ARRAY=rep("LARGE", nrow(dtaFROH.1000.100KB.700K)))
colnames(dtaFROH.1000.100KB.700K) = c("Simu_ID", "individuals", "Froh_array", "NROH_array", "SROH_array", "Replicate", "Froh_WGS", "NROH_WGS", "SROH_WGS", "ARRAY")
dtaFROH.1000.100KB.ARRAY = rbind(dtaFROH.1000.100KB.50K,dtaFROH.1000.100KB.700K)
rm(dtaFROH.1000.100KB.50K);rm(dtaFROH.1000.100KB.700K)
#read the WGS data to merge wwith TRUE IBD
dtaWGS = read.table("./SmallPop/Analyses/WGS_PLINK_FROH.txt", header = T)[,c(1,2,6)]
#change indv to match
dtaWGS$individuals = paste0("tsk_", dtaWGS$individuals)
colnames(dtaWGS)[1] = "Simu_ID"
#merge with sub data default param
dtaFROH.1000.100KB.ARRAY = merge(dtaFROH.1000.100KB.ARRAY,dtaWGS, by = c("Simu_ID","individuals"))
dtaFROH.1000.100KB.ARRAY$ARRAY = factor(dtaFROH.1000.100KB.ARRAY$ARRAY, levels = c("SMALL","LARGE"))

#Read the data Ne 10'000 50K
dtaFROH.10000.100KB.50K = read.table("./LargePop/Analyses/SmallArray_PLINK_defaultPARAM_FROH.txt", header = T)
dtaFROH.10000.100KB.50K = cbind(dtaFROH.10000.100KB.50K, ARRAY=rep("SMALL", nrow(dtaFROH.10000.100KB.50K)))
colnames(dtaFROH.10000.100KB.50K) = c("Simu_ID", "individuals", "Froh_array", "NROH_array", "SROH_array", "Replicate", "Froh_WGS", "NROH_WGS", "SROH_WGS", "ARRAY")
#Read the data Ne 10'000 700K
dtaFROH.10000.100KB.700K = read.table("./LargePop/Analyses/LargeArray_PLINK_defaultPARAM_FROH.txt", header = T)
dtaFROH.10000.100KB.700K = cbind(dtaFROH.10000.100KB.700K, ARRAY=rep("LARGE", nrow(dtaFROH.10000.100KB.700K)))
colnames(dtaFROH.10000.100KB.700K) = c("Simu_ID", "individuals", "Froh_array", "NROH_array", "SROH_array", "Replicate", "Froh_WGS", "NROH_WGS", "SROH_WGS", "ARRAY")
dtaFROH.10000.100KB.ARRAY = rbind(dtaFROH.10000.100KB.50K,dtaFROH.10000.100KB.700K)
rm(dtaFROH.10000.100KB.50K);rm(dtaFROH.10000.100KB.700K)
#read the WGS data to merge wwith TRUE IBD
dtaWGS = read.table("./LargePop/Analyses/WGS_PLINK_FROH.txt", header = T)[,c(1,2,6)]
#change indv to match
dtaWGS$individuals = paste0("tsk_", dtaWGS$individuals)
colnames(dtaWGS)[1] = "Simu_ID"
#merge with sub data default param
dtaFROH.10000.100KB.ARRAY = merge(dtaFROH.10000.100KB.ARRAY,dtaWGS, by = c("Simu_ID","individuals"))
dtaFROH.10000.100KB.ARRAY$ARRAY = factor(dtaFROH.10000.100KB.ARRAY$ARRAY, levels = c("SMALL","LARGE"))

#PLOT
pdf("~/PhD/ROH/Which_data/Manuscript/SUBMITTED_VERSION/Review_round_I/FIGURES/SUPPLEMENTARY_FIGURES/FIGURES5.pdf", width = 40, height = 20)

#Layout for having all plots
layout(matrix(c(rep(0,5),0,1,0,2,0,rep(0,5)), 3, 5, byrow = T), heights = c(.5,1,.5), widths = c(.5,1,1,1,.5))

#Par for margin and outer margings
par(mar = c(3.5,3.5,3.5,3.5), oma = c(.2,.2,.2,.2))

#### PANEL A PLINK SMALL NE ####

#Plot PLINK RADs SMALL POP EMPTY PLOT
plot(dtaFROH.1000.100KB.RAD$Froh_sub_classical ~ dtaFROH.1000.100KB.RAD$Froh_TRUE_IBD_100GEN, pch = 19,
     ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(0,0.9), ylim = c(0,0.9), cex.axis = 4, type = "n", frame = F)

#Colors RAD SMALL NE
RADcolours = c("#82EEFD", "#0492C2", "#016064")
names(RADcolours) = sort(unique(dtaFROH.1000.100KB.RAD$RAD_WINDOW_NB))
#Shape points SMALL NE
POINTSshapes = c(15,16,17)
names(POINTSshapes) = sort(unique(dtaFROH.1000.100KB.RAD$RAD_WINDOW_NB))

#Loop through the RAD subset SMALL NE
for (sub in sort(unique(dtaFROH.1000.100KB.RAD$RAD_WINDOW_NB))) {
  
  #Add the points SMALL NE
  points(dtaFROH.1000.100KB.RAD$Froh_sub_classical[dtaFROH.1000.100KB.RAD$RAD_WINDOW_NB == sub] ~
           dtaFROH.1000.100KB.RAD$Froh_TRUE_IBD_100GEN[dtaFROH.1000.100KB.RAD$RAD_WINDOW_NB == sub],
         col = add.alpha(RADcolours[names(RADcolours) == sub],1), pch = POINTSshapes[names(POINTSshapes) == sub], cex = 3)
  
}

#Add the points ARRAYs
points(dtaFROH.1000.100KB.ARRAY$Froh_array ~ dtaFROH.1000.100KB.ARRAY$Froh_TRUE_IBD_100GEN, col = add.alpha(c("#2832C2", "#0A1172")[dtaFROH.1000.100KB.ARRAY$ARRAY],.8),
       pch = c(15,16)[dtaFROH.1000.100KB.ARRAY$ARRAY], cex = 3)

## ADD AXES
axis(1,xlab = NULL, at=c(0:9)*0.1, padj = 1.5, cex.axis= 8, lwd.ticks = 3, srt = 45, tck = -0.03)
axis(2,xlab = NULL, at=c(0:9)*0.1, hadj = 1.5, cex.axis= 8, lwd.ticks = 3, las = 1, tck = -0.03)

#Add axes titles
mtext(expression(F["IBD"]), side = 1, outer = FALSE, line = 20, cex = 6)
mtext(expression(F["ROH"]), side = 2, outer = FALSE, line = 20, cex = 6)

#Add perfect line (must be below everything else)
abline(0,1, lwd = 4)

#Add panel legend
mtext("A", side = 3, outer = FALSE, line = 20, cex = 10, at = -.35)
#Add panel lab
mtext("PLINK", side = 3, outer = FALSE, line = 20, cex = 8, at = .4)


#### PANEL B PLINK LARGE NE ####

#Plot PLINK RADs SMALL POP EMPTY PLOT
plot(dtaFROH.10000.100KB.RAD$FROH_sub ~ dtaFROH.10000.100KB.RAD$Froh_TRUE_IBD_100GEN, pch = 19,
     ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(0,0.9), ylim = c(0,0.9), cex.axis = 4, type = "n", frame = F)

#Colors RAD SMALL NE
RADcolours = c("#FFBF00", "#FD6A02", "#FDA172")
names(RADcolours) = sort(unique(dtaFROH.10000.100KB.RAD$RAD_WINDOW_NB))
#Shape points SMALL NE
POINTSshapes = c(1,4,8)
names(POINTSshapes) = sort(unique(dtaFROH.10000.100KB.RAD$RAD_WINDOW_NB))

#Loop through the RAD subset SMALL NE
for (sub in sort(unique(dtaFROH.10000.100KB.RAD$RAD_WINDOW_NB))) {
  
  #Add the points SMALL NE
  points(dtaFROH.10000.100KB.RAD$FROH_sub[dtaFROH.10000.100KB.RAD$RAD_WINDOW_NB == sub] ~
           dtaFROH.10000.100KB.RAD$Froh_TRUE_IBD_100GEN[dtaFROH.10000.100KB.RAD$RAD_WINDOW_NB == sub],
         col = add.alpha(RADcolours[names(RADcolours) == sub],1), pch = POINTSshapes[names(POINTSshapes) == sub], cex = 3)
  
}

#Add the points ARRAYs
points(dtaFROH.10000.100KB.ARRAY$Froh_array ~ dtaFROH.10000.100KB.ARRAY$Froh_TRUE_IBD_100GEN, col = add.alpha(c("#CC7722", "#793802")[dtaFROH.10000.100KB.ARRAY$ARRAY],.8),
       pch = c(15,16)[dtaFROH.10000.100KB.ARRAY$ARRAY], cex = 3)

## ADD AXES
axis(1,xlab = NULL, at=c(0:9)*0.1, padj = 1.5, cex.axis= 8, lwd.ticks = 3, srt = 45, tck = -0.03)
axis(2,xlab = NULL, at=c(0:9)*0.1, hadj = 1.5, cex.axis= 8, lwd.ticks = 3, las = 1, tck = -0.03)

#Add axes titles
mtext(expression(F["IBD"]), side = 1, outer = FALSE, line = 20, cex = 6)
mtext(expression(F["ROH"]), side = 2, outer = FALSE, line = 20, cex = 6)

#Add perfect line (must be below everything else)
abline(0,1, lwd = 4)

#Add panel legend
mtext("B", side = 3, outer = FALSE, line = 20, cex = 10, at = -.5)
#Add panel lab
mtext("PLINK", side = 3, outer = FALSE, line = 20, cex = 8, at = .4)

dev.off()

#################################################################################
################## FIGURE S6 PLINK DEFAULT PARAMETERS ROHs DIST #################
#################################################################################

#### Read the data ROHsDist PLINK100KB Ne 1'000 ####

dtaROHSDIST.1000.PLINK.RAD = read.table("./SmallPop/Analyses/RADSEQ_PLINK_defaultPARAM_ROHsDistribution.txt", header = T)[,-2]
#Take only the perc we are interested in
dtaROHSDIST.1000.PLINK.RAD = dtaROHSDIST.1000.PLINK.RAD[dtaROHSDIST.1000.PLINK.RAD$NB_WIND %in% c("180000", "240000", "600000"),]
colnames(dtaROHSDIST.1000.PLINK.RAD) = c("SimuID","SUB_method","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.RAD = aggregate(dtaROHSDIST.1000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.RAD$SUB_method), FUN = mean)
#Calculate sd
ROHsDist.SD.1000.PLINK.RAD = as.data.frame(as.matrix(aggregate(dtaROHSDIST.1000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.RAD$SUB_method), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.RAD = merge(dtaROHSDIST.plot.1000.PLINK.RAD, ROHsDist.SD.1000.PLINK.RAD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.RAD) = c("CLASS", "PERC", "x", "sd")
#Pass sd to numeric
dtaROHSDIST.plot.1000.PLINK.RAD$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.RAD$sd))
#Order the df
dtaROHSDIST.plot.1000.PLINK.RAD = dtaROHSDIST.plot.1000.PLINK.RAD[order(dtaROHSDIST.plot.1000.PLINK.RAD$CLASS, dtaROHSDIST.plot.1000.PLINK.RAD$PERC),]
dtaROHSDIST.plot.1000.PLINK.RAD$PERC = factor(dtaROHSDIST.plot.1000.PLINK.RAD$PERC)
#rm RAD and SD
rm(dtaROHSDIST.1000.PLINK.RAD); rm(ROHsDist.SD.1000.PLINK.RAD)

#Read the ne 1'000 50k array PLINK
dtaROHSDIST.1000.PLINK.50k = read.table("./SmallPop/Analyses/SmallArray_PLINK_defaultPARAM_ROHsDistribution.txt", header = T)
dtaROHSDIST.1000.PLINK.50k$SUB_method = "SMALL_ARRAY"
#Read the ne 1'000 700k array PLINK
dtaROHSDIST.1000.PLINK.700k = read.table("./SmallPop/Analyses/LargeArray_PLINK_defaultPARAM_ROHsDistribution.txt", header = T)
dtaROHSDIST.1000.PLINK.700k$SUB_method = "LARGE_ARRAY"
dtaROHSDIST.1000.PLINK.ARRAY = unique(rbind(dtaROHSDIST.1000.PLINK.50k,dtaROHSDIST.1000.PLINK.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.ARRAY = aggregate(dtaROHSDIST.1000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.ARRAY$SUB_method), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.ARRAY = aggregate(dtaROHSDIST.1000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.ARRAY$SUB_method, SimuIDs=dtaROHSDIST.1000.PLINK.ARRAY$SimuID, REP=dtaROHSDIST.1000.PLINK.ARRAY$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.ARRAY = merge(dtaROHSDIST.plot.1000.PLINK.ARRAY, ROHsDist.SD.1000.PLINK.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.ARRAY = dtaROHSDIST.plot.1000.PLINK.ARRAY[order(dtaROHSDIST.plot.1000.PLINK.ARRAY$CLASS, dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC),]
dtaROHSDIST.plot.1000.PLINK.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.50k);rm(dtaROHSDIST.1000.PLINK.700k);rm(dtaROHSDIST.1000.PLINK.ARRAY);rm(ROHsDist.MeanforCI.1000.PLINK.ARRAY);rm(ROHsDist.SD.1000.PLINK.ARRAY)

#Read the ne 1'000 TRUE IBD PLINK
dtaROHSDIST.1000.PLINK.IBD = read.table("./SmallPop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.IBD = aggregate(dtaROHSDIST.1000.PLINK.IBD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.IBD$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.IBD$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.IBD$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.IBD$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.IBD$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.IBD = aggregate(dtaROHSDIST.1000.PLINK.IBD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.IBD$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.IBD$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.PLINK.IBD$SimID, REP=dtaROHSDIST.1000.PLINK.IBD$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.IBD = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.IBD$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.IBD$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.IBD$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.IBD = merge(dtaROHSDIST.plot.1000.PLINK.IBD, ROHsDist.SD.1000.PLINK.IBD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.IBD) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.IBD = dtaROHSDIST.plot.1000.PLINK.IBD[order(dtaROHSDIST.plot.1000.PLINK.IBD$CLASS, dtaROHSDIST.plot.1000.PLINK.IBD$PERC),]
dtaROHSDIST.plot.1000.PLINK.IBD$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.IBD$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.IBD);rm(ROHsDist.MeanforCI.1000.PLINK.IBD);rm(ROHsDist.SD.1000.PLINK.IBD)

#Merge ALL and remove duplicated WGS lines
dtaROHSDIST.1000.PLINK = unique(rbind(dtaROHSDIST.plot.1000.PLINK.RAD, dtaROHSDIST.plot.1000.PLINK.ARRAY, dtaROHSDIST.plot.1000.PLINK.IBD))
#Change factor levels
dtaROHSDIST.1000.PLINK$PERC = factor(dtaROHSDIST.1000.PLINK$PERC, levels = c("180000", "240000", "600000","SMALL_ARRAY","LARGE_ARRAY","TRUE_IBD_100"))
#Order the df
dtaROHSDIST.1000.PLINK = dtaROHSDIST.1000.PLINK[order(dtaROHSDIST.1000.PLINK$CLASS, dtaROHSDIST.1000.PLINK$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.1000.PLINK.RAD); rm(dtaROHSDIST.plot.1000.PLINK.ARRAY,dtaROHSDIST.plot.1000.PLINK.IBD)
#Divide TRUE IBD by 1000 (because in basepair while PLINK in kb)
dtaROHSDIST.1000.PLINK$x[dtaROHSDIST.1000.PLINK$PERC %in% c("TRUE_IBD_100")] = dtaROHSDIST.1000.PLINK$x[dtaROHSDIST.1000.PLINK$PERC %in% c("TRUE_IBD_100")]/1000

#### Read the data ROHsDist PLINK100KB Ne 10'000 ####

dtaROHSDIST.10000.PLINK.RAD = read.table("./LargePop/Analyses/RADSEQ_PLINK_defaultPARAM_ROHsDistribution.txt", header = T)
#Take only the perc we are interested in
dtaROHSDIST.10000.PLINK.RAD = dtaROHSDIST.10000.PLINK.RAD[dtaROHSDIST.10000.PLINK.RAD$NB_WIND %in% c("15000", "25000", "60000"),]
colnames(dtaROHSDIST.10000.PLINK.RAD) = c("SimuID","SUB_method","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.RAD = aggregate(dtaROHSDIST.10000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.RAD$SUB_method), FUN = mean)
#Calculate sd
ROHsDist.SD.10000.PLINK.RAD = as.data.frame(as.matrix(aggregate(dtaROHSDIST.10000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.RAD$SUB_method), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.RAD = merge(dtaROHSDIST.plot.10000.PLINK.RAD, ROHsDist.SD.10000.PLINK.RAD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.RAD) = c("CLASS", "PERC", "x", "sd")
#Pass sd to numeric
dtaROHSDIST.plot.10000.PLINK.RAD$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.RAD$sd))
#Order the df
dtaROHSDIST.plot.10000.PLINK.RAD = dtaROHSDIST.plot.10000.PLINK.RAD[order(dtaROHSDIST.plot.10000.PLINK.RAD$CLASS, dtaROHSDIST.plot.10000.PLINK.RAD$PERC),]
dtaROHSDIST.plot.10000.PLINK.RAD$PERC = factor(dtaROHSDIST.plot.10000.PLINK.RAD$PERC)
#rm RAD and SD
rm(dtaROHSDIST.10000.PLINK.RAD); rm(ROHsDist.SD.10000.PLINK.RAD)

#Read the ne 1'000 50k array PLINK
dtaROHSDIST.10000.PLINK.50k = read.table("./LargePop/Analyses/SmallArray_PLINK_defaultPARAM_ROHsDistribution.txt", header = T)
dtaROHSDIST.10000.PLINK.50k$SUB_method = "SMALL_ARRAY"
#Read the ne 1'000 700k array PLINK
dtaROHSDIST.10000.PLINK.700k = read.table("./LargePop/Analyses/LargeArray_PLINK_defaultPARAM_ROHsDistribution.txt", header = T)
dtaROHSDIST.10000.PLINK.700k$SUB_method = "LARGE_ARRAY"
dtaROHSDIST.10000.PLINK.ARRAY = unique(rbind(dtaROHSDIST.10000.PLINK.50k,dtaROHSDIST.10000.PLINK.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.ARRAY = aggregate(dtaROHSDIST.10000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.ARRAY$SUB_method), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.ARRAY = aggregate(dtaROHSDIST.10000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.ARRAY$SUB_method, SimuIDs=dtaROHSDIST.10000.PLINK.ARRAY$SimuID, REP=dtaROHSDIST.10000.PLINK.ARRAY$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.ARRAY = merge(dtaROHSDIST.plot.10000.PLINK.ARRAY, ROHsDist.SD.10000.PLINK.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.ARRAY = dtaROHSDIST.plot.10000.PLINK.ARRAY[order(dtaROHSDIST.plot.10000.PLINK.ARRAY$CLASS, dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC),]
dtaROHSDIST.plot.10000.PLINK.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.50k);rm(dtaROHSDIST.10000.PLINK.700k);rm(dtaROHSDIST.10000.PLINK.ARRAY);rm(ROHsDist.MeanforCI.10000.PLINK.ARRAY);rm(ROHsDist.SD.10000.PLINK.ARRAY)

#Read the ne 1'000 TRUE IBD PLINK
dtaROHSDIST.10000.PLINK.IBD = read.table("./LargePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.IBD = aggregate(dtaROHSDIST.10000.PLINK.IBD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.IBD$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.IBD$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.IBD$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.IBD$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.IBD$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.IBD = aggregate(dtaROHSDIST.10000.PLINK.IBD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.IBD$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.IBD$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.PLINK.IBD$SimID, REP=dtaROHSDIST.10000.PLINK.IBD$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.IBD = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.IBD$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.IBD$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.IBD$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.IBD = merge(dtaROHSDIST.plot.10000.PLINK.IBD, ROHsDist.SD.10000.PLINK.IBD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.IBD) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.IBD = dtaROHSDIST.plot.10000.PLINK.IBD[order(dtaROHSDIST.plot.10000.PLINK.IBD$CLASS, dtaROHSDIST.plot.10000.PLINK.IBD$PERC),]
dtaROHSDIST.plot.10000.PLINK.IBD$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.IBD$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.IBD);rm(ROHsDist.MeanforCI.10000.PLINK.IBD);rm(ROHsDist.SD.10000.PLINK.IBD)

#Merge ALL and remove duplicated WGS lines
dtaROHSDIST.10000.PLINK = unique(rbind(dtaROHSDIST.plot.10000.PLINK.RAD, dtaROHSDIST.plot.10000.PLINK.ARRAY, dtaROHSDIST.plot.10000.PLINK.IBD))
#Change factor levels
dtaROHSDIST.10000.PLINK$PERC = factor(dtaROHSDIST.10000.PLINK$PERC, levels = c("15000", "25000", "60000","SMALL_ARRAY","LARGE_ARRAY","TRUE_IBD_100"))
#Order the df
dtaROHSDIST.10000.PLINK = dtaROHSDIST.10000.PLINK[order(dtaROHSDIST.10000.PLINK$CLASS, dtaROHSDIST.10000.PLINK$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.10000.PLINK.RAD); rm(dtaROHSDIST.plot.10000.PLINK.ARRAY,dtaROHSDIST.plot.10000.PLINK.IBD)
#Divide TRUE IBD by 10000 (because in basepair while PLINK in kb)
dtaROHSDIST.10000.PLINK$x[dtaROHSDIST.10000.PLINK$PERC %in% c("TRUE_IBD_100")] = dtaROHSDIST.10000.PLINK$x[dtaROHSDIST.10000.PLINK$PERC %in% c("TRUE_IBD_100")]/1000

#Define ROHs classes
classesROH = c("< 2", "2 - 4","4 - 6","6 - 10","10 - 16", "> 16")

pdf("~/PhD/ROH/Which_data/Manuscript/SUBMITTED_VERSION/Review_round_I/FIGURES/SUPPLEMENTARY_FIGURES/FIGURES6.pdf", width = 55, height = 23)

#Layout for having all plots
layout(matrix(c(rep(0,5),0,1,0,2,0,rep(0,5)), 3, 5, byrow = T), heights = c(.5,1,.8), widths = c(.5,2,1,2,1))

#Par for margin and outer margings
par(mar = c(3.5,3.5,3.5,3.5), oma = c(.2,10,.2,.2))

#### PANNEL A: ROHsDist Ne 1'000 PLINK 100KB RAD + ARRAYS ####

#Define colors for subsampled dataset
subCOLORS = c("#82EEFD", "#0492C2", "#016064", "#2832C2", "#0A1172")

#Create empty plot
plot(dtaROHSDIST.1000.PLINK$CLASS, dtaROHSDIST.1000.PLINK$x/1000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0,600))

#Loop through classes to draw boxplots
for(class in sort(unique(dtaROHSDIST.1000.PLINK$CLASS))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.35, class-0.21,class-0.07, class+0.07, class+0.21,class+0.35)
  
  #Loop through reduced representation
  for(reduced in 1:5){
    
    #subset the data
    dta_sub = dtaROHSDIST.1000.PLINK[dtaROHSDIST.1000.PLINK$CLASS == class & dtaROHSDIST.1000.PLINK$PERC == levels(dtaROHSDIST.1000.PLINK$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$CLASS == class & dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_100"],
         ytop = dta_sub$x/1000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.07), x1 = (polyXcoord[reduced] + 0.07), y0 = ((dta_sub$x - dta_sub$sd)/1000), y1 = ((dta_sub$x + dta_sub$sd)/1000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.1000.PLINK$CLASS == class],
             y1 = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.1000.PLINK$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add the y axis
axis(2, cex.axis = 4, hadj = 1.5, las = 1, cex.axis = 5, lwd.ticks = 5, tck = -0.03)
#Add y axis title
mtext(text = "Mean Mb ± SD (among individuals)", side = 2, cex = 3.5, line = 20)
#Add x axis title
mtext(text = "ROHs Length Classes [Mb]", side = 1, cex = 4, line = 15)
#Add panel lab
mtext(text = "A", side = 3, cex = 7, line = 15, at = -2)
mtext(text = "PLINK", side = 3, cex = 5, line = 15, at = 3)

#### PANNEL B: ROHsDist Ne 10'000 PLINK 100KB RAD + ARRAYS ####

subCOLORS = c("#FFBF00", "#FD6A02", "#FDA172","#CC7722", "#793802", "grey18")

#Create empty plot
plot(dtaROHSDIST.10000.PLINK$CLASS, dtaROHSDIST.10000.PLINK$x/1000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0,600))

#Loop through classes to draw boxplots
for(class in sort(unique(dtaROHSDIST.10000.PLINK$CLASS))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.35, class-0.21,class-0.07, class+0.07, class+0.21,class+0.35)
  
  #Loop through reduced representation
  for(reduced in 1:5){
    
    #subset the data
    dta_sub = dtaROHSDIST.10000.PLINK[dtaROHSDIST.10000.PLINK$CLASS == class & dtaROHSDIST.10000.PLINK$PERC == levels(dtaROHSDIST.10000.PLINK$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$CLASS == class & dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_100"],
         ytop = dta_sub$x/1000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.07), x1 = (polyXcoord[reduced] + 0.07), y0 = ((dta_sub$x - dta_sub$sd)/1000), y1 = ((dta_sub$x + dta_sub$sd)/1000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.10000.PLINK$CLASS == class],
             y1 = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.10000.PLINK$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add the y axis
axis(2, cex.axis = 4, hadj = 1.5, las = 1, cex.axis = 5, lwd.ticks = 5, tck = -0.03)
#Add y axis title
mtext(text = "Mean Mb ± SD (among individuals)", side = 2, cex = 3.5, line = 20)
#Add x axis title
mtext(text = "ROHs Length Classes [Mb]", side = 1, cex = 4, line = 15)
#Add panel lab
mtext(text = "B", side = 3, cex = 7, line = 15, at = -2)
mtext(text = "PLINK", side = 3, cex = 5, line = 15, at = 3)

#Close the plot
dev.off()

#######################################################
############# FIGURE S7 TFPN Rates 1000 GEN ###########
#######################################################

## 1K PLINK

#Read the data TFPN Rates Ne 1'000 PLINK 100KB RAD
dtaTFPN.1000.100KB.RAD = read.table("./SmallPop/Analyses/RADSEQ_PLINK_TFPN_Rates_GEN1000.txt", header = T)
dtaTFPN.1000.100KB.RAD = dtaTFPN.1000.100KB.RAD[dtaTFPN.1000.100KB.RAD$NB_RAD_FRAG %in% c(180000, 240000, 600000),]
#Define Percentage genome sequenced Ne 1'000
PercSeq.1000.100KB.RAD = sort(round((as.numeric(as.character(unique(dtaTFPN.1000.100KB.RAD$NB_RAD_FRAG)))*500)/30000000, digits = 0))
#Get a vector with mean Rates per RAD_Frag
TP.1000.100KB.RAD = aggregate(dtaTFPN.1000.100KB.RAD$TP, by = list(dtaTFPN.1000.100KB.RAD$NB_RAD_FRAG), FUN = mean)
TN.1000.100KB.RAD = aggregate(dtaTFPN.1000.100KB.RAD$TN, by = list(dtaTFPN.1000.100KB.RAD$NB_RAD_FRAG), FUN = mean)
FP.1000.100KB.RAD  = aggregate(dtaTFPN.1000.100KB.RAD$FP, by = list(dtaTFPN.1000.100KB.RAD$NB_RAD_FRAG), FUN = mean)
FN.1000.100KB.RAD  = aggregate(dtaTFPN.1000.100KB.RAD$FN, by = list(dtaTFPN.1000.100KB.RAD$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one df
Rates.1000.100KB.RAD = rbind(cbind(TP.1000.100KB.RAD, Rate=rep("TP", nrow(TP.1000.100KB.RAD))), cbind(TN.1000.100KB.RAD, Rate=rep("TN", nrow(TN.1000.100KB.RAD))), cbind(FP.1000.100KB.RAD, Rate=rep("FP", nrow(FP.1000.100KB.RAD))),
                             cbind(FN.1000.100KB.RAD, Rate=rep("FN", nrow(FN.1000.100KB.RAD))))
#Change names
colnames(Rates.1000.100KB.RAD) = c("PercSeq", "Value", "Rate")
#Order df per Perc of genome sequenced
Rates.1000.100KB.RAD = Rates.1000.100KB.RAD[order(Rates.1000.100KB.RAD$PercSeq),]
#Create the datafarme
stackedRates.1000.PLINK100KB.RAD = t(data.frame(Rates.1000.100KB.RAD$Value[Rates.1000.100KB.RAD$Rate == "TN"], Rates.1000.100KB.RAD$Value[Rates.1000.100KB.RAD$Rate == "TP"],
                                                Rates.1000.100KB.RAD$Value[Rates.1000.100KB.RAD$Rate == "FN"], Rates.1000.100KB.RAD$Value[Rates.1000.100KB.RAD$Rate == "FP"]))
row.names(stackedRates.1000.PLINK100KB.RAD) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.1000.PLINK100KB.RAD) = PercSeq.1000.100KB.RAD

## Read the data TFPN Rates Ne 1'000 PLINK 100KB ARRAY ##

PercSeq.ARRAY = c("SMALL\nARRAY", "LARGE\nARRAY")

#Read the data
dtaTFPN.1000.100KB.50k = read.table("./SmallPop/Analyses/SmallArray_PLINK_TFPN_Rates_GEN1000.txt", header = T)
dtaTFPN.1000.100KB.700k = read.table("./SmallPop/Analyses/LargeArray_PLINK_TFPN_Rates_GEN1000.txt", header = T)
dtaTFPN.1000.100KB.ARRAY = rbind(dtaTFPN.1000.100KB.50k, dtaTFPN.1000.100KB.700k)
#Get a vector with mean Rates per RAD_Frag
TP.1000.100KB.ARRAY = aggregate(dtaTFPN.1000.100KB.ARRAY$TP, by = list(dtaTFPN.1000.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
TN.1000.100KB.ARRAY = aggregate(dtaTFPN.1000.100KB.ARRAY$TN, by = list(dtaTFPN.1000.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
FP.1000.100KB.ARRAY = aggregate(dtaTFPN.1000.100KB.ARRAY$FP, by = list(dtaTFPN.1000.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
FN.1000.100KB.ARRAY = aggregate(dtaTFPN.1000.100KB.ARRAY$FN, by = list(dtaTFPN.1000.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one df
Rates.1000.100KB.ARRAY = rbind(cbind(TP.1000.100KB.ARRAY, Rate=rep("TP", nrow(TP.1000.100KB.ARRAY))), cbind(TN.1000.100KB.ARRAY, Rate=rep("TN", nrow(TN.1000.100KB.ARRAY))),
                               cbind(FP.1000.100KB.ARRAY, Rate=rep("FP", nrow(FP.1000.100KB.ARRAY))), cbind(FN.1000.100KB.ARRAY, Rate=rep("FN", nrow(FN.1000.100KB.ARRAY))))
#Change names
colnames(Rates.1000.100KB.ARRAY) = c("PercSeq", "Value", "Rate")
Rates.1000.100KB.ARRAY$PercSeq = factor(Rates.1000.100KB.ARRAY$PercSeq, levels = c("SmallArray", "LargeArray"))
#Order data
Rates.1000.100KB.ARRAY = Rates.1000.100KB.ARRAY[order(Rates.1000.100KB.ARRAY$PercSeq),]
#Create the datafarme
stackedRates.1000.100KB.ARRAY = t(data.frame(Rates.1000.100KB.ARRAY$Value[Rates.1000.100KB.ARRAY$Rate == "TN"], Rates.1000.100KB.ARRAY$Value[Rates.1000.100KB.ARRAY$Rate == "TP"],
                                             Rates.1000.100KB.ARRAY$Value[Rates.1000.100KB.ARRAY$Rate == "FN"], Rates.1000.100KB.ARRAY$Value[Rates.1000.100KB.ARRAY$Rate == "FP"]))
row.names(stackedRates.1000.100KB.ARRAY) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.1000.100KB.ARRAY) = PercSeq.ARRAY

## Read the data TFPN Rates Ne 1'000 PLINK 100KB WGS ##

PercSeq.WGS = "WGS"

#Read the data
dtaTFPN.1000.100KB.WGS = read.table("./SmallPop/Analyses/WGS_PLINK_TFPN_Rates_GEN1000.txt", header = F)
colnames(dtaTFPN.1000.100KB.WGS) = c("SimID","NB_RAD_FRAG","INDIV","TP","TN","FP","FN" )
#Get a vector with mean Rates per RAD_Frag
TP.1000.100KB.WGS = aggregate(dtaTFPN.1000.100KB.WGS$TP, by = list(dtaTFPN.1000.100KB.WGS$NB_RAD_FRAG), FUN = mean)
TN.1000.100KB.WGS = aggregate(dtaTFPN.1000.100KB.WGS$TN, by = list(dtaTFPN.1000.100KB.WGS$NB_RAD_FRAG), FUN = mean)
FP.1000.100KB.WGS = aggregate(dtaTFPN.1000.100KB.WGS$FP, by = list(dtaTFPN.1000.100KB.WGS$NB_RAD_FRAG), FUN = mean)
FN.1000.100KB.WGS = aggregate(dtaTFPN.1000.100KB.WGS$FN, by = list(dtaTFPN.1000.100KB.WGS$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one dfWGS
Rates.1000.100KB.WGS = rbind(cbind(TP.1000.100KB.WGS, Rate=rep("TP", nrow(TP.1000.100KB.WGS))), cbind(TN.1000.100KB.WGS, Rate=rep("TN", nrow(TN.1000.100KB.WGS))),
                             cbind(FP.1000.100KB.WGS, Rate=rep("FP", nrow(FP.1000.100KB.WGS))), cbind(FN.1000.100KB.WGS, Rate=rep("FN", nrow(FN.1000.100KB.WGS))))
#Change names
colnames(Rates.1000.100KB.WGS) = c("PercSeq", "Value", "Rate")
#Create the datafarme
stackedRates.1000.100KB.WGS = t(data.frame(Rates.1000.100KB.WGS$Value[Rates.1000.100KB.WGS$Rate == "TN"], Rates.1000.100KB.WGS$Value[Rates.1000.100KB.WGS$Rate == "TP"],
                                           Rates.1000.100KB.WGS$Value[Rates.1000.100KB.WGS$Rate == "FN"], Rates.1000.100KB.WGS$Value[Rates.1000.100KB.WGS$Rate == "FP"]))
row.names(stackedRates.1000.100KB.WGS) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.1000.100KB.WGS) = PercSeq.WGS

#Merge WGS, arrays and RAD
stackedRates.1000.100KB = cbind(stackedRates.1000.PLINK100KB.RAD, stackedRates.1000.100KB.ARRAY, stackedRates.1000.100KB.WGS)

#rm the rest
rm(dtaTFPN.1000.100KB.RAD); rm(dtaTFPN.1000.100KB.ARRAY); rm(dtaTFPN.1000.100KB.50k); rm(dtaTFPN.1000.100KB.700k); rm(dtaTFPN.1000.100KB.WGS)
rm(FN.1000.100KB.RAD); rm(FP.1000.100KB.RAD); rm(TN.1000.100KB.RAD); rm(TP.1000.100KB.RAD)
rm(FN.1000.100KB.ARRAY); rm(FP.1000.100KB.ARRAY); rm(TN.1000.100KB.ARRAY); rm(TP.1000.100KB.ARRAY)
rm(FN.1000.100KB.WGS); rm(FP.1000.100KB.WGS); rm(TN.1000.100KB.WGS); rm(TP.1000.100KB.WGS)
rm(Rates.1000.100KB.RAD);rm(Rates.1000.100KB.ARRAY);rm(Rates.1000.100KB.WGS)
rm(stackedRates.1000.100KB.ARRAY);rm(stackedRates.1000.PLINK100KB.RAD);rm(stackedRates.1000.100KB.WGS)

## 10K PLINK

#Read the data TFPN Rates Ne 10'000 PLINK 100KB RAD
dtaTFPN.10000.100KB.RAD = read.table("./LargePop/Analyses/RADSEQ_PLINK_TFPN_Rates_GEN1000.txt", header = T)
dtaTFPN.10000.100KB.RAD = dtaTFPN.10000.100KB.RAD[dtaTFPN.10000.100KB.RAD$NB_RAD_FRAG %in% c(15000, 25000, 60000),]
#Define Percentage genome sequenced Ne 1'000
PercSeq.10000.100KB.RAD = sort(round((as.numeric(as.character(unique(dtaTFPN.10000.100KB.RAD$NB_RAD_FRAG)))*500)/30000000, digits = 3))
PercSeq.10000.100KB.RAD[2] = 0.4
#Get a vector with mean Rates per RAD_Frag
TP.10000.100KB.RAD = aggregate(dtaTFPN.10000.100KB.RAD$TP, by = list(dtaTFPN.10000.100KB.RAD$NB_RAD_FRAG), FUN = mean)
TN.10000.100KB.RAD = aggregate(dtaTFPN.10000.100KB.RAD$TN, by = list(dtaTFPN.10000.100KB.RAD$NB_RAD_FRAG), FUN = mean)
FP.10000.100KB.RAD  = aggregate(dtaTFPN.10000.100KB.RAD$FP, by = list(dtaTFPN.10000.100KB.RAD$NB_RAD_FRAG), FUN = mean)
FN.10000.100KB.RAD  = aggregate(dtaTFPN.10000.100KB.RAD$FN, by = list(dtaTFPN.10000.100KB.RAD$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one df
Rates.10000.100KB.RAD = rbind(cbind(TP.10000.100KB.RAD, Rate=rep("TP", nrow(TP.10000.100KB.RAD))), cbind(TN.10000.100KB.RAD, Rate=rep("TN", nrow(TN.10000.100KB.RAD))), cbind(FP.10000.100KB.RAD, Rate=rep("FP", nrow(FP.10000.100KB.RAD))),
                              cbind(FN.10000.100KB.RAD, Rate=rep("FN", nrow(FN.10000.100KB.RAD))))
#Change names
colnames(Rates.10000.100KB.RAD) = c("PercSeq", "Value", "Rate")
#Order df per Perc of genome sequenced
Rates.10000.100KB.RAD = Rates.10000.100KB.RAD[order(Rates.10000.100KB.RAD$PercSeq),]
#Create the datafarme
stackedRates.10000.PLINK100KB.RAD = t(data.frame(Rates.10000.100KB.RAD$Value[Rates.10000.100KB.RAD$Rate == "TN"], Rates.10000.100KB.RAD$Value[Rates.10000.100KB.RAD$Rate == "TP"],
                                                 Rates.10000.100KB.RAD$Value[Rates.10000.100KB.RAD$Rate == "FN"], Rates.10000.100KB.RAD$Value[Rates.10000.100KB.RAD$Rate == "FP"]))
row.names(stackedRates.10000.PLINK100KB.RAD) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.10000.PLINK100KB.RAD) = PercSeq.10000.100KB.RAD

## Read the data TFPN Rates Ne 1'000 PLINK 100KB ARRAY ##

PercSeq.ARRAY = c("SMALL\nARRAY", "LARGE\nARRAY")

#Read the data
dtaTFPN.10000.100KB.50k = read.table("./LargePop/Analyses/SmallArray_PLINK_TFPN_Rates_GEN1000.txt", header = T)
dtaTFPN.10000.100KB.700k = read.table("./LargePop/Analyses/LargeArray_PLINK_TFPN_Rates_GEN1000.txt", header = T)
dtaTFPN.10000.100KB.ARRAY = rbind(dtaTFPN.10000.100KB.50k, dtaTFPN.10000.100KB.700k)
#Get a vector with mean Rates per RAD_Frag
TP.10000.100KB.ARRAY = aggregate(dtaTFPN.10000.100KB.ARRAY$TP, by = list(dtaTFPN.10000.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
TN.10000.100KB.ARRAY = aggregate(dtaTFPN.10000.100KB.ARRAY$TN, by = list(dtaTFPN.10000.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
FP.10000.100KB.ARRAY = aggregate(dtaTFPN.10000.100KB.ARRAY$FP, by = list(dtaTFPN.10000.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
FN.10000.100KB.ARRAY = aggregate(dtaTFPN.10000.100KB.ARRAY$FN, by = list(dtaTFPN.10000.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one df
Rates.10000.100KB.ARRAY = rbind(cbind(TP.10000.100KB.ARRAY, Rate=rep("TP", nrow(TP.10000.100KB.ARRAY))), cbind(TN.10000.100KB.ARRAY, Rate=rep("TN", nrow(TN.10000.100KB.ARRAY))),
                                cbind(FP.10000.100KB.ARRAY, Rate=rep("FP", nrow(FP.10000.100KB.ARRAY))), cbind(FN.10000.100KB.ARRAY, Rate=rep("FN", nrow(FN.10000.100KB.ARRAY))))
#Change names
colnames(Rates.10000.100KB.ARRAY) = c("PercSeq", "Value", "Rate")
Rates.10000.100KB.ARRAY$PercSeq = factor(Rates.10000.100KB.ARRAY$PercSeq, levels = c("SmallArray", "LargeArray"))
#Order data
Rates.10000.100KB.ARRAY = Rates.10000.100KB.ARRAY[order(Rates.10000.100KB.ARRAY$PercSeq),]
#Create the datafarme
stackedRates.10000.100KB.ARRAY = t(data.frame(Rates.10000.100KB.ARRAY$Value[Rates.10000.100KB.ARRAY$Rate == "TN"], Rates.10000.100KB.ARRAY$Value[Rates.10000.100KB.ARRAY$Rate == "TP"],
                                              Rates.10000.100KB.ARRAY$Value[Rates.10000.100KB.ARRAY$Rate == "FN"], Rates.10000.100KB.ARRAY$Value[Rates.10000.100KB.ARRAY$Rate == "FP"]))
row.names(stackedRates.10000.100KB.ARRAY) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.10000.100KB.ARRAY) = PercSeq.ARRAY

## Read the data TFPN Rates Ne 1'000 PLINK 100KB WGS ##

PercSeq.WGS = "WGS"

#Read the data
dtaTFPN.10000.100KB.WGS = read.table("./LargePop/Analyses/WGS_PLINK_TFPN_Rates_GEN1000.txt", header = T)
#Get a vector with mean Rates per RAD_Frag
TP.10000.100KB.WGS = aggregate(dtaTFPN.10000.100KB.WGS$TP, by = list(dtaTFPN.10000.100KB.WGS$NB_RAD_FRAG), FUN = mean)
TN.10000.100KB.WGS = aggregate(dtaTFPN.10000.100KB.WGS$TN, by = list(dtaTFPN.10000.100KB.WGS$NB_RAD_FRAG), FUN = mean)
FP.10000.100KB.WGS = aggregate(dtaTFPN.10000.100KB.WGS$FP, by = list(dtaTFPN.10000.100KB.WGS$NB_RAD_FRAG), FUN = mean)
FN.10000.100KB.WGS = aggregate(dtaTFPN.10000.100KB.WGS$FN, by = list(dtaTFPN.10000.100KB.WGS$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one dfWGS
Rates.10000.100KB.WGS = rbind(cbind(TP.10000.100KB.WGS, Rate=rep("TP", nrow(TP.10000.100KB.WGS))), cbind(TN.10000.100KB.WGS, Rate=rep("TN", nrow(TN.10000.100KB.WGS))),
                              cbind(FP.10000.100KB.WGS, Rate=rep("FP", nrow(FP.10000.100KB.WGS))), cbind(FN.10000.100KB.WGS, Rate=rep("FN", nrow(FN.10000.100KB.WGS))))
#Change names
colnames(Rates.10000.100KB.WGS) = c("PercSeq", "Value", "Rate")
#Create the datafarme
stackedRates.10000.100KB.WGS = t(data.frame(Rates.10000.100KB.WGS$Value[Rates.10000.100KB.WGS$Rate == "TN"], Rates.10000.100KB.WGS$Value[Rates.10000.100KB.WGS$Rate == "TP"],
                                            Rates.10000.100KB.WGS$Value[Rates.10000.100KB.WGS$Rate == "FN"], Rates.10000.100KB.WGS$Value[Rates.10000.100KB.WGS$Rate == "FP"]))
row.names(stackedRates.10000.100KB.WGS) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.10000.100KB.WGS) = PercSeq.WGS

#Merge WGS, arrays and RAD
stackedRates.10000.100KB = cbind(stackedRates.10000.PLINK100KB.RAD, stackedRates.10000.100KB.ARRAY, stackedRates.10000.100KB.WGS)

#rm the rest
rm(dtaTFPN.10000.100KB.RAD); rm(dtaTFPN.10000.100KB.ARRAY); rm(dtaTFPN.10000.100KB.50k); rm(dtaTFPN.10000.100KB.700k); rm(dtaTFPN.10000.100KB.WGS)
rm(FN.10000.100KB.RAD); rm(FP.10000.100KB.RAD); rm(TN.10000.100KB.RAD); rm(TP.10000.100KB.RAD)
rm(FN.10000.100KB.ARRAY); rm(FP.10000.100KB.ARRAY); rm(TN.10000.100KB.ARRAY); rm(TP.10000.100KB.ARRAY)
rm(FN.10000.100KB.WGS); rm(FP.10000.100KB.WGS); rm(TN.10000.100KB.WGS); rm(TP.10000.100KB.WGS)
rm(Rates.10000.100KB.RAD);rm(Rates.10000.100KB.ARRAY);rm(Rates.10000.100KB.WGS)
rm(stackedRates.10000.100KB.ARRAY);rm(stackedRates.10000.PLINK100KB.RAD);rm(stackedRates.10000.100KB.WGS)

## 1K RZooRoH

#Read the data TFPN Rates Ne 1'000 RzooRoH 3R RAD
dtaTFPN.1000.Rzoo.RAD = read.table("./SmallPop/Analyses/RADSEQ_RZooRoH_TFPN_Rates_GEN1000.txt", header = F)
colnames(dtaTFPN.1000.Rzoo.RAD) = c("SimID","NB_RAD_FRAG", "REP","INDIV","TP","TN","FP","FN")
dtaTFPN.1000.Rzoo.RAD = dtaTFPN.1000.Rzoo.RAD[dtaTFPN.1000.Rzoo.RAD$NB_RAD_FRAG %in% c(3000, 9000, 120000),]
dtaTFPN.1000.Rzoo.RAD[,2] = as.numeric(dtaTFPN.1000.Rzoo.RAD[,2])
dtaTFPN.1000.Rzoo.RAD[,3] = as.numeric(dtaTFPN.1000.Rzoo.RAD[,3])
dtaTFPN.1000.Rzoo.RAD[,4] = as.numeric(dtaTFPN.1000.Rzoo.RAD[,4])
dtaTFPN.1000.Rzoo.RAD[,5] = as.numeric(dtaTFPN.1000.Rzoo.RAD[,5])
dtaTFPN.1000.Rzoo.RAD[,6] = as.numeric(dtaTFPN.1000.Rzoo.RAD[,6])
dtaTFPN.1000.Rzoo.RAD[,7] = as.numeric(dtaTFPN.1000.Rzoo.RAD[,7])
dtaTFPN.1000.Rzoo.RAD[,8] = as.numeric(dtaTFPN.1000.Rzoo.RAD[,8])
#Define Percentage genome sequenced Ne 1'000
PercSeq.1000.Rzoo.RAD = sort(round((as.numeric(as.character(unique(dtaTFPN.1000.Rzoo.RAD$NB_RAD_FRAG)))*500)/30000000, digits = 3))
#Get a vector with mean Rates per RAD_Frag
TP.1000.Rzoo.RAD = aggregate(dtaTFPN.1000.Rzoo.RAD$TP, by = list(dtaTFPN.1000.Rzoo.RAD$NB_RAD_FRAG), FUN = mean)
TN.1000.Rzoo.RAD = aggregate(dtaTFPN.1000.Rzoo.RAD$TN, by = list(dtaTFPN.1000.Rzoo.RAD$NB_RAD_FRAG), FUN = mean)
FP.1000.Rzoo.RAD  = aggregate(dtaTFPN.1000.Rzoo.RAD$FP, by = list(dtaTFPN.1000.Rzoo.RAD$NB_RAD_FRAG), FUN = mean)
FN.1000.Rzoo.RAD  = aggregate(dtaTFPN.1000.Rzoo.RAD$FN, by = list(dtaTFPN.1000.Rzoo.RAD$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one df
Rates.1000.Rzoo.RAD = rbind(cbind(TP.1000.Rzoo.RAD, Rate=rep("TP", nrow(TP.1000.Rzoo.RAD))), cbind(TN.1000.Rzoo.RAD, Rate=rep("TN", nrow(TN.1000.Rzoo.RAD))), cbind(FP.1000.Rzoo.RAD, Rate=rep("FP", nrow(FP.1000.Rzoo.RAD))),
                            cbind(FN.1000.Rzoo.RAD, Rate=rep("FN", nrow(FN.1000.Rzoo.RAD))))
#Change names
colnames(Rates.1000.Rzoo.RAD) = c("PercSeq", "Value", "Rate")
#Order df per Perc of genome sequenced
Rates.1000.Rzoo.RAD = Rates.1000.Rzoo.RAD[order(Rates.1000.Rzoo.RAD$PercSeq),]
#Create the datafarme
stackedRates.1000.Rzoo.RAD = t(data.frame(Rates.1000.Rzoo.RAD$Value[Rates.1000.Rzoo.RAD$Rate == "TN"], Rates.1000.Rzoo.RAD$Value[Rates.1000.Rzoo.RAD$Rate == "TP"],
                                          Rates.1000.Rzoo.RAD$Value[Rates.1000.Rzoo.RAD$Rate == "FN"], Rates.1000.Rzoo.RAD$Value[Rates.1000.Rzoo.RAD$Rate == "FP"]))
row.names(stackedRates.1000.Rzoo.RAD) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.1000.Rzoo.RAD) = PercSeq.1000.Rzoo.RAD

## Read the data TFPN Rates Ne 1'000 RzooROH 3R ARRAY ##

#Read the data
dtaTFPN.1000.Rzoo.50k = read.table("./SmallPop/Analyses/SmallArray_RZooRoH_TFPN_Rates_GEN1000.txt", header = T)
dtaTFPN.1000.Rzoo.700k = read.table("./SmallPop/Analyses/LargeArray_RZooRoH_TFPN_Rates_GEN1000.txt", header = T)
dtaTFPN.1000.Rzoo.ARRAY = rbind(dtaTFPN.1000.Rzoo.50k, dtaTFPN.1000.Rzoo.700k)
#Get a vector with mean Rates per RAD_Frag
TP.1000.Rzoo.ARRAY = aggregate(dtaTFPN.1000.Rzoo.ARRAY$TP, by = list(dtaTFPN.1000.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
TN.1000.Rzoo.ARRAY = aggregate(dtaTFPN.1000.Rzoo.ARRAY$TN, by = list(dtaTFPN.1000.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
FP.1000.Rzoo.ARRAY = aggregate(dtaTFPN.1000.Rzoo.ARRAY$FP, by = list(dtaTFPN.1000.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
FN.1000.Rzoo.ARRAY = aggregate(dtaTFPN.1000.Rzoo.ARRAY$FN, by = list(dtaTFPN.1000.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one df
Rates.1000.Rzoo.ARRAY = rbind(cbind(TP.1000.Rzoo.ARRAY, Rate=rep("TP", nrow(TP.1000.Rzoo.ARRAY))), cbind(TN.1000.Rzoo.ARRAY, Rate=rep("TN", nrow(TN.1000.Rzoo.ARRAY))),
                              cbind(FP.1000.Rzoo.ARRAY, Rate=rep("FP", nrow(FP.1000.Rzoo.ARRAY))), cbind(FN.1000.Rzoo.ARRAY, Rate=rep("FN", nrow(FN.1000.Rzoo.ARRAY))))
#Change names
colnames(Rates.1000.Rzoo.ARRAY) = c("PercSeq", "Value", "Rate")
Rates.1000.Rzoo.ARRAY$PercSeq = factor(Rates.1000.Rzoo.ARRAY$PercSeq, levels = c("SmallArray", "LargeArray"))
#Order data
Rates.1000.Rzoo.ARRAY = Rates.1000.Rzoo.ARRAY[order(Rates.1000.Rzoo.ARRAY$PercSeq),]
#Create the datafarme
stackedRates.1000.Rzoo.ARRAY = t(data.frame(Rates.1000.Rzoo.ARRAY$Value[Rates.1000.Rzoo.ARRAY$Rate == "TN"], Rates.1000.Rzoo.ARRAY$Value[Rates.1000.Rzoo.ARRAY$Rate == "TP"],
                                            Rates.1000.Rzoo.ARRAY$Value[Rates.1000.Rzoo.ARRAY$Rate == "FN"], Rates.1000.Rzoo.ARRAY$Value[Rates.1000.Rzoo.ARRAY$Rate == "FP"]))
row.names(stackedRates.1000.Rzoo.ARRAY) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.1000.Rzoo.ARRAY) = PercSeq.ARRAY

## Read the data TFPN Rates Ne 1'000 RZooROH WGS ##

PercSeq.WGS = "WGS"

#Read the data
dtaTFPN.1000.Rzoo.WGS = read.table("./SmallPop/Analyses/WGS_RZooRoH_TFPN_Rates_GEN1000.txt", header = T)
#Get a vector with mean Rates per RAD_Frag
TP.1000.Rzoo.WGS = aggregate(dtaTFPN.1000.Rzoo.WGS$TP, by = list(dtaTFPN.1000.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
TN.1000.Rzoo.WGS = aggregate(dtaTFPN.1000.Rzoo.WGS$TN, by = list(dtaTFPN.1000.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
FP.1000.Rzoo.WGS = aggregate(dtaTFPN.1000.Rzoo.WGS$FP, by = list(dtaTFPN.1000.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
FN.1000.Rzoo.WGS = aggregate(dtaTFPN.1000.Rzoo.WGS$FN, by = list(dtaTFPN.1000.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one dfWGS
Rates.1000.Rzoo.WGS = rbind(cbind(TP.1000.Rzoo.WGS, Rate=rep("TP", nrow(TP.1000.Rzoo.WGS))), cbind(TN.1000.Rzoo.WGS, Rate=rep("TN", nrow(TN.1000.Rzoo.WGS))),
                            cbind(FP.1000.Rzoo.WGS, Rate=rep("FP", nrow(FP.1000.Rzoo.WGS))), cbind(FN.1000.Rzoo.WGS, Rate=rep("FN", nrow(FN.1000.Rzoo.WGS))))
#Change names
colnames(Rates.1000.Rzoo.WGS) = c("PercSeq", "Value", "Rate")
#Create the datafarme
stackedRates.1000.Rzoo.WGS = t(data.frame(Rates.1000.Rzoo.WGS$Value[Rates.1000.Rzoo.WGS$Rate == "TN"], Rates.1000.Rzoo.WGS$Value[Rates.1000.Rzoo.WGS$Rate == "TP"],
                                          Rates.1000.Rzoo.WGS$Value[Rates.1000.Rzoo.WGS$Rate == "FN"], Rates.1000.Rzoo.WGS$Value[Rates.1000.Rzoo.WGS$Rate == "FP"]))
row.names(stackedRates.1000.Rzoo.WGS) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.1000.Rzoo.WGS) = PercSeq.WGS

#Merge WGS, arrays and RAD
stackedRates.1000.Rzoo = cbind(stackedRates.1000.Rzoo.RAD, stackedRates.1000.Rzoo.ARRAY, stackedRates.1000.Rzoo.WGS)

#rm the rest
rm(dtaTFPN.1000.Rzoo.RAD); rm(dtaTFPN.1000.Rzoo.ARRAY); rm(dtaTFPN.1000.Rzoo.50k); rm(dtaTFPN.1000.Rzoo.700k); rm(dtaTFPN.1000.Rzoo.WGS)
rm(FN.1000.Rzoo.RAD); rm(FP.1000.Rzoo.RAD); rm(TN.1000.Rzoo.RAD); rm(TP.1000.Rzoo.RAD)
rm(FN.1000.Rzoo.ARRAY); rm(FP.1000.Rzoo.ARRAY); rm(TN.1000.Rzoo.ARRAY); rm(TP.1000.Rzoo.ARRAY)
rm(FN.1000.Rzoo.WGS); rm(FP.1000.Rzoo.WGS); rm(TN.1000.Rzoo.WGS); rm(TP.1000.Rzoo.WGS)
rm(Rates.1000.Rzoo.RAD);rm(Rates.1000.Rzoo.ARRAY);rm(Rates.1000.Rzoo.WGS)
rm(stackedRates.1000.Rzoo.ARRAY);rm(stackedRates.1000.Rzoo.RAD);rm(stackedRates.1000.Rzoo.WGS)


## 10K RZooRoH

#Read the data TFPN Rates Ne 1'000 RzooRoH 3R RAD
dtaTFPN.10000.Rzoo.RAD = read.table("./LargePop/Analyses/RADSEQ_RZooRoH_TFPN_Rates_GEN1000.txt", header = T)
dtaTFPN.10000.Rzoo.RAD = dtaTFPN.10000.Rzoo.RAD[dtaTFPN.10000.Rzoo.RAD$NB_RAD_FRAG %in% c(250, 500, 7500),]
#Define Percentage genome sequenced Ne 1'000
PercSeq.10000.Rzoo.RAD = sort(round((as.numeric(as.character(unique(dtaTFPN.10000.Rzoo.RAD$NB_RAD_FRAG)))*500)/30000000, digits = 3))
#Get a vector with mean Rates per RAD_Frag
TP.10000.Rzoo.RAD = aggregate(dtaTFPN.10000.Rzoo.RAD$TP, by = list(dtaTFPN.10000.Rzoo.RAD$NB_RAD_FRAG), FUN = mean)
TN.10000.Rzoo.RAD = aggregate(dtaTFPN.10000.Rzoo.RAD$TN, by = list(dtaTFPN.10000.Rzoo.RAD$NB_RAD_FRAG), FUN = mean)
FP.10000.Rzoo.RAD  = aggregate(dtaTFPN.10000.Rzoo.RAD$FP, by = list(dtaTFPN.10000.Rzoo.RAD$NB_RAD_FRAG), FUN = mean)
FN.10000.Rzoo.RAD  = aggregate(dtaTFPN.10000.Rzoo.RAD$FN, by = list(dtaTFPN.10000.Rzoo.RAD$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one df
Rates.10000.Rzoo.RAD = rbind(cbind(TP.10000.Rzoo.RAD, Rate=rep("TP", nrow(TP.10000.Rzoo.RAD))), cbind(TN.10000.Rzoo.RAD, Rate=rep("TN", nrow(TN.10000.Rzoo.RAD))), cbind(FP.10000.Rzoo.RAD, Rate=rep("FP", nrow(FP.10000.Rzoo.RAD))),
                             cbind(FN.10000.Rzoo.RAD, Rate=rep("FN", nrow(FN.10000.Rzoo.RAD))))
#Change names
colnames(Rates.10000.Rzoo.RAD) = c("PercSeq", "Value", "Rate")
#Order df per Perc of genome sequenced
Rates.10000.Rzoo.RAD = Rates.10000.Rzoo.RAD[order(Rates.10000.Rzoo.RAD$PercSeq),]
#Create the datafarme
stackedRates.10000.Rzoo.RAD = t(data.frame(Rates.10000.Rzoo.RAD$Value[Rates.10000.Rzoo.RAD$Rate == "TN"], Rates.10000.Rzoo.RAD$Value[Rates.10000.Rzoo.RAD$Rate == "TP"],
                                           Rates.10000.Rzoo.RAD$Value[Rates.10000.Rzoo.RAD$Rate == "FN"], Rates.10000.Rzoo.RAD$Value[Rates.10000.Rzoo.RAD$Rate == "FP"]))
row.names(stackedRates.10000.Rzoo.RAD) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.10000.Rzoo.RAD) = PercSeq.10000.Rzoo.RAD

## Read the data TFPN Rates Ne 1'000 RzooROH 3R ARRAY ##

#Read the data
dtaTFPN.10000.Rzoo.50k = read.table("./LargePop/Analyses/SmallArray_RZooRoH_TFPN_Rates_GEN1000.txt", header = T)
dtaTFPN.10000.Rzoo.700k = read.table("./LargePop/Analyses/LargeArray_RZooRoH_TFPN_Rates_GEN1000.txt", header = T)
dtaTFPN.10000.Rzoo.ARRAY = rbind(dtaTFPN.10000.Rzoo.50k, dtaTFPN.10000.Rzoo.700k)
#Get a vector with mean Rates per RAD_Frag
TP.10000.Rzoo.ARRAY = aggregate(dtaTFPN.10000.Rzoo.ARRAY$TP, by = list(dtaTFPN.10000.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
TN.10000.Rzoo.ARRAY = aggregate(dtaTFPN.10000.Rzoo.ARRAY$TN, by = list(dtaTFPN.10000.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
FP.10000.Rzoo.ARRAY = aggregate(dtaTFPN.10000.Rzoo.ARRAY$FP, by = list(dtaTFPN.10000.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
FN.10000.Rzoo.ARRAY = aggregate(dtaTFPN.10000.Rzoo.ARRAY$FN, by = list(dtaTFPN.10000.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one df
Rates.10000.Rzoo.ARRAY = rbind(cbind(TP.10000.Rzoo.ARRAY, Rate=rep("TP", nrow(TP.10000.Rzoo.ARRAY))), cbind(TN.10000.Rzoo.ARRAY, Rate=rep("TN", nrow(TN.10000.Rzoo.ARRAY))),
                               cbind(FP.10000.Rzoo.ARRAY, Rate=rep("FP", nrow(FP.10000.Rzoo.ARRAY))), cbind(FN.10000.Rzoo.ARRAY, Rate=rep("FN", nrow(FN.10000.Rzoo.ARRAY))))
#Change names
colnames(Rates.10000.Rzoo.ARRAY) = c("PercSeq", "Value", "Rate")
Rates.10000.Rzoo.ARRAY$PercSeq = factor(Rates.10000.Rzoo.ARRAY$PercSeq, levels = c("SmallArray", "LargeArray"))
#Order data
Rates.10000.Rzoo.ARRAY = Rates.10000.Rzoo.ARRAY[order(Rates.10000.Rzoo.ARRAY$PercSeq),]
#Create the datafarme
stackedRates.10000.Rzoo.ARRAY = t(data.frame(Rates.10000.Rzoo.ARRAY$Value[Rates.10000.Rzoo.ARRAY$Rate == "TN"], Rates.10000.Rzoo.ARRAY$Value[Rates.10000.Rzoo.ARRAY$Rate == "TP"],
                                             Rates.10000.Rzoo.ARRAY$Value[Rates.10000.Rzoo.ARRAY$Rate == "FN"], Rates.10000.Rzoo.ARRAY$Value[Rates.10000.Rzoo.ARRAY$Rate == "FP"]))
row.names(stackedRates.10000.Rzoo.ARRAY) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.10000.Rzoo.ARRAY) = PercSeq.ARRAY

## Read the data TFPN Rates Ne 1'000 RZooROH WGS ##

PercSeq.WGS = "WGS"

#Read the data
dtaTFPN.10000.Rzoo.WGS = read.table("./LargePop/Analyses/WGS_RZooRoH_TFPN_Rates_GEN1000.txt", header = T)
#Get a vector with mean Rates per RAD_Frag
TP.10000.Rzoo.WGS = aggregate(dtaTFPN.10000.Rzoo.WGS$TP, by = list(dtaTFPN.10000.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
TN.10000.Rzoo.WGS = aggregate(dtaTFPN.10000.Rzoo.WGS$TN, by = list(dtaTFPN.10000.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
FP.10000.Rzoo.WGS = aggregate(dtaTFPN.10000.Rzoo.WGS$FP, by = list(dtaTFPN.10000.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
FN.10000.Rzoo.WGS = aggregate(dtaTFPN.10000.Rzoo.WGS$FN, by = list(dtaTFPN.10000.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one dfWGS
Rates.10000.Rzoo.WGS = rbind(cbind(TP.10000.Rzoo.WGS, Rate=rep("TP", nrow(TP.10000.Rzoo.WGS))), cbind(TN.10000.Rzoo.WGS, Rate=rep("TN", nrow(TN.10000.Rzoo.WGS))),
                             cbind(FP.10000.Rzoo.WGS, Rate=rep("FP", nrow(FP.10000.Rzoo.WGS))), cbind(FN.10000.Rzoo.WGS, Rate=rep("FN", nrow(FN.10000.Rzoo.WGS))))
#Change names
colnames(Rates.10000.Rzoo.WGS) = c("PercSeq", "Value", "Rate")
#Create the datafarme
stackedRates.10000.Rzoo.WGS = t(data.frame(Rates.10000.Rzoo.WGS$Value[Rates.10000.Rzoo.WGS$Rate == "TN"], Rates.10000.Rzoo.WGS$Value[Rates.10000.Rzoo.WGS$Rate == "TP"],
                                           Rates.10000.Rzoo.WGS$Value[Rates.10000.Rzoo.WGS$Rate == "FN"], Rates.10000.Rzoo.WGS$Value[Rates.10000.Rzoo.WGS$Rate == "FP"]))
row.names(stackedRates.10000.Rzoo.WGS) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.10000.Rzoo.WGS) = PercSeq.WGS

#Merge WGS, arrays and RAD
stackedRates.10000.Rzoo = cbind(stackedRates.10000.Rzoo.RAD, stackedRates.10000.Rzoo.ARRAY, stackedRates.10000.Rzoo.WGS)

#rm the rest
rm(dtaTFPN.10000.Rzoo.RAD); rm(dtaTFPN.10000.Rzoo.ARRAY); rm(dtaTFPN.10000.Rzoo.50k); rm(dtaTFPN.10000.Rzoo.700k); rm(dtaTFPN.10000.Rzoo.WGS)
rm(FN.10000.Rzoo.RAD); rm(FP.10000.Rzoo.RAD); rm(TN.10000.Rzoo.RAD); rm(TP.10000.Rzoo.RAD)
rm(FN.10000.Rzoo.ARRAY); rm(FP.10000.Rzoo.ARRAY); rm(TN.10000.Rzoo.ARRAY); rm(TP.10000.Rzoo.ARRAY)
rm(FN.10000.Rzoo.WGS); rm(FP.10000.Rzoo.WGS); rm(TN.10000.Rzoo.WGS); rm(TP.10000.Rzoo.WGS)
rm(Rates.10000.Rzoo.RAD);rm(Rates.10000.Rzoo.ARRAY);rm(Rates.10000.Rzoo.WGS)
rm(stackedRates.10000.Rzoo.ARRAY);rm(stackedRates.10000.Rzoo.RAD);rm(stackedRates.10000.Rzoo.WGS)

########################

pdf("~/PhD/ROH/Which_data/Manuscript/SUBMITTED_VERSION/Review_round_I/FIGURES/SUPPLEMENTARY_FIGURES/FIGURES7.pdf", width = 100, height = 50)

#Layout for having all plots
layout(matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5)), 5, 5, byrow = T), heights = c(.5,1,1,1,.4), widths = c(.5,2,1.5,2,1))

#Par for margin and outer margings
par(mar = c(3.5,3.5,3.5,3.5), oma = c(.2,10,.2,.2))

#### PANNEL A: TFPN Rates Ne 1'000 PLINK 100KB ####

#Barplot
barplot(stackedRates.1000.100KB, col = c("#17A398", "#78E33B", "firebrick4", "#F55014"), space = c(.3,.3,.3,1.5,.3,1.5), yaxt = 'n',
        xaxt = 'n', axes = F, ylab = NA, xlab = NA)

#Add the X axis
axis(1, at = c(seq(0.8, to = (3*(1+0.3)-0.5), length.out = 3), seq((3*(1+0.3)-0.5 + 2.5), to = (3*(1+0.3)-0.5 + 3.8), length.out = 2),3*(1+0.3)-0.5 + 3.8 + 2.5), labels = colnames(stackedRates.1000.100KB), cex.axis = 7, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add thre y axis
axis(2, at = seq(0,10,1)/10, labels = seq(0,10,1)/10, cex.axis = 7, hadj = 1.5, las = 1, lwd.ticks = 5, tck = -0.03)

#Add x axis name
mtext(text = "Fraction of genome", side = 2, cex = 7, line = 30)
#Add y axis name
mtext(text = "RAD-sequencing\n(% of genome sequenced)", side = 1, cex = 5, line = 45, at = 2.1)
#Add y axis name
mtext(text = "SNP ARRAYs", side = 1, cex = 5, line = 45, at = 6.5)

#### PANNEL B: TFPN Rates Ne 10'000 PLINK 100KB RAD ####

#Barplot
barplot(stackedRates.1000.Rzoo, col = c("#17A398", "#78E33B", "firebrick4", "#F55014"), space = c(.3,.3,.3,1.5,.3,1.5), yaxt = 'n',
        xaxt = 'n', axes = F, ylab = NA, xlab = NA)

#Add the X axis
axis(1, at = c(seq(0.8, to = (3*(1+0.3)-0.5), length.out = 3), seq((3*(1+0.3)-0.5 + 2.5), to = (3*(1+0.3)-0.5 + 3.8), length.out = 2),3*(1+0.3)-0.5 + 3.8 + 2.5), labels = colnames(stackedRates.1000.Rzoo), cex.axis = 7, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add thre y axis
axis(2, at = seq(0,10,1)/10, labels = seq(0,10,1)/10, cex.axis = 7, hadj = 1.5, las = 1, lwd.ticks = 5, tck = -0.03)

#Add x axis name
mtext(text = "Fraction of genome", side = 2, cex = 7, line = 30)
#Add y axis name
mtext(text = "RAD-sequencing\n(% of genome sequenced)", side = 1, cex = 5, line = 45, at = 2.1)
#Add y axis name
mtext(text = "SNP ARRAYs", side = 1, cex = 5, line = 45, at = 6.5)

#### PANNEL C: TFPN Rates Ne 1'000 RzooRoH 3R RAD ####

#Barplot
barplot(stackedRates.10000.100KB, col = c("#17A398", "#78E33B", "firebrick4", "#F55014"), space = c(.3,.3,.3,1.5,.3,1.5), yaxt = 'n',
        xaxt = 'n', axes = F, ylab = NA, xlab = NA)

#Add the X axis
axis(1, at = c(seq(0.8, to = (3*(1+0.3)-0.5), length.out = 3), seq((3*(1+0.3)-0.5 + 2.5), to = (3*(1+0.3)-0.5 + 3.8), length.out = 2),3*(1+0.3)-0.5 + 3.8 + 2.5), labels = colnames(stackedRates.10000.100KB), cex.axis = 7, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add thre y axis
axis(2, at = seq(0,10,1)/10, labels = seq(0,10,1)/10, cex.axis = 7, hadj = 1.5, las = 1, lwd.ticks = 5, tck = -0.03)

#Add x axis name
mtext(text = "Fraction of genome", side = 2, cex = 7, line = 30)
#Add y axis name
mtext(text = "RAD-sequencing\n(% of genome sequenced)", side = 1, cex = 5, line = 45, at = 2.1)
#Add y axis name
mtext(text = "SNP ARRAYs", side = 1, cex = 5, line = 45, at = 6.5)

#### PANNEL D: TFPN Rates Ne 10'000 RzooRoH 3R RAD ####

#Barplot
barplot(stackedRates.10000.Rzoo, col = c("#17A398", "#78E33B", "firebrick4", "#F55014"), space = c(.3,.3,.3,1.5,.3,1.5), yaxt = 'n',
        xaxt = 'n', axes = F, ylab = NA, xlab = NA)

#Add the X axis
axis(1, at = c(seq(0.8, to = (3*(1+0.3)-0.5), length.out = 3), seq((3*(1+0.3)-0.5 + 2.5), to = (3*(1+0.3)-0.5 + 3.8), length.out = 2),3*(1+0.3)-0.5 + 3.8 + 2.5), labels = colnames(stackedRates.10000.Rzoo), cex.axis = 7, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add thre y axis
axis(2, at = seq(0,10,1)/10, labels = seq(0,10,1)/10, cex.axis = 7, hadj = 1.5, las = 1, lwd.ticks = 5, tck = -0.03)

#Add x axis name
mtext(text = "Fraction of genome", side = 2, cex = 7, line = 30)
#Add y axis name
mtext(text = "RAD-sequencing\n(% of genome sequenced)", side = 1, cex = 5, line = 45, at = 2.1)
#Add y axis name
mtext(text = "SNP ARRAYs", side = 1, cex = 5, line = 45, at = 6.5)

#CLose the plot
dev.off()


#######################################################################
################## FIGURE S8 PLINK TESTING PARAMETERS #################
#######################################################################

paramset = 12

#FROH

setwd("~/PhD/ROH/Which_data/Data/REVIEWed/")

## SmallPop

#Read the data FROH ne 1'000
dtaFROH.1000.100KB.RAD = read.table(paste0("./SmallPop/TRIAL_PLINK_PARAM/RADSEQ_PLINK_FROH_paramset", paramset, ".txt"), header = T)
#Subsample the 4 PERCENTAGES
dtaFROH.1000.100KB.RAD = dtaFROH.1000.100KB.RAD[dtaFROH.1000.100KB.RAD$WINDOWS.NB %in% c(180000, 240000, 600000),]
#SNPs ARRAYs
dtaFROH.1000.100KB.ARRAY = rbind(read.table("./SmallPop/Analyses/SmallArray_PLINK_FROH.txt", header = T),read.table("./SmallPop/Analyses/LargeArray_PLINK_FROH.txt", header = T))
dtaFROH.1000.100KB.ARRAY$WINDOWS.NB = factor(dtaFROH.1000.100KB.ARRAY$WINDOWS.NB, levels = c("SMALL_ARRAY","LARGE_ARRAY"))
#WGS
dtaFROH.1000.100KB.WGS = read.table("./SmallPop/Analyses/WGS_PLINK_FROH.txt", header = T)

# ROHs DIST

dtaROHSDIST.1000.PLINK.RAD = read.table(paste0("./SmallPop/TRIAL_PLINK_PARAM/RADSEQ_PLINK_ROHsDistribution_paramset", paramset,".txt"), header = T)
#Take only the perc we are interested in
dtaROHSDIST.1000.PLINK.RAD = dtaROHSDIST.1000.PLINK.RAD[dtaROHSDIST.1000.PLINK.RAD$WINDOWS.NB %in% c("180000", "240000", "600000"),]
colnames(dtaROHSDIST.1000.PLINK.RAD) = c("SimuID","SUB_method","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.RAD = aggregate(dtaROHSDIST.1000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.RAD$SUB_method), FUN = mean)
#Calculate sd
ROHsDist.SD.1000.PLINK.RAD = as.data.frame(as.matrix(aggregate(dtaROHSDIST.1000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.RAD$SUB_method), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.RAD = merge(dtaROHSDIST.plot.1000.PLINK.RAD, ROHsDist.SD.1000.PLINK.RAD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.RAD) = c("CLASS", "PERC", "x", "sd")
#Pass sd to numeric
dtaROHSDIST.plot.1000.PLINK.RAD$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.RAD$sd))
#Order the df
dtaROHSDIST.plot.1000.PLINK.RAD = dtaROHSDIST.plot.1000.PLINK.RAD[order(dtaROHSDIST.plot.1000.PLINK.RAD$CLASS, dtaROHSDIST.plot.1000.PLINK.RAD$PERC),]
dtaROHSDIST.plot.1000.PLINK.RAD$PERC = factor(dtaROHSDIST.plot.1000.PLINK.RAD$PERC)
#rm RAD and SD
rm(dtaROHSDIST.1000.PLINK.RAD); rm(ROHsDist.SD.1000.PLINK.RAD)

#Read the ne 1'000 50k array PLINK
dtaROHSDIST.1000.PLINK.50k = read.table("./SmallPop/Analyses/SmallArray_PLINK_ROHsDistributions.txt", header = T)
#Read the ne 1'000 700k array PLINK
dtaROHSDIST.1000.PLINK.700k = read.table("./SmallPop/Analyses/LargeArray_PLINK_ROHsDistributions.txt", header = T)
dtaROHSDIST.1000.PLINK.ARRAY = unique(rbind(dtaROHSDIST.1000.PLINK.50k,dtaROHSDIST.1000.PLINK.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.ARRAY = aggregate(dtaROHSDIST.1000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.ARRAY$SUB_TECH), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.ARRAY = aggregate(dtaROHSDIST.1000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.ARRAY$SUB_TECH, SimuIDs=dtaROHSDIST.1000.PLINK.ARRAY$SimuID, REP=dtaROHSDIST.1000.PLINK.ARRAY$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.ARRAY = merge(dtaROHSDIST.plot.1000.PLINK.ARRAY, ROHsDist.SD.1000.PLINK.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.ARRAY = dtaROHSDIST.plot.1000.PLINK.ARRAY[order(dtaROHSDIST.plot.1000.PLINK.ARRAY$CLASS, dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC),]
dtaROHSDIST.plot.1000.PLINK.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.50k);rm(dtaROHSDIST.1000.PLINK.700k);rm(dtaROHSDIST.1000.PLINK.ARRAY);rm(ROHsDist.MeanforCI.1000.PLINK.ARRAY);rm(ROHsDist.SD.1000.PLINK.ARRAY)

#Read the ne 1'000 WGS PLINK
dtaROHSDIST.1000.PLINK.WGS = read.table("./SmallPop/Analyses//WGS_PLINK_ROHsDistribution.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.WGS = aggregate(dtaROHSDIST.1000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.WGS$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.WGS$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.WGS$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.WGS$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.WGS = aggregate(dtaROHSDIST.1000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.WGS$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.PLINK.WGS$SimID, REP=dtaROHSDIST.1000.PLINK.WGS$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.WGS = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.WGS$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.WGS$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.WGS$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.WGS = merge(dtaROHSDIST.plot.1000.PLINK.WGS, ROHsDist.SD.1000.PLINK.WGS, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.WGS) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.WGS = dtaROHSDIST.plot.1000.PLINK.WGS[order(dtaROHSDIST.plot.1000.PLINK.WGS$CLASS, dtaROHSDIST.plot.1000.PLINK.WGS$PERC),]
dtaROHSDIST.plot.1000.PLINK.WGS$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.WGS$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.WGS);rm(ROHsDist.MeanforCI.1000.PLINK.WGS);rm(ROHsDist.SD.1000.PLINK.WGS)

#Read the ne 1'000 TRUE IBD 100GEN PLINK
dtaROHSDIST.1000.PLINK.TRUEIBD100 = read.table("./SmallPop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD100$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD100$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.PLINK.TRUEIBD100$SimID, REP=dtaROHSDIST.1000.PLINK.TRUEIBD100$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.TRUEIBD100 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100 = merge(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100, ROHsDist.SD.1000.PLINK.TRUEIBD100, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100 = dtaROHSDIST.plot.1000.PLINK.TRUEIBD100[order(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$CLASS, dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC),]
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.TRUEIBD100);rm(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100);rm(ROHsDist.SD.1000.PLINK.TRUEIBD100)
#divide by /1000 to match PLINK (which is in KB)
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$x = dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$x/1000

#Read the ne 1'000 TRUE IBD 1000GEN PLINK
dtaROHSDIST.1000.PLINK.TRUEIBD1000 = read.table("./SmallPop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD1000$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD1000$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.PLINK.TRUEIBD1000$SimID, REP=dtaROHSDIST.1000.PLINK.TRUEIBD1000$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.TRUEIBD1000 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000 = merge(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000, ROHsDist.SD.1000.PLINK.TRUEIBD1000, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000 = dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000[order(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$CLASS, dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC),]
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.TRUEIBD1000);rm(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000);rm(ROHsDist.SD.1000.PLINK.TRUEIBD1000)
#divide by /1000 to match PLINK (which is in KB)
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$x = dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$x/1000

#Merge ALL and remove duplicated WGS lines
dtaROHSDIST.1000.PLINK = unique(rbind(dtaROHSDIST.plot.1000.PLINK.RAD, dtaROHSDIST.plot.1000.PLINK.ARRAY, dtaROHSDIST.plot.1000.PLINK.WGS, dtaROHSDIST.plot.1000.PLINK.TRUEIBD100, dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000))
#Change factor levels
dtaROHSDIST.1000.PLINK$PERC = factor(dtaROHSDIST.1000.PLINK$PERC, levels = c("180000", "240000", "600000", "SMALL_ARRAY", "LARGE_ARRAY", "WGS","TRUE_IBD_100", "TRUE_IBD_1000"))
#Order the df
dtaROHSDIST.1000.PLINK = dtaROHSDIST.1000.PLINK[order(dtaROHSDIST.1000.PLINK$CLASS, dtaROHSDIST.1000.PLINK$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.1000.PLINK.RAD,dtaROHSDIST.plot.1000.PLINK.ARRAY); rm(dtaROHSDIST.plot.1000.PLINK.WGS,dtaROHSDIST.plot.1000.PLINK.TRUEIBD100,dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000)

## LargePop

#Read the data FROH ne 1'000
dtaFROH.10000.100KB.RAD = read.table(paste0("./LargePop/TRIAL_PLINK_PARAM/RADSEQ_PLINK_FROH_paramset", paramset, ".txt"), header = T)
#Subsample the 4 PERCENTAGES
dtaFROH.10000.100KB.RAD = dtaFROH.10000.100KB.RAD[dtaFROH.10000.100KB.RAD$WINDOWS.NB %in% c(15000, 25000, 60000),]
#SNPs ARRAYs
dtaFROH.10000.100KB.ARRAY = rbind(read.table("./LargePop/Analyses/SmallArray_PLINK_FROH.txt", header = T),read.table("./LargePop/Analyses/LargeArray_PLINK_FROH.txt", header = T))
dtaFROH.10000.100KB.ARRAY$WINDOWS.NB = factor(dtaFROH.10000.100KB.ARRAY$WINDOWS.NB, levels = c("SMALL_ARRAY","LARGE_ARRAY"))
#WGS
dtaFROH.10000.100KB.WGS = read.table("./LargePop/Analyses/WGS_PLINK_FROH.txt", header = T)

# ROHs DIST

dtaROHSDIST.10000.PLINK.RAD = read.table(paste0("./LargePop/TRIAL_PLINK_PARAM/RADSEQ_PLINK_ROHsDistribution_paramset", paramset,".txt"), header = T)
#Take only the perc we are interested in
dtaROHSDIST.10000.PLINK.RAD = dtaROHSDIST.10000.PLINK.RAD[dtaROHSDIST.10000.PLINK.RAD$WINDOWS.NB %in% c("15000", "25000", "60000"),]
colnames(dtaROHSDIST.10000.PLINK.RAD) = c("SimuID","SUB_method","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.RAD = aggregate(dtaROHSDIST.10000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.RAD$SUB_method), FUN = mean)
#Calculate sd
ROHsDist.SD.10000.PLINK.RAD = as.data.frame(as.matrix(aggregate(dtaROHSDIST.10000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.RAD$SUB_method), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.RAD = merge(dtaROHSDIST.plot.10000.PLINK.RAD, ROHsDist.SD.10000.PLINK.RAD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.RAD) = c("CLASS", "PERC", "x", "sd")
#Pass sd to numeric
dtaROHSDIST.plot.10000.PLINK.RAD$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.RAD$sd))
#Order the df
dtaROHSDIST.plot.10000.PLINK.RAD = dtaROHSDIST.plot.10000.PLINK.RAD[order(dtaROHSDIST.plot.10000.PLINK.RAD$CLASS, dtaROHSDIST.plot.10000.PLINK.RAD$PERC),]
dtaROHSDIST.plot.10000.PLINK.RAD$PERC = factor(dtaROHSDIST.plot.10000.PLINK.RAD$PERC)
#rm RAD and SD
rm(dtaROHSDIST.10000.PLINK.RAD); rm(ROHsDist.SD.10000.PLINK.RAD)

#Read the ne 1'000 50k array PLINK
dtaROHSDIST.10000.PLINK.50k = read.table("./LargePop/Analyses/SmallArray_PLINK_ROHsDistributions.txt", header = T)
#Read the ne 1'000 700k array PLINK
dtaROHSDIST.10000.PLINK.700k = read.table("./LargePop/Analyses/LargeArray_PLINK_ROHsDistributions.txt", header = T)
dtaROHSDIST.10000.PLINK.ARRAY = unique(rbind(dtaROHSDIST.10000.PLINK.50k,dtaROHSDIST.10000.PLINK.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.ARRAY = aggregate(dtaROHSDIST.10000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.ARRAY$SUB_TECH), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.ARRAY = aggregate(dtaROHSDIST.10000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.ARRAY$SUB_TECH, SimuIDs=dtaROHSDIST.10000.PLINK.ARRAY$SimuID, REP=dtaROHSDIST.10000.PLINK.ARRAY$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.ARRAY = merge(dtaROHSDIST.plot.10000.PLINK.ARRAY, ROHsDist.SD.10000.PLINK.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.ARRAY = dtaROHSDIST.plot.10000.PLINK.ARRAY[order(dtaROHSDIST.plot.10000.PLINK.ARRAY$CLASS, dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC),]
dtaROHSDIST.plot.10000.PLINK.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.50k);rm(dtaROHSDIST.10000.PLINK.700k);rm(dtaROHSDIST.10000.PLINK.ARRAY);rm(ROHsDist.MeanforCI.10000.PLINK.ARRAY);rm(ROHsDist.SD.10000.PLINK.ARRAY)

#Read the ne 1'000 WGS PLINK
dtaROHSDIST.10000.PLINK.WGS = read.table("./LargePop/Analyses//WGS_PLINK_ROHsDistribution.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.WGS = aggregate(dtaROHSDIST.10000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.WGS$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.WGS$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.WGS$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.WGS$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.WGS = aggregate(dtaROHSDIST.10000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.WGS$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.PLINK.WGS$SimID, REP=dtaROHSDIST.10000.PLINK.WGS$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.WGS = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.WGS$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.WGS$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.WGS$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.WGS = merge(dtaROHSDIST.plot.10000.PLINK.WGS, ROHsDist.SD.10000.PLINK.WGS, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.WGS) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.WGS = dtaROHSDIST.plot.10000.PLINK.WGS[order(dtaROHSDIST.plot.10000.PLINK.WGS$CLASS, dtaROHSDIST.plot.10000.PLINK.WGS$PERC),]
dtaROHSDIST.plot.10000.PLINK.WGS$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.WGS$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.WGS);rm(ROHsDist.MeanforCI.10000.PLINK.WGS);rm(ROHsDist.SD.10000.PLINK.WGS)

#Read the ne 1'000 TRUE IBD 100GEN PLINK
dtaROHSDIST.10000.PLINK.TRUEIBD100 = read.table("./LargePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD100$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD100$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.PLINK.TRUEIBD100$SimID, REP=dtaROHSDIST.10000.PLINK.TRUEIBD100$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.TRUEIBD100 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100 = merge(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100, ROHsDist.SD.10000.PLINK.TRUEIBD100, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100 = dtaROHSDIST.plot.10000.PLINK.TRUEIBD100[order(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$CLASS, dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC),]
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.TRUEIBD100);rm(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100);rm(ROHsDist.SD.10000.PLINK.TRUEIBD100)
#divide by /10000 to match PLINK (which is in KB)
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$x = dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$x/1000

#Read the ne 1'000 TRUE IBD 10000GEN PLINK
dtaROHSDIST.10000.PLINK.TRUEIBD10000 = read.table("./LargePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD10000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD10000$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD10000$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.TRUEIBD10000 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD10000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD10000$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD10000$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.PLINK.TRUEIBD10000$SimID, REP=dtaROHSDIST.10000.PLINK.TRUEIBD10000$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.TRUEIBD10000 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD10000$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD10000$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD10000$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000 = merge(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000, ROHsDist.SD.10000.PLINK.TRUEIBD10000, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000 = dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000[order(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$CLASS, dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$PERC),]
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.TRUEIBD10000);rm(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD10000);rm(ROHsDist.SD.10000.PLINK.TRUEIBD10000)
#divide by /10000 to match PLINK (which is in KB)
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$x = dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$x/1000

#Merge ALL and remove duplicated WGS lines
dtaROHSDIST.10000.PLINK = unique(rbind(dtaROHSDIST.plot.10000.PLINK.RAD, dtaROHSDIST.plot.10000.PLINK.ARRAY, dtaROHSDIST.plot.10000.PLINK.WGS, dtaROHSDIST.plot.10000.PLINK.TRUEIBD100, dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000))
#Change factor levels
dtaROHSDIST.10000.PLINK$PERC = factor(dtaROHSDIST.10000.PLINK$PERC, levels = c("15000", "25000", "60000", "SMALL_ARRAY", "LARGE_ARRAY", "WGS","TRUE_IBD_100", "TRUE_IBD_1000"))
#Order the df
dtaROHSDIST.10000.PLINK = dtaROHSDIST.10000.PLINK[order(dtaROHSDIST.10000.PLINK$CLASS, dtaROHSDIST.10000.PLINK$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.10000.PLINK.RAD,dtaROHSDIST.plot.10000.PLINK.ARRAY); rm(dtaROHSDIST.plot.10000.PLINK.WGS,dtaROHSDIST.plot.10000.PLINK.TRUEIBD100,dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000)


add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

namess = c("1_PLINK: WINDOWSize ~ DENS; NO MIN; NO MAX; GAP 200, Het 1","2_PLINK: WINDOWSize ~ DENS; MIN 10; MAX 100; GAP 200, Het 1",
           "3_PLINK: WINDOWSize ~ DENS; MIN 30; MAX 100; GAP 200, Het 1","3_PLINK: WINDOWSize ~ DENS; MIN 20; MAX 100; GAP 200, Het 1",
           "5_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 500, Het 2", "6_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 750, Het 2",
           "7_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; NO MIN; NO MAX, Het 2",
           "8_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; MIN 200; MAX 1000, Het 2",
           "9_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; MIN 500; MAX 1000, Het 2",
           "10_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 1000, Het 2","11_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 1500, Het 2",
           "12_PLINK: WINDOWSize ~ DENS; NO MIN; NO MAX; GAP 1000, Het 2")

name = namess[paramset]


#PLOT
pdf(paste0("~/PhD/ROH/Which_data/Manuscript/SUBMITTED_VERSION/Review_round_I/FIGURES/SUPPLEMENTARY_FIGURES/FIGURES8.pdf"), width = 45, height = 38)

#Layout for having all plots
layout(matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5)), 5, 5, byrow = T), heights = c(.5,1,1.5,1,.5), widths = c(.5,2,1,2,1))

#Par for margin and outer margings
par(mar = c(3.5,3.5,3.5,3.5), oma = c(.2,10,.2,.2))

#### PANEL A: FROH Small Ne ####

#Plot PLINK RADs SMALL POP EMPTY PLOT
plot(dtaFROH.1000.100KB.RAD$FROH_sub ~ dtaFROH.1000.100KB.RAD$Froh_TRUE_IBD_100GEN, pch = 19,
     ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(0,0.9), ylim = c(0,0.9), cex.axis = 4, type = "n", frame = F)

#Colors RAD SMALL NE
RADcolours = c("#82EEFD", "#0492C2", "#016064")
names(RADcolours) = sort(unique(dtaFROH.1000.100KB.RAD$WINDOWS.NB))
#Shape points SMALL NE
POINTSshapes = c(15,16,17)
names(POINTSshapes) = sort(unique(dtaFROH.1000.100KB.RAD$WINDOWS.NB))

#Loop through the RAD subset SMALL NE
for (sub in sort(unique(dtaFROH.1000.100KB.RAD$WINDOWS.NB))) {
  
  #Add the points SMALL NE
  points(dtaFROH.1000.100KB.RAD$FROH_sub[dtaFROH.1000.100KB.RAD$WINDOWS.NB == sub] ~
           dtaFROH.1000.100KB.RAD$Froh_TRUE_IBD_100GEN[dtaFROH.1000.100KB.RAD$WINDOWS.NB == sub],
         col = add.alpha(RADcolours[names(RADcolours) == sub],.5), pch = POINTSshapes[names(POINTSshapes) == sub], cex = 3)
  
}

## ADD AXES
axis(1,xlab = NULL, at=c(0:9)*0.1, padj = 1.5, cex.axis= 5, lwd.ticks = 3, srt = 45, tck = -0.03)
axis(2,xlab = NULL, at=c(0:9)*0.1, hadj = 1.5, cex.axis= 5, lwd.ticks = 3, las = 1, tck = -0.03)

#Add axes titles
mtext(expression(F[IBD]), side = 1, outer = FALSE, line = 15, cex = 4)
mtext(expression(F[ROH]), side = 2, outer = FALSE, line = 15, cex = 4)

#Add perfect line (must be below everything else)
abline(0,1, lwd = 4)

#Add panel lab
mtext(text = "A", side = 3, cex = 10, line = 18, at = -.4)
mtext(text = "PLINK", side = 3, cex = 8, line = 18, at = .3)

#### PANEL B: FROH Large Ne ####

#Plot PLINK RADs SMALL POP EMPTY PLOT
plot(dtaFROH.10000.100KB.RAD$FROH_sub ~ dtaFROH.10000.100KB.RAD$Froh_TRUE_IBD_100GEN, pch = 19,
     ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(0,0.9), ylim = c(0,0.9), cex.axis = 4, type = "n", frame = F)

#Colors RAD SMALL NE
RADcolours = c("#FFBF00", "#FD6A02", "#FDA172")
names(RADcolours) = sort(unique(dtaFROH.10000.100KB.RAD$WINDOWS.NB))
#Shape points SMALL NE
POINTSshapes = c(15,16,17)
names(POINTSshapes) = sort(unique(dtaFROH.10000.100KB.RAD$WINDOWS.NB))

#Loop through the RAD subset SMALL NE
for (sub in sort(unique(dtaFROH.10000.100KB.RAD$WINDOWS.NB))) {
  
  #Add the points SMALL NE
  points(dtaFROH.10000.100KB.RAD$FROH_sub[dtaFROH.10000.100KB.RAD$WINDOWS.NB == sub] ~
           dtaFROH.10000.100KB.RAD$Froh_TRUE_IBD_100GEN[dtaFROH.10000.100KB.RAD$WINDOWS.NB == sub],
         col = add.alpha(RADcolours[names(RADcolours) == sub],.5), pch = POINTSshapes[names(POINTSshapes) == sub], cex = 3)
  
}

## ADD AXES
axis(1,xlab = NULL, at=c(0:9)*0.1, padj = 1.5, cex.axis= 5, lwd.ticks = 3, srt = 45, tck = -0.03)
axis(2,xlab = NULL, at=c(0:9)*0.1, hadj = 1.5, cex.axis= 5, lwd.ticks = 3, las = 1, tck = -0.03)

#Add axes titles
mtext(expression(F[IBD]), side = 1, outer = FALSE, line = 15, cex = 4)
mtext(expression(F[ROH]), side = 2, outer = FALSE, line = 15, cex = 4)

#Add perfect line (must be below everything else)
abline(0,1, lwd = 4)

#Add panel lab
mtext(text = "B", side = 3, cex = 10, line = 18, at = -.4)
mtext(text = "PLINK", side = 3, cex = 8, line = 18, at = .3)

# ROHs Dist

#Define ROHs classes
classesROH = c("< 2", "2 - 4","4 - 6","6 - 10","10 - 16", "> 16")


namess = c("1_PLINK: WINDOWSize ~ DENS; NO MIN; NO MAX; GAP 200, Het 1","2_PLINK: WINDOWSize ~ DENS; MIN 10; MAX 100; GAP 200, Het 1",
           "3_PLINK: WINDOWSize ~ DENS; MIN 30; MAX 100; GAP 200, Het 1","3_PLINK: WINDOWSize ~ DENS; MIN 20; MAX 100; GAP 200, Het 1",
           "5_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 500, Het 2", "6_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 750, Het 2",
           "7_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; NO MIN; NO MAX, Het 2",
           "8_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; MIN 200; MAX 1000, Het 2",
           "9_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; MIN 500; MAX 1000, Het 2",
           "10_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 1000, Het 2","11_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 1500, Het 2",
           "12_PLINK: WINDOWSize ~ DENS; NO MIN; NO MAX; GAP 1000, Het 2")

name = namess[paramset]

#### PANNEL C: ROHsDist SMALL NE ####

#Define colors for subsampled dataset
subCOLORS = c("#82EEFD", "#0492C2", "#016064", "#2832C2", "#0A1172", "grey18")

#Create empty plot
plot(dtaROHSDIST.1000.PLINK$CLASS, dtaROHSDIST.1000.PLINK$x/1000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0,1500))

#Loop through classes to draw boxplots
for(class in sort(unique(dtaROHSDIST.1000.PLINK$CLASS))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.36,class-0.18, class, class+0.18, class+0.36)
  
  #Loop through reduced representation
  for(reduced in 1:3){
    
    #subset the data
    dta_sub = dtaROHSDIST.1000.PLINK[dtaROHSDIST.1000.PLINK$CLASS == class & dtaROHSDIST.1000.PLINK$PERC == levels(dtaROHSDIST.1000.PLINK$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$CLASS == class & dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_100"],
         ytop = dta_sub$x/1000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.09), x1 = (polyXcoord[reduced] + 0.09), y0 = ((dta_sub$x - dta_sub$sd)/1000), y1 = ((dta_sub$x + dta_sub$sd)/1000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.1000.PLINK$CLASS == class],
             y1 = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.1000.PLINK$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add the y axis
axis(2, cex.axis = 4, hadj = 1.5, las = 1, cex.axis = 5, lwd.ticks = 5, tck = -0.03)

#Add axes titles
mtext("Mean Mb ± SD (among individuals)", side = 1, outer = FALSE, line = 15, cex = 4)
mtext("ROHs Length Classes [Mb]", side = 2, outer = FALSE, line = 15, cex = 4)

#Add panel lab
mtext(text = "C", side = 3, cex = 10, line = 18, at = -2)
mtext(text = "PLINK", side = 3, cex = 8, line = 18, at = 3)

#### PANNEL D: ROHsDist LARGE NE ####

#Define colors for subsampled dataset
subCOLORS = c("#FFBF00", "#FD6A02", "#FDA172","#CC7722", "#793802", "grey18")

#Create empty plot
plot(dtaROHSDIST.10000.PLINK$CLASS, dtaROHSDIST.10000.PLINK$x/1000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0,1600))

#Loop through classes to draw boxplots
for(class in sort(unique(dtaROHSDIST.10000.PLINK$CLASS))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.36,class-0.18, class, class+0.18, class+0.36)
  
  #Loop through reduced representation
  for(reduced in 1:3){
    
    #subset the data
    dta_sub = dtaROHSDIST.10000.PLINK[dtaROHSDIST.10000.PLINK$CLASS == class & dtaROHSDIST.10000.PLINK$PERC == levels(dtaROHSDIST.10000.PLINK$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$CLASS == class & dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_100"],
         ytop = dta_sub$x/1000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.09), x1 = (polyXcoord[reduced] + 0.09), y0 = ((dta_sub$x - dta_sub$sd)/1000), y1 = ((dta_sub$x + dta_sub$sd)/1000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.10000.PLINK$CLASS == class],
             y1 = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.10000.PLINK$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add the y axis
axis(2, cex.axis = 4, hadj = 1.5, las = 1, cex.axis = 5, lwd.ticks = 5, tck = -0.03)

#Add axes titles
mtext("Mean Mb ± SD (among individuals)", side = 1, outer = FALSE, line = 15, cex = 4)
mtext("ROHs Length Classes [Mb]", side = 2, outer = FALSE, line = 15, cex = 4)

#Add panel lab
mtext(text = "D", side = 3, cex = 10, line = 18, at = -2)
mtext(text = "PLINK", side = 3, cex = 8, line = 18, at = 3)


dev.off()


#######################################################################
################## FIGURE S9 PLINK TESTING PARAMETERS #################
#######################################################################

paramset = 8

#FROH

setwd("~/PhD/ROH/Which_data/Data/REVIEWed/")

## SmallPop

#Read the data FROH ne 1'000
dtaFROH.1000.100KB.RAD = read.table(paste0("./SmallPop/TRIAL_PLINK_PARAM/RADSEQ_PLINK_FROH_paramset", paramset, ".txt"), header = T)
#Subsample the 4 PERCENTAGES
dtaFROH.1000.100KB.RAD = dtaFROH.1000.100KB.RAD[dtaFROH.1000.100KB.RAD$WINDOWS.NB %in% c(180000, 240000, 600000),]
#SNPs ARRAYs
dtaFROH.1000.100KB.ARRAY = rbind(read.table("./SmallPop/Analyses/SmallArray_PLINK_FROH.txt", header = T),read.table("./SmallPop/Analyses/LargeArray_PLINK_FROH.txt", header = T))
dtaFROH.1000.100KB.ARRAY$WINDOWS.NB = factor(dtaFROH.1000.100KB.ARRAY$WINDOWS.NB, levels = c("SMALL_ARRAY","LARGE_ARRAY"))
#WGS
dtaFROH.1000.100KB.WGS = read.table("./SmallPop/Analyses/WGS_PLINK_FROH.txt", header = T)

# ROHs DIST

dtaROHSDIST.1000.PLINK.RAD = read.table(paste0("./SmallPop/TRIAL_PLINK_PARAM/RADSEQ_PLINK_ROHsDistribution_paramset", paramset,".txt"), header = T)
#Take only the perc we are interested in
dtaROHSDIST.1000.PLINK.RAD = dtaROHSDIST.1000.PLINK.RAD[dtaROHSDIST.1000.PLINK.RAD$WINDOWS.NB %in% c("180000", "240000", "600000"),]
colnames(dtaROHSDIST.1000.PLINK.RAD) = c("SimuID","SUB_method","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.RAD = aggregate(dtaROHSDIST.1000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.RAD$SUB_method), FUN = mean)
#Calculate sd
ROHsDist.SD.1000.PLINK.RAD = as.data.frame(as.matrix(aggregate(dtaROHSDIST.1000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.RAD$SUB_method), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.RAD = merge(dtaROHSDIST.plot.1000.PLINK.RAD, ROHsDist.SD.1000.PLINK.RAD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.RAD) = c("CLASS", "PERC", "x", "sd")
#Pass sd to numeric
dtaROHSDIST.plot.1000.PLINK.RAD$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.RAD$sd))
#Order the df
dtaROHSDIST.plot.1000.PLINK.RAD = dtaROHSDIST.plot.1000.PLINK.RAD[order(dtaROHSDIST.plot.1000.PLINK.RAD$CLASS, dtaROHSDIST.plot.1000.PLINK.RAD$PERC),]
dtaROHSDIST.plot.1000.PLINK.RAD$PERC = factor(dtaROHSDIST.plot.1000.PLINK.RAD$PERC)
#rm RAD and SD
rm(dtaROHSDIST.1000.PLINK.RAD); rm(ROHsDist.SD.1000.PLINK.RAD)

#Read the ne 1'000 50k array PLINK
dtaROHSDIST.1000.PLINK.50k = read.table("./SmallPop/Analyses/SmallArray_PLINK_ROHsDistributions.txt", header = T)
#Read the ne 1'000 700k array PLINK
dtaROHSDIST.1000.PLINK.700k = read.table("./SmallPop/Analyses/LargeArray_PLINK_ROHsDistributions.txt", header = T)
dtaROHSDIST.1000.PLINK.ARRAY = unique(rbind(dtaROHSDIST.1000.PLINK.50k,dtaROHSDIST.1000.PLINK.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.ARRAY = aggregate(dtaROHSDIST.1000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.ARRAY$SUB_TECH), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.ARRAY = aggregate(dtaROHSDIST.1000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.ARRAY$SUB_TECH, SimuIDs=dtaROHSDIST.1000.PLINK.ARRAY$SimuID, REP=dtaROHSDIST.1000.PLINK.ARRAY$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.ARRAY = merge(dtaROHSDIST.plot.1000.PLINK.ARRAY, ROHsDist.SD.1000.PLINK.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.ARRAY = dtaROHSDIST.plot.1000.PLINK.ARRAY[order(dtaROHSDIST.plot.1000.PLINK.ARRAY$CLASS, dtaROHSDIST.plot.1000.PLINK.ARRAY$PERC),]
dtaROHSDIST.plot.1000.PLINK.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.50k);rm(dtaROHSDIST.1000.PLINK.700k);rm(dtaROHSDIST.1000.PLINK.ARRAY);rm(ROHsDist.MeanforCI.1000.PLINK.ARRAY);rm(ROHsDist.SD.1000.PLINK.ARRAY)

#Read the ne 1'000 WGS PLINK
dtaROHSDIST.1000.PLINK.WGS = read.table("./SmallPop/Analyses//WGS_PLINK_ROHsDistribution.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.WGS = aggregate(dtaROHSDIST.1000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.WGS$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.WGS$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.WGS$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.WGS$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.WGS = aggregate(dtaROHSDIST.1000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.WGS$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.PLINK.WGS$SimID, REP=dtaROHSDIST.1000.PLINK.WGS$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.WGS = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.WGS$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.WGS$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.WGS$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.WGS = merge(dtaROHSDIST.plot.1000.PLINK.WGS, ROHsDist.SD.1000.PLINK.WGS, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.WGS) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.WGS = dtaROHSDIST.plot.1000.PLINK.WGS[order(dtaROHSDIST.plot.1000.PLINK.WGS$CLASS, dtaROHSDIST.plot.1000.PLINK.WGS$PERC),]
dtaROHSDIST.plot.1000.PLINK.WGS$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.WGS$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.WGS);rm(ROHsDist.MeanforCI.1000.PLINK.WGS);rm(ROHsDist.SD.1000.PLINK.WGS)

#Read the ne 1'000 TRUE IBD 100GEN PLINK
dtaROHSDIST.1000.PLINK.TRUEIBD100 = read.table("./SmallPop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD100$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD100$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.PLINK.TRUEIBD100$SimID, REP=dtaROHSDIST.1000.PLINK.TRUEIBD100$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.TRUEIBD100 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100 = merge(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100, ROHsDist.SD.1000.PLINK.TRUEIBD100, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100 = dtaROHSDIST.plot.1000.PLINK.TRUEIBD100[order(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$CLASS, dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$PERC),]
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.TRUEIBD100);rm(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD100);rm(ROHsDist.SD.1000.PLINK.TRUEIBD100)
#divide by /1000 to match PLINK (which is in KB)
dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$x = dtaROHSDIST.plot.1000.PLINK.TRUEIBD100$x/1000

#Read the ne 1'000 TRUE IBD 1000GEN PLINK
dtaROHSDIST.1000.PLINK.TRUEIBD1000 = read.table("./SmallPop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD1000$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC = factor(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC), levels = unique(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000 = aggregate(dtaROHSDIST.1000.PLINK.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.1000.PLINK.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.1000.PLINK.TRUEIBD1000$WINDOWS.NB, SimuIDs=dtaROHSDIST.1000.PLINK.TRUEIBD1000$SimID, REP=dtaROHSDIST.1000.PLINK.TRUEIBD1000$REP), FUN = mean)
ROHsDist.SD.1000.PLINK.TRUEIBD1000 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000$x, by = list(CLASS=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000$CLASS, PERC=ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000 = merge(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000, ROHsDist.SD.1000.PLINK.TRUEIBD1000, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000 = dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000[order(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$CLASS, dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$PERC),]
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$sd = as.numeric(as.character(dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.1000.PLINK.TRUEIBD1000);rm(ROHsDist.MeanforCI.1000.PLINK.TRUEIBD1000);rm(ROHsDist.SD.1000.PLINK.TRUEIBD1000)
#divide by /1000 to match PLINK (which is in KB)
dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$x = dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000$x/1000

#Merge ALL and remove duplicated WGS lines
dtaROHSDIST.1000.PLINK = unique(rbind(dtaROHSDIST.plot.1000.PLINK.RAD, dtaROHSDIST.plot.1000.PLINK.ARRAY, dtaROHSDIST.plot.1000.PLINK.WGS, dtaROHSDIST.plot.1000.PLINK.TRUEIBD100, dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000))
#Change factor levels
dtaROHSDIST.1000.PLINK$PERC = factor(dtaROHSDIST.1000.PLINK$PERC, levels = c("180000", "240000", "600000", "SMALL_ARRAY", "LARGE_ARRAY", "WGS","TRUE_IBD_100", "TRUE_IBD_1000"))
#Order the df
dtaROHSDIST.1000.PLINK = dtaROHSDIST.1000.PLINK[order(dtaROHSDIST.1000.PLINK$CLASS, dtaROHSDIST.1000.PLINK$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.1000.PLINK.RAD,dtaROHSDIST.plot.1000.PLINK.ARRAY); rm(dtaROHSDIST.plot.1000.PLINK.WGS,dtaROHSDIST.plot.1000.PLINK.TRUEIBD100,dtaROHSDIST.plot.1000.PLINK.TRUEIBD1000)

## LargePop

#Read the data FROH ne 1'000
dtaFROH.10000.100KB.RAD = read.table(paste0("./LargePop/TRIAL_PLINK_PARAM/RADSEQ_PLINK_FROH_paramset", paramset, ".txt"), header = T)
#Subsample the 4 PERCENTAGES
dtaFROH.10000.100KB.RAD = dtaFROH.10000.100KB.RAD[dtaFROH.10000.100KB.RAD$WINDOWS.NB %in% c(15000, 25000, 60000),]
#SNPs ARRAYs
dtaFROH.10000.100KB.ARRAY = rbind(read.table("./LargePop/Analyses/SmallArray_PLINK_FROH.txt", header = T),read.table("./LargePop/Analyses/LargeArray_PLINK_FROH.txt", header = T))
dtaFROH.10000.100KB.ARRAY$WINDOWS.NB = factor(dtaFROH.10000.100KB.ARRAY$WINDOWS.NB, levels = c("SMALL_ARRAY","LARGE_ARRAY"))
#WGS
dtaFROH.10000.100KB.WGS = read.table("./LargePop/Analyses/WGS_PLINK_FROH.txt", header = T)

# ROHs DIST

dtaROHSDIST.10000.PLINK.RAD = read.table(paste0("./LargePop/TRIAL_PLINK_PARAM/RADSEQ_PLINK_ROHsDistribution_paramset", paramset,".txt"), header = T)
#Take only the perc we are interested in
dtaROHSDIST.10000.PLINK.RAD = dtaROHSDIST.10000.PLINK.RAD[dtaROHSDIST.10000.PLINK.RAD$WINDOWS.NB %in% c("15000", "25000", "60000"),]
colnames(dtaROHSDIST.10000.PLINK.RAD) = c("SimuID","SUB_method","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.RAD = aggregate(dtaROHSDIST.10000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.RAD$SUB_method), FUN = mean)
#Calculate sd
ROHsDist.SD.10000.PLINK.RAD = as.data.frame(as.matrix(aggregate(dtaROHSDIST.10000.PLINK.RAD$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.RAD$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.RAD$SUB_method), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.RAD = merge(dtaROHSDIST.plot.10000.PLINK.RAD, ROHsDist.SD.10000.PLINK.RAD, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.RAD) = c("CLASS", "PERC", "x", "sd")
#Pass sd to numeric
dtaROHSDIST.plot.10000.PLINK.RAD$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.RAD$sd))
#Order the df
dtaROHSDIST.plot.10000.PLINK.RAD = dtaROHSDIST.plot.10000.PLINK.RAD[order(dtaROHSDIST.plot.10000.PLINK.RAD$CLASS, dtaROHSDIST.plot.10000.PLINK.RAD$PERC),]
dtaROHSDIST.plot.10000.PLINK.RAD$PERC = factor(dtaROHSDIST.plot.10000.PLINK.RAD$PERC)
#rm RAD and SD
rm(dtaROHSDIST.10000.PLINK.RAD); rm(ROHsDist.SD.10000.PLINK.RAD)

#Read the ne 1'000 50k array PLINK
dtaROHSDIST.10000.PLINK.50k = read.table("./LargePop/Analyses/SmallArray_PLINK_ROHsDistributions.txt", header = T)
#Read the ne 1'000 700k array PLINK
dtaROHSDIST.10000.PLINK.700k = read.table("./LargePop/Analyses/LargeArray_PLINK_ROHsDistributions.txt", header = T)
dtaROHSDIST.10000.PLINK.ARRAY = unique(rbind(dtaROHSDIST.10000.PLINK.50k,dtaROHSDIST.10000.PLINK.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.ARRAY = aggregate(dtaROHSDIST.10000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.ARRAY$SUB_TECH), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.ARRAY = aggregate(dtaROHSDIST.10000.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.ARRAY$SUB_TECH, SimuIDs=dtaROHSDIST.10000.PLINK.ARRAY$SimuID, REP=dtaROHSDIST.10000.PLINK.ARRAY$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.ARRAY = merge(dtaROHSDIST.plot.10000.PLINK.ARRAY, ROHsDist.SD.10000.PLINK.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.ARRAY = dtaROHSDIST.plot.10000.PLINK.ARRAY[order(dtaROHSDIST.plot.10000.PLINK.ARRAY$CLASS, dtaROHSDIST.plot.10000.PLINK.ARRAY$PERC),]
dtaROHSDIST.plot.10000.PLINK.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.50k);rm(dtaROHSDIST.10000.PLINK.700k);rm(dtaROHSDIST.10000.PLINK.ARRAY);rm(ROHsDist.MeanforCI.10000.PLINK.ARRAY);rm(ROHsDist.SD.10000.PLINK.ARRAY)

#Read the ne 1'000 WGS PLINK
dtaROHSDIST.10000.PLINK.WGS = read.table("./LargePop/Analyses//WGS_PLINK_ROHsDistribution.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.WGS = aggregate(dtaROHSDIST.10000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.WGS$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.WGS$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.WGS$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.WGS$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.WGS = aggregate(dtaROHSDIST.10000.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.WGS$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.PLINK.WGS$SimID, REP=dtaROHSDIST.10000.PLINK.WGS$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.WGS = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.WGS$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.WGS$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.WGS$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.WGS = merge(dtaROHSDIST.plot.10000.PLINK.WGS, ROHsDist.SD.10000.PLINK.WGS, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.WGS) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.WGS = dtaROHSDIST.plot.10000.PLINK.WGS[order(dtaROHSDIST.plot.10000.PLINK.WGS$CLASS, dtaROHSDIST.plot.10000.PLINK.WGS$PERC),]
dtaROHSDIST.plot.10000.PLINK.WGS$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.WGS$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.WGS);rm(ROHsDist.MeanforCI.10000.PLINK.WGS);rm(ROHsDist.SD.10000.PLINK.WGS)

#Read the ne 1'000 TRUE IBD 100GEN PLINK
dtaROHSDIST.10000.PLINK.TRUEIBD100 = read.table("./LargePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD100$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD100$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.PLINK.TRUEIBD100$SimID, REP=dtaROHSDIST.10000.PLINK.TRUEIBD100$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.TRUEIBD100 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100 = merge(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100, ROHsDist.SD.10000.PLINK.TRUEIBD100, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100 = dtaROHSDIST.plot.10000.PLINK.TRUEIBD100[order(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$CLASS, dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$PERC),]
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.TRUEIBD100);rm(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD100);rm(ROHsDist.SD.10000.PLINK.TRUEIBD100)
#divide by /10000 to match PLINK (which is in KB)
dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$x = dtaROHSDIST.plot.10000.PLINK.TRUEIBD100$x/1000

#Read the ne 1'000 TRUE IBD 10000GEN PLINK
dtaROHSDIST.10000.PLINK.TRUEIBD10000 = read.table("./LargePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD10000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD10000$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD10000$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$PERC = factor(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$PERC), levels = unique(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.10000.PLINK.TRUEIBD10000 = aggregate(dtaROHSDIST.10000.PLINK.TRUEIBD10000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.10000.PLINK.TRUEIBD10000$ROHs_CLASS, PERC=dtaROHSDIST.10000.PLINK.TRUEIBD10000$WINDOWS.NB, SimuIDs=dtaROHSDIST.10000.PLINK.TRUEIBD10000$SimID, REP=dtaROHSDIST.10000.PLINK.TRUEIBD10000$REP), FUN = mean)
ROHsDist.SD.10000.PLINK.TRUEIBD10000 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD10000$x, by = list(CLASS=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD10000$CLASS, PERC=ROHsDist.MeanforCI.10000.PLINK.TRUEIBD10000$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000 = merge(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000, ROHsDist.SD.10000.PLINK.TRUEIBD10000, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000 = dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000[order(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$CLASS, dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$PERC),]
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$sd = as.numeric(as.character(dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.10000.PLINK.TRUEIBD10000);rm(ROHsDist.MeanforCI.10000.PLINK.TRUEIBD10000);rm(ROHsDist.SD.10000.PLINK.TRUEIBD10000)
#divide by /10000 to match PLINK (which is in KB)
dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$x = dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000$x/1000

#Merge ALL and remove duplicated WGS lines
dtaROHSDIST.10000.PLINK = unique(rbind(dtaROHSDIST.plot.10000.PLINK.RAD, dtaROHSDIST.plot.10000.PLINK.ARRAY, dtaROHSDIST.plot.10000.PLINK.WGS, dtaROHSDIST.plot.10000.PLINK.TRUEIBD100, dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000))
#Change factor levels
dtaROHSDIST.10000.PLINK$PERC = factor(dtaROHSDIST.10000.PLINK$PERC, levels = c("15000", "25000", "60000", "SMALL_ARRAY", "LARGE_ARRAY", "WGS","TRUE_IBD_100", "TRUE_IBD_1000"))
#Order the df
dtaROHSDIST.10000.PLINK = dtaROHSDIST.10000.PLINK[order(dtaROHSDIST.10000.PLINK$CLASS, dtaROHSDIST.10000.PLINK$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.10000.PLINK.RAD,dtaROHSDIST.plot.10000.PLINK.ARRAY); rm(dtaROHSDIST.plot.10000.PLINK.WGS,dtaROHSDIST.plot.10000.PLINK.TRUEIBD100,dtaROHSDIST.plot.10000.PLINK.TRUEIBD10000)


add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

namess = c("1_PLINK: WINDOWSize ~ DENS; NO MIN; NO MAX; GAP 200, Het 1","2_PLINK: WINDOWSize ~ DENS; MIN 10; MAX 100; GAP 200, Het 1",
           "3_PLINK: WINDOWSize ~ DENS; MIN 30; MAX 100; GAP 200, Het 1","3_PLINK: WINDOWSize ~ DENS; MIN 20; MAX 100; GAP 200, Het 1",
           "5_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 500, Het 2", "6_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 750, Het 2",
           "7_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; NO MIN; NO MAX, Het 2",
           "8_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; MIN 200; MAX 1000, Het 2",
           "9_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; MIN 500; MAX 1000, Het 2",
           "10_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 1000, Het 2","11_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 1500, Het 2",
           "12_PLINK: WINDOWSize ~ DENS; NO MIN; NO MAX; GAP 1000, Het 2")

name = namess[paramset]

#PLOT
pdf(paste0("~/PhD/ROH/Which_data/Manuscript/SUBMITTED_VERSION/Review_round_I/FIGURES/SUPPLEMENTARY_FIGURES/FIGURES9.pdf"), width = 55, height = 38)

#Layout for having all plots
layout(matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5)), 5, 5, byrow = T), heights = c(.5,1,.8,1,.5), widths = c(.5,2,1,2,1))

#Par for margin and outer margings
par(mar = c(3.5,3.5,3.5,3.5), oma = c(.2,10,.2,.2))

#### PANEL A: FROH Small Ne ####

#Plot PLINK RADs SMALL POP EMPTY PLOT
plot(dtaFROH.1000.100KB.RAD$FROH_sub ~ dtaFROH.1000.100KB.RAD$Froh_TRUE_IBD_100GEN, pch = 19,
     ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(0,0.9), ylim = c(0,0.9), cex.axis = 4, type = "n", frame = F)

#Colors RAD SMALL NE
RADcolours = c("#82EEFD", "#0492C2", "#016064")
names(RADcolours) = sort(unique(dtaFROH.1000.100KB.RAD$WINDOWS.NB))
#Shape points SMALL NE
POINTSshapes = c(15,16,17)
names(POINTSshapes) = sort(unique(dtaFROH.1000.100KB.RAD$WINDOWS.NB))

#Loop through the RAD subset SMALL NE
for (sub in sort(unique(dtaFROH.1000.100KB.RAD$WINDOWS.NB))) {
  
  #Add the points SMALL NE
  points(dtaFROH.1000.100KB.RAD$FROH_sub[dtaFROH.1000.100KB.RAD$WINDOWS.NB == sub] ~
           dtaFROH.1000.100KB.RAD$Froh_TRUE_IBD_100GEN[dtaFROH.1000.100KB.RAD$WINDOWS.NB == sub],
         col = add.alpha(RADcolours[names(RADcolours) == sub],.5), pch = POINTSshapes[names(POINTSshapes) == sub], cex = 3)
  
}

## ADD AXES
axis(1,xlab = NULL, at=c(0:9)*0.1, padj = 1.5, cex.axis= 8, lwd.ticks = 3, srt = 45, tck = -0.03, lwd = 3)
axis(2,xlab = NULL, at=c(0:9)*0.1, hadj = 1.5, cex.axis= 8, lwd.ticks = 3, las = 1, tck = -0.03, lwd = 3)

#Add axes titles
mtext(expression(F[IBD]), side = 1, outer = FALSE, line = 20, cex = 8)
mtext(expression(F[ROH]), side = 2, outer = FALSE, line = 20, cex = 8)

#Add perfect line (must be below everything else)
abline(0,1, lwd = 4)

#Add panel lab
mtext(text = "A", side = 3, cex = 10, line = 18, at = -.4)
mtext(text = "PLINK", side = 3, cex = 8, line = 18, at = .3)

#### PANEL B: FROH Large Ne ####

#Plot PLINK RADs SMALL POP EMPTY PLOT
plot(dtaFROH.10000.100KB.RAD$FROH_sub ~ dtaFROH.10000.100KB.RAD$Froh_TRUE_IBD_100GEN, pch = 19,
     ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(0,0.9), ylim = c(0,0.9), type = "n", frame = F)

#Colors RAD SMALL NE
RADcolours = c("#FFBF00", "#FD6A02", "#FDA172")
names(RADcolours) = sort(unique(dtaFROH.10000.100KB.RAD$WINDOWS.NB))
#Shape points SMALL NE
POINTSshapes = c(15,16,17)
names(POINTSshapes) = sort(unique(dtaFROH.10000.100KB.RAD$WINDOWS.NB))

#Loop through the RAD subset SMALL NE
for (sub in sort(unique(dtaFROH.10000.100KB.RAD$WINDOWS.NB))) {
  
  #Add the points SMALL NE
  points(dtaFROH.10000.100KB.RAD$FROH_sub[dtaFROH.10000.100KB.RAD$WINDOWS.NB == sub] ~
           dtaFROH.10000.100KB.RAD$Froh_TRUE_IBD_100GEN[dtaFROH.10000.100KB.RAD$WINDOWS.NB == sub],
         col = add.alpha(RADcolours[names(RADcolours) == sub],.5), pch = POINTSshapes[names(POINTSshapes) == sub], cex = 3)
  
}

## ADD AXES
axis(1,xlab = NULL, at=c(0:9)*0.1, padj = 1.5, cex.axis= 8, lwd.ticks = 3, srt = 45, tck = -0.03, lwd = 3)
axis(2,xlab = NULL, at=c(0:9)*0.1, hadj = 1.5, cex.axis= 8, lwd.ticks = 3, las = 1, tck = -0.03, lwd = 3)

#Add axes titles
mtext(expression(F[IBD]), side = 1, outer = FALSE, line = 20, cex = 8)
mtext(expression(F[ROH]), side = 2, outer = FALSE, line = 20, cex = 8)

#Add perfect line (must be below everything else)
abline(0,1, lwd = 4)

#Add panel lab
mtext(text = "B", side = 3, cex = 10, line = 18, at = -.4)
mtext(text = "PLINK", side = 3, cex = 8, line = 18, at = .3)

# ROHs Dist

#Define ROHs classes
classesROH = c("< 2", "2 - 4","4 - 6","6 - 10","10 - 16", "> 16")


namess = c("1_PLINK: WINDOWSize ~ DENS; NO MIN; NO MAX; GAP 200, Het 1","2_PLINK: WINDOWSize ~ DENS; MIN 10; MAX 100; GAP 200, Het 1",
           "3_PLINK: WINDOWSize ~ DENS; MIN 30; MAX 100; GAP 200, Het 1","3_PLINK: WINDOWSize ~ DENS; MIN 20; MAX 100; GAP 200, Het 1",
           "5_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 500, Het 2", "6_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 750, Het 2",
           "7_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; NO MIN; NO MAX, Het 2",
           "8_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; MIN 200; MAX 1000, Het 2",
           "9_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP ~ DENSITY; MIN 500; MAX 1000, Het 2",
           "10_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 1000, Het 2","11_PLINK: WINDOWSize ~ DENS: MIN 30; MAX 100; GAP 1500, Het 2",
           "12_PLINK: WINDOWSize ~ DENS; NO MIN; NO MAX; GAP 1000, Het 2")

name = namess[paramset]

#### PANNEL C: ROHsDist SMALL NE ####

#Define colors for subsampled dataset
subCOLORS = c("#82EEFD", "#0492C2", "#016064", "#2832C2", "#0A1172", "grey18")

#Create empty plot
plot(dtaROHSDIST.1000.PLINK$CLASS, dtaROHSDIST.1000.PLINK$x/1000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0,600))

#Loop through classes to draw boxplots
for(class in sort(unique(dtaROHSDIST.1000.PLINK$CLASS))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.36,class-0.18, class, class+0.18, class+0.36)
  
  #Loop through reduced representation
  for(reduced in 1:3){
    
    #subset the data
    dta_sub = dtaROHSDIST.1000.PLINK[dtaROHSDIST.1000.PLINK$CLASS == class & dtaROHSDIST.1000.PLINK$PERC == levels(dtaROHSDIST.1000.PLINK$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$CLASS == class & dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_100"],
         ytop = dta_sub$x/1000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.09), x1 = (polyXcoord[reduced] + 0.09), y0 = ((dta_sub$x - dta_sub$sd)/1000), y1 = ((dta_sub$x + dta_sub$sd)/1000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.1000.PLINK$CLASS == class],
             y1 = (dtaROHSDIST.1000.PLINK$x/1000)[dtaROHSDIST.1000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.1000.PLINK$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03, lwd = 3)
#Add the y axis
axis(2, cex.axis = 4, hadj = 1.5, las = 1, cex.axis = 6, lwd.ticks = 5, tck = -0.03, lwd = 3)

#Add axes titles
mtext("Mean Mb ± SD (among individuals)", side = 1, outer = FALSE, line = 20, cex = 6)
mtext("ROHs Length Classes [Mb]", side = 2, outer = FALSE, line = 20, cex = 6)

#Add panel lab
mtext(text = "C", side = 3, cex = 10, line = 18, at = -2)
mtext(text = "PLINK", side = 3, cex = 8, line = 18, at = 3)

#### PANNEL D: ROHsDist LARGE NE ####

#Define colors for subsampled dataset
subCOLORS = c("#FFBF00", "#FD6A02", "#FDA172","#CC7722", "#793802", "grey18")

#Create empty plot
plot(dtaROHSDIST.10000.PLINK$CLASS, dtaROHSDIST.10000.PLINK$x/1000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0,600))

#Loop through classes to draw boxplots
for(class in sort(unique(dtaROHSDIST.10000.PLINK$CLASS))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.36,class-0.18, class, class+0.18, class+0.36)
  
  #Loop through reduced representation
  for(reduced in 1:3){
    
    #subset the data
    dta_sub = dtaROHSDIST.10000.PLINK[dtaROHSDIST.10000.PLINK$CLASS == class & dtaROHSDIST.10000.PLINK$PERC == levels(dtaROHSDIST.10000.PLINK$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$CLASS == class & dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_100"],
         ytop = dta_sub$x/1000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.09), x1 = (polyXcoord[reduced] + 0.09), y0 = ((dta_sub$x - dta_sub$sd)/1000), y1 = ((dta_sub$x + dta_sub$sd)/1000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.10000.PLINK$CLASS == class],
             y1 = (dtaROHSDIST.10000.PLINK$x/1000)[dtaROHSDIST.10000.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.10000.PLINK$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03, lwd = 3)
#Add the y axis
axis(2, cex.axis = 4, hadj = 1.5, las = 1, cex.axis = 6, lwd.ticks = 5, tck = -0.03, lwd = 3)

#Add axes titles
mtext("Mean Mb ± SD (among individuals)", side = 1, outer = FALSE, line = 20, cex = 6)
mtext("ROHs Length Classes [Mb]", side = 2, outer = FALSE, line = 20, cex = 6)

#Add panel lab
mtext(text = "D", side = 3, cex = 10, line = 18, at = -2)
mtext(text = "PLINK", side = 3, cex = 8, line = 18, at = 3)


dev.off()

###################################################
############ FIGURE S10 Cattle SIMULATIONS ########
###################################################

##Read the data Ne 1'000 50K
dtaFROH.Cat.100KB.50K = read.table("./CattlePop/Analyses/SmallArray_PLINK_FROH.txt", header = T)
colnames(dtaFROH.Cat.100KB.50K)[c(8,9)] = c("Froh_TRUE_IBD_100GEN","Froh_TRUE_IBD_1000GEN")
dtaFROH.Cat.100KB.50K = cbind(dtaFROH.Cat.100KB.50K, ARRAY=rep("SMALL", nrow(dtaFROH.Cat.100KB.50K)))
#Read the data Ne 1'000 700K
dtaFROH.Cat.100KB.700K = read.table("./CattlePop/Analyses/LargeArray_PLINK_FROH.txt", header = T)
dtaFROH.Cat.100KB.700K = cbind(dtaFROH.Cat.100KB.700K, ARRAY=rep("LARGE", nrow(dtaFROH.Cat.100KB.700K)))
dtaFROH.Cat.100KB.ARRAY = rbind(dtaFROH.Cat.100KB.50K,dtaFROH.Cat.100KB.700K)
rm(dtaFROH.Cat.100KB.50K);rm(dtaFROH.Cat.100KB.700K)
dtaFROH.Cat.100KB.ARRAY$ARRAY = factor(dtaFROH.Cat.100KB.ARRAY$ARRAY, levels = c("SMALL","LARGE"))

#Read the data Ne 1'000 50K
dtaFROH.Cat.Rzoo.50K = read.table("./CattlePop/Analyses/SmallArray_RZooRoH_FROH.txt", header = T)
colnames(dtaFROH.Cat.Rzoo.50K) = c("Simu_ID","individuals", "ARRAY", "Froh_array", "Replicate", "Froh_trueIBD100Gen", "Froh_trueIBD1000Gen")
#change model/array into "SMALL"
dtaFROH.Cat.Rzoo.50K$ARRAY = "SMALL"
#Read the data Ne 1'000 700K
dtaFROH.Cat.Rzoo.700K = read.table("./CattlePop/Analyses/LargeArray_RZooRoH_FROH.txt", header = T)
colnames(dtaFROH.Cat.Rzoo.700K) =  c("Simu_ID","individuals", "ARRAY", "Froh_array", "Replicate", "Froh_trueIBD100Gen", "Froh_trueIBD1000Gen")
#change model/array into "SMALL"
dtaFROH.Cat.Rzoo.700K$ARRAY = "LARGE"
dtaFROH.Cat.Rzoo.ARRAY = rbind(dtaFROH.Cat.Rzoo.50K,dtaFROH.Cat.Rzoo.700K)
rm(dtaFROH.Cat.Rzoo.50K);rm(dtaFROH.Cat.Rzoo.700K)
dtaFROH.Cat.Rzoo.ARRAY$ARRAY = factor(dtaFROH.Cat.Rzoo.ARRAY$ARRAY, levels = c("SMALL", "LARGE"))

#Read the data Ne 1'000 WGS PLINK
dtaFROH.Cat.100KB.WGS = read.table("./CattlePop/Analyses/WGS_PLINK_FROH_newparamset5.txt", header = T)
#Read the data Ne 10'000 WGS RZooRoH
dtaFROH.Cat.RZoo.WGS = read.table("./CattlePop/Analyses/WGS_RZooRoH_FROH.txt", header = T)

## ROHs DIST

#Read the ne 1'000 50k array PLINK
dtaROHSDIST.Cat.PLINK.50k = read.table("./CattlePop/Analyses/SmallArray_PLINK_ROHsDistributions.txt", header = T)
#Read the ne 1'000 700k array PLINK
dtaROHSDIST.Cat.PLINK.700k = read.table("./CattlePop/Analyses/LargeArray_PLINK_ROHsDistributions.txt", header = T)
dtaROHSDIST.Cat.PLINK.ARRAY = unique(rbind(dtaROHSDIST.Cat.PLINK.50k,dtaROHSDIST.Cat.PLINK.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.Cat.PLINK.ARRAY = aggregate(dtaROHSDIST.Cat.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.Cat.PLINK.ARRAY$SUB_TECH), FUN = mean)
dtaROHSDIST.plot.Cat.PLINK.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.Cat.PLINK.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.Cat.PLINK.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.Cat.PLINK.ARRAY = aggregate(dtaROHSDIST.Cat.PLINK.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.PLINK.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.Cat.PLINK.ARRAY$SUB_TECH, SimuIDs=dtaROHSDIST.Cat.PLINK.ARRAY$SimuID, REP=dtaROHSDIST.Cat.PLINK.ARRAY$REP), FUN = mean)
ROHsDist.SD.Cat.PLINK.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.Cat.PLINK.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.Cat.PLINK.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.Cat.PLINK.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.Cat.PLINK.ARRAY = merge(dtaROHSDIST.plot.Cat.PLINK.ARRAY, ROHsDist.SD.Cat.PLINK.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.Cat.PLINK.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.Cat.PLINK.ARRAY = dtaROHSDIST.plot.Cat.PLINK.ARRAY[order(dtaROHSDIST.plot.Cat.PLINK.ARRAY$CLASS, dtaROHSDIST.plot.Cat.PLINK.ARRAY$PERC),]
dtaROHSDIST.plot.Cat.PLINK.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.Cat.PLINK.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.Cat.PLINK.50k);rm(dtaROHSDIST.Cat.PLINK.700k);rm(dtaROHSDIST.Cat.PLINK.ARRAY);rm(ROHsDist.MeanforCI.Cat.PLINK.ARRAY);rm(ROHsDist.SD.Cat.PLINK.ARRAY)

#Read the ne 1'000 WGS PLINK
dtaROHSDIST.Cat.PLINK.WGS = read.table("./CattlePop/Analyses//WGS_PLINK_ROHsDistribution_newparamset5.txt", header = F)
colnames(dtaROHSDIST.Cat.PLINK.WGS) = c("SimuID","SUB_TECH","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.Cat.PLINK.WGS = aggregate(dtaROHSDIST.Cat.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.PLINK.WGS$ROHs_CLASS), FUN = mean)
dtaROHSDIST.plot.Cat.PLINK.WGS = as.data.frame(cbind(CLASS=dtaROHSDIST.plot.Cat.PLINK.WGS$CLASS,PERC=rep("WGS",nrow(dtaROHSDIST.plot.Cat.PLINK.WGS)),x = dtaROHSDIST.plot.Cat.PLINK.WGS$x))
dtaROHSDIST.plot.Cat.PLINK.WGS$PERC = factor(as.character(dtaROHSDIST.plot.Cat.PLINK.WGS$PERC), levels = unique(dtaROHSDIST.plot.Cat.PLINK.WGS$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.Cat.PLINK.WGS = aggregate(dtaROHSDIST.Cat.PLINK.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.PLINK.WGS$ROHs_CLASS, PERC=dtaROHSDIST.Cat.PLINK.WGS$SUB_TECH), FUN = mean)
ROHsDist.SD.Cat.PLINK.WGS = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.Cat.PLINK.WGS$x, by = list(CLASS=ROHsDist.MeanforCI.Cat.PLINK.WGS$CLASS, PERC=ROHsDist.MeanforCI.Cat.PLINK.WGS$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.Cat.PLINK.WGS = merge(dtaROHSDIST.plot.Cat.PLINK.WGS, ROHsDist.SD.Cat.PLINK.WGS, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.Cat.PLINK.WGS) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.Cat.PLINK.WGS = dtaROHSDIST.plot.Cat.PLINK.WGS[order(dtaROHSDIST.plot.Cat.PLINK.WGS$CLASS, dtaROHSDIST.plot.Cat.PLINK.WGS$PERC),]
dtaROHSDIST.plot.Cat.PLINK.WGS$sd = as.numeric(as.character(dtaROHSDIST.plot.Cat.PLINK.WGS$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.Cat.PLINK.WGS);rm(ROHsDist.MeanforCI.Cat.PLINK.WGS);rm(ROHsDist.SD.Cat.PLINK.WGS)

#Read the ne 1'000 TRUE IBD 100GEN PLINK
dtaROHSDIST.Cat.PLINK.TRUEIBD100 = read.table("./CattlePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.Cat.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.Cat.PLINK.TRUEIBD100$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100$PERC = factor(as.character(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100$PERC), levels = unique(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.Cat.PLINK.TRUEIBD100 = aggregate(dtaROHSDIST.Cat.PLINK.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.PLINK.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.Cat.PLINK.TRUEIBD100$WINDOWS.NB, SimuIDs=dtaROHSDIST.Cat.PLINK.TRUEIBD100$SimID, REP=dtaROHSDIST.Cat.PLINK.TRUEIBD100$REP), FUN = mean)
ROHsDist.SD.Cat.PLINK.TRUEIBD100 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.Cat.PLINK.TRUEIBD100$x, by = list(CLASS=ROHsDist.MeanforCI.Cat.PLINK.TRUEIBD100$CLASS, PERC=ROHsDist.MeanforCI.Cat.PLINK.TRUEIBD100$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100 = merge(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100, ROHsDist.SD.Cat.PLINK.TRUEIBD100, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100 = dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100[order(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100$CLASS, dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100$PERC),]
dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100$sd = as.numeric(as.character(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100$sd))
#rm 50k, 700k and others
rm(ROHsDist.MeanforCI.Cat.PLINK.TRUEIBD100);rm(ROHsDist.SD.Cat.PLINK.TRUEIBD100)

#Read the ne 1'000 TRUE IBD CatGEN PLINK
dtaROHSDIST.Cat.PLINK.TRUEIBD1000 = read.table("./CattlePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt", header = F)
colnames(dtaROHSDIST.Cat.PLINK.TRUEIBD1000) = colnames(dtaROHSDIST.Cat.PLINK.TRUEIBD100)
#Calculate the mean per individual
dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000 = aggregate(dtaROHSDIST.Cat.PLINK.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.PLINK.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.Cat.PLINK.TRUEIBD1000$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000$PERC = factor(as.character(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000$PERC), levels = unique(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.Cat.PLINK.TRUEIBD1000 = aggregate(dtaROHSDIST.Cat.PLINK.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.PLINK.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.Cat.PLINK.TRUEIBD1000$WINDOWS.NB, SimuIDs=dtaROHSDIST.Cat.PLINK.TRUEIBD1000$SimID, REP=dtaROHSDIST.Cat.PLINK.TRUEIBD1000$REP), FUN = mean)
ROHsDist.SD.Cat.PLINK.TRUEIBD1000 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.Cat.PLINK.TRUEIBD1000$x, by = list(CLASS=ROHsDist.MeanforCI.Cat.PLINK.TRUEIBD1000$CLASS, PERC=ROHsDist.MeanforCI.Cat.PLINK.TRUEIBD1000$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000 = merge(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000, ROHsDist.SD.Cat.PLINK.TRUEIBD1000, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000 = dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000[order(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000$CLASS, dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000$PERC),]
dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000$sd = as.numeric(as.character(dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.Cat.PLINK.TRUEIBD100);rm(dtaROHSDIST.Cat.PLINK.TRUEIBD1000);rm(ROHsDist.MeanforCI.Cat.PLINK.TRUEIBD1000);rm(ROHsDist.SD.Cat.PLINK.TRUEIBD1000)

#Merge ALL and remove duplicated WGS lines
dtaROHSDIST.Cat.PLINK = unique(rbind(dtaROHSDIST.plot.Cat.PLINK.ARRAY, dtaROHSDIST.plot.Cat.PLINK.WGS, dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100, dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000))
#Change factor levels
dtaROHSDIST.Cat.PLINK$PERC = factor(dtaROHSDIST.Cat.PLINK$PERC, levels = c("SMALL_ARRAY","LARGE_ARRAY","WGS","TRUE_IBD_100", "TRUE_IBD_1000"))
#Order the df
dtaROHSDIST.Cat.PLINK = dtaROHSDIST.Cat.PLINK[order(dtaROHSDIST.Cat.PLINK$CLASS, dtaROHSDIST.Cat.PLINK$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.Cat.PLINK.ARRAY,dtaROHSDIST.plot.Cat.PLINK.WGS,dtaROHSDIST.plot.Cat.PLINK.TRUEIBD100,dtaROHSDIST.plot.Cat.PLINK.TRUEIBD1000)
#Divide TRUE IBD by Cat (because in basepair while PLINK in kb)
dtaROHSDIST.Cat.PLINK$x[dtaROHSDIST.Cat.PLINK$PERC %in% c("TRUE_IBD_100", "TRUE_IBD_1000")] = as.numeric(as.character(dtaROHSDIST.Cat.PLINK$x[dtaROHSDIST.Cat.PLINK$PERC %in% c("TRUE_IBD_100", "TRUE_IBD_1000")]))/1000
dtaROHSDIST.Cat.PLINK$x = as.numeric(as.character(dtaROHSDIST.Cat.PLINK$x))
## RZOO

#Read the ne 1'000 50k array Rzoo
dtaROHSDIST.Cat.Rzoo.50k = read.table("./CattlePop/Analyses/SmallArray_RZooRoH_ROHsDistributions.txt", header = T)
#Read the ne 1'000 700k array Rzoo
dtaROHSDIST.Cat.Rzoo.700k = read.table("./CattlePop/Analyses/LargeArray_RZooRoH_ROHsDistributions.txt", header = F)
colnames(dtaROHSDIST.Cat.Rzoo.700k) = colnames(dtaROHSDIST.Cat.Rzoo.50k)
dtaROHSDIST.Cat.Rzoo.ARRAY = unique(rbind(dtaROHSDIST.Cat.Rzoo.50k,dtaROHSDIST.Cat.Rzoo.700k))
#Calculate the mean per individual
dtaROHSDIST.plot.Cat.Rzoo.ARRAY = aggregate(dtaROHSDIST.Cat.Rzoo.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.Rzoo.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.Cat.Rzoo.ARRAY$SUB_TECH), FUN = mean)
dtaROHSDIST.plot.Cat.Rzoo.ARRAY$PERC = factor(as.character(dtaROHSDIST.plot.Cat.Rzoo.ARRAY$PERC), levels = unique(dtaROHSDIST.plot.Cat.Rzoo.ARRAY$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.Cat.Rzoo.ARRAY = aggregate(dtaROHSDIST.Cat.Rzoo.ARRAY$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.Rzoo.ARRAY$ROHs_CLASS, PERC=dtaROHSDIST.Cat.Rzoo.ARRAY$SUB_TECH, SimuIDs=dtaROHSDIST.Cat.Rzoo.ARRAY$SimuID, REP=dtaROHSDIST.Cat.Rzoo.ARRAY$REP), FUN = mean)
ROHsDist.SD.Cat.Rzoo.ARRAY = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.Cat.Rzoo.ARRAY$x, by = list(CLASS=ROHsDist.MeanforCI.Cat.Rzoo.ARRAY$CLASS, PERC=ROHsDist.MeanforCI.Cat.Rzoo.ARRAY$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.Cat.Rzoo.ARRAY = merge(dtaROHSDIST.plot.Cat.Rzoo.ARRAY, ROHsDist.SD.Cat.Rzoo.ARRAY, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.Cat.Rzoo.ARRAY) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.Cat.Rzoo.ARRAY = dtaROHSDIST.plot.Cat.Rzoo.ARRAY[order(dtaROHSDIST.plot.Cat.Rzoo.ARRAY$CLASS, dtaROHSDIST.plot.Cat.Rzoo.ARRAY$PERC),]
dtaROHSDIST.plot.Cat.Rzoo.ARRAY$sd = as.numeric(as.character(dtaROHSDIST.plot.Cat.Rzoo.ARRAY$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.Cat.Rzoo.50k);rm(dtaROHSDIST.Cat.Rzoo.700k);rm(dtaROHSDIST.Cat.Rzoo.ARRAY);rm(ROHsDist.MeanforCI.Cat.Rzoo.ARRAY);rm(ROHsDist.SD.Cat.Rzoo.ARRAY)

#Read the ne 1'000 WGS Rzoo
dtaROHSDIST.Cat.Rzoo.WGS = read.table("./CattlePop/Analyses/WGS_RZooRoH_ROHsDistribution.txt", header = F)
colnames(dtaROHSDIST.Cat.Rzoo.WGS) = c("SimuID","SUB_TECH","REP","ROHs_CLASS","Mean_NROH_per_ind","Mean_Tot_Kb_per_ind")
#Calculate the mean per individual
dtaROHSDIST.plot.Cat.Rzoo.WGS = aggregate(dtaROHSDIST.Cat.Rzoo.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.Rzoo.WGS$ROHs_CLASS), FUN = mean)
dtaROHSDIST.plot.Cat.Rzoo.WGS = as.data.frame(cbind(CLASS=dtaROHSDIST.plot.Cat.Rzoo.WGS$CLASS,PERC=rep("WGS",nrow(dtaROHSDIST.plot.Cat.Rzoo.WGS)),x = dtaROHSDIST.plot.Cat.Rzoo.WGS$x))
dtaROHSDIST.plot.Cat.Rzoo.WGS$PERC = factor(as.character(dtaROHSDIST.plot.Cat.Rzoo.WGS$PERC), levels = unique(dtaROHSDIST.plot.Cat.Rzoo.WGS$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.Cat.Rzoo.WGS = aggregate(dtaROHSDIST.Cat.Rzoo.WGS$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.Rzoo.WGS$ROHs_CLASS, PERC=dtaROHSDIST.Cat.Rzoo.WGS$SUB_TECH), FUN = mean)
ROHsDist.SD.Cat.Rzoo.WGS = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.Cat.Rzoo.WGS$x, by = list(CLASS=ROHsDist.MeanforCI.Cat.Rzoo.WGS$CLASS, PERC=ROHsDist.MeanforCI.Cat.Rzoo.WGS$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.Cat.Rzoo.WGS = merge(dtaROHSDIST.plot.Cat.Rzoo.WGS, ROHsDist.SD.Cat.Rzoo.WGS, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.Cat.Rzoo.WGS) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.Cat.Rzoo.WGS = dtaROHSDIST.plot.Cat.Rzoo.WGS[order(dtaROHSDIST.plot.Cat.Rzoo.WGS$CLASS, dtaROHSDIST.plot.Cat.Rzoo.WGS$PERC),]
dtaROHSDIST.plot.Cat.Rzoo.WGS$sd = as.numeric(as.character(dtaROHSDIST.plot.Cat.Rzoo.WGS$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.Cat.Rzoo.WGS);rm(ROHsDist.MeanforCI.Cat.Rzoo.WGS);rm(ROHsDist.SD.Cat.Rzoo.WGS)

#Read the ne 1'000 TRUE IBD 100GEN Rzoo
dtaROHSDIST.Cat.Rzoo.TRUEIBD100 = read.table("./CattlePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_100GEN.txt", header = T)
#Calculate the mean per individual
dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100 = aggregate(dtaROHSDIST.Cat.Rzoo.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.Rzoo.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.Cat.Rzoo.TRUEIBD100$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100$PERC = factor(as.character(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100$PERC), levels = unique(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.Cat.Rzoo.TRUEIBD100 = aggregate(dtaROHSDIST.Cat.Rzoo.TRUEIBD100$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.Rzoo.TRUEIBD100$ROHs_CLASS, PERC=dtaROHSDIST.Cat.Rzoo.TRUEIBD100$WINDOWS.NB, SimuIDs=dtaROHSDIST.Cat.Rzoo.TRUEIBD100$SimID, REP=dtaROHSDIST.Cat.Rzoo.TRUEIBD100$REP), FUN = mean)
ROHsDist.SD.Cat.Rzoo.TRUEIBD100 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.Cat.Rzoo.TRUEIBD100$x, by = list(CLASS=ROHsDist.MeanforCI.Cat.Rzoo.TRUEIBD100$CLASS, PERC=ROHsDist.MeanforCI.Cat.Rzoo.TRUEIBD100$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100 = merge(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100, ROHsDist.SD.Cat.Rzoo.TRUEIBD100, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100 = dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100[order(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100$CLASS, dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100$PERC),]
dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100$sd = as.numeric(as.character(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100$sd))
#rm 50k, 700k and others
rm(ROHsDist.MeanforCI.Cat.Rzoo.TRUEIBD100);rm(ROHsDist.SD.Cat.Rzoo.TRUEIBD100)

#Read the ne 1'000 TRUE IBD CatGEN Rzoo
dtaROHSDIST.Cat.Rzoo.TRUEIBD1000 = read.table("./CattlePop/Analyses/TRUE_IBD_SEGMENTS_DISTRIBUTION_1000GEN.txt", header = F)
colnames(dtaROHSDIST.Cat.Rzoo.TRUEIBD1000) = colnames(dtaROHSDIST.Cat.Rzoo.TRUEIBD100)
#Calculate the mean per individual
dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000 = aggregate(dtaROHSDIST.Cat.Rzoo.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.Rzoo.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.Cat.Rzoo.TRUEIBD1000$WINDOWS.NB), FUN = mean)
dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000$PERC = factor(as.character(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000$PERC), levels = unique(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000$PERC))
#Calculate mean per replicate per Simu for CI
ROHsDist.MeanforCI.Cat.Rzoo.TRUEIBD1000 = aggregate(dtaROHSDIST.Cat.Rzoo.TRUEIBD1000$Mean_Tot_Kb_per_ind, by = list(CLASS=dtaROHSDIST.Cat.Rzoo.TRUEIBD1000$ROHs_CLASS, PERC=dtaROHSDIST.Cat.Rzoo.TRUEIBD1000$WINDOWS.NB, SimuIDs=dtaROHSDIST.Cat.Rzoo.TRUEIBD1000$SimID, REP=dtaROHSDIST.Cat.Rzoo.TRUEIBD1000$REP), FUN = mean)
ROHsDist.SD.Cat.Rzoo.TRUEIBD1000 = as.data.frame(as.matrix(aggregate(ROHsDist.MeanforCI.Cat.Rzoo.TRUEIBD1000$x, by = list(CLASS=ROHsDist.MeanforCI.Cat.Rzoo.TRUEIBD1000$CLASS, PERC=ROHsDist.MeanforCI.Cat.Rzoo.TRUEIBD1000$PERC), FUN = sd)))
#Merge CI and ROHDist df
dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000 = merge(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000, ROHsDist.SD.Cat.Rzoo.TRUEIBD1000, by = c("CLASS", "PERC"))
colnames(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000) = c("CLASS", "PERC", "x", "sd")
#Order the df
dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000 = dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000[order(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000$CLASS, dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000$PERC),]
dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000$sd = as.numeric(as.character(dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000$sd))
#rm 50k, 700k and others
rm(dtaROHSDIST.Cat.Rzoo.TRUEIBD100);rm(dtaROHSDIST.Cat.Rzoo.TRUEIBD1000);rm(ROHsDist.MeanforCI.Cat.Rzoo.TRUEIBD1000);rm(ROHsDist.SD.Cat.Rzoo.TRUEIBD1000)

dtaROHSDIST.Cat.Rzoo = unique(rbind(dtaROHSDIST.plot.Cat.Rzoo.ARRAY, dtaROHSDIST.plot.Cat.Rzoo.WGS, dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100, dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000))
#Change factor levels
dtaROHSDIST.Cat.Rzoo$PERC = factor(dtaROHSDIST.Cat.Rzoo$PERC, levels = c("SMALL_ARRAY","LARGE_ARRAY","WGS","TRUE_IBD_100", "TRUE_IBD_1000"))
#Order the df
dtaROHSDIST.Cat.Rzoo = dtaROHSDIST.Cat.Rzoo[order(dtaROHSDIST.Cat.Rzoo$CLASS, dtaROHSDIST.Cat.Rzoo$PERC),]
#RM OTHER DF
rm(dtaROHSDIST.plot.Cat.Rzoo.ARRAY,dtaROHSDIST.plot.Cat.Rzoo.WGS,dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD100,dtaROHSDIST.plot.Cat.Rzoo.TRUEIBD1000)
#Divide TRUE IBD by Cat (because in basepair while Rzoo in kb)
dtaROHSDIST.Cat.Rzoo$x[dtaROHSDIST.Cat.Rzoo$PERC %in% c("TRUE_IBD_100", "TRUE_IBD_1000")] = as.numeric(as.character(dtaROHSDIST.Cat.Rzoo$x[dtaROHSDIST.Cat.Rzoo$PERC %in% c("TRUE_IBD_100", "TRUE_IBD_1000")]))/1000
dtaROHSDIST.Cat.Rzoo$x = as.numeric(as.character(dtaROHSDIST.Cat.Rzoo$x))

## TFPN

## 1K PLINK

## Read the data TFPN Rates Ne 1'000 PLINK 100KB ARRAY ##

PercSeq.ARRAY = c("SMALL\nARRAY", "LARGE\nARRAY")

#Read the data
dtaTFPN.Cat.100KB.50k = read.table("./CattlePop/Analyses/SmallArray_PLINK_TFPN_Rates_GEN100.txt", header = T)
dtaTFPN.Cat.100KB.700k = read.table("./CattlePop/Analyses/LargeArray_PLINK_TFPN_Rates_GEN100.txt", header = T)
dtaTFPN.Cat.100KB.ARRAY = rbind(dtaTFPN.Cat.100KB.50k, dtaTFPN.Cat.100KB.700k)
#Get a vector with mean Rates per RAD_Frag
TP.Cat.100KB.ARRAY = aggregate(dtaTFPN.Cat.100KB.ARRAY$TP, by = list(dtaTFPN.Cat.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
TN.Cat.100KB.ARRAY = aggregate(dtaTFPN.Cat.100KB.ARRAY$TN, by = list(dtaTFPN.Cat.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
FP.Cat.100KB.ARRAY = aggregate(dtaTFPN.Cat.100KB.ARRAY$FP, by = list(dtaTFPN.Cat.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
FN.Cat.100KB.ARRAY = aggregate(dtaTFPN.Cat.100KB.ARRAY$FN, by = list(dtaTFPN.Cat.100KB.ARRAY$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one df
Rates.Cat.100KB.ARRAY = rbind(cbind(TP.Cat.100KB.ARRAY, Rate=rep("TP", nrow(TP.Cat.100KB.ARRAY))), cbind(TN.Cat.100KB.ARRAY, Rate=rep("TN", nrow(TN.Cat.100KB.ARRAY))),
                              cbind(FP.Cat.100KB.ARRAY, Rate=rep("FP", nrow(FP.Cat.100KB.ARRAY))), cbind(FN.Cat.100KB.ARRAY, Rate=rep("FN", nrow(FN.Cat.100KB.ARRAY))))
#Change names
colnames(Rates.Cat.100KB.ARRAY) = c("PercSeq", "Value", "Rate")
#Create the datafarme
stackedRates.Cat.100KB.ARRAY = t(data.frame(Rates.Cat.100KB.ARRAY$Value[Rates.Cat.100KB.ARRAY$Rate == "TN"], Rates.Cat.100KB.ARRAY$Value[Rates.Cat.100KB.ARRAY$Rate == "TP"],
                                            Rates.Cat.100KB.ARRAY$Value[Rates.Cat.100KB.ARRAY$Rate == "FN"], Rates.Cat.100KB.ARRAY$Value[Rates.Cat.100KB.ARRAY$Rate == "FP"]))
row.names(stackedRates.Cat.100KB.ARRAY) = c("TN", "TP", "FN", "FP")
stackedRates.Cat.100KB.ARRAY = as.matrix(as.data.frame(cbind(stackedRates.Cat.100KB.ARRAY[,2],stackedRates.Cat.100KB.ARRAY[,1])))
colnames(stackedRates.Cat.100KB.ARRAY) = PercSeq.ARRAY

## Read the data TFPN Rates Ne 1'000 PLINK 100KB WGS ##

PercSeq.WGS = "WGS"

#Read the data
dtaTFPN.Cat.100KB.WGS = read.table("./CattlePop/Analyses/WGS_PLINK_TFPN_Rates_GEN100.txt", header = T)
#Get a vector with mean Rates per RAD_Frag
TP.Cat.100KB.WGS = aggregate(dtaTFPN.Cat.100KB.WGS$TP, by = list(dtaTFPN.Cat.100KB.WGS$NB_RAD_FRAG), FUN = mean)
TN.Cat.100KB.WGS = aggregate(dtaTFPN.Cat.100KB.WGS$TN, by = list(dtaTFPN.Cat.100KB.WGS$NB_RAD_FRAG), FUN = mean)
FP.Cat.100KB.WGS = aggregate(dtaTFPN.Cat.100KB.WGS$FP, by = list(dtaTFPN.Cat.100KB.WGS$NB_RAD_FRAG), FUN = mean)
FN.Cat.100KB.WGS = aggregate(dtaTFPN.Cat.100KB.WGS$FN, by = list(dtaTFPN.Cat.100KB.WGS$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one dfWGS
Rates.Cat.100KB.WGS = rbind(cbind(TP.Cat.100KB.WGS, Rate=rep("TP", nrow(TP.Cat.100KB.WGS))), cbind(TN.Cat.100KB.WGS, Rate=rep("TN", nrow(TN.Cat.100KB.WGS))),
                            cbind(FP.Cat.100KB.WGS, Rate=rep("FP", nrow(FP.Cat.100KB.WGS))), cbind(FN.Cat.100KB.WGS, Rate=rep("FN", nrow(FN.Cat.100KB.WGS))))
#Change names
colnames(Rates.Cat.100KB.WGS) = c("PercSeq", "Value", "Rate")
#Create the datafarme
stackedRates.Cat.100KB.WGS = t(data.frame(Rates.Cat.100KB.WGS$Value[Rates.Cat.100KB.WGS$Rate == "TN"], Rates.Cat.100KB.WGS$Value[Rates.Cat.100KB.WGS$Rate == "TP"],
                                          Rates.Cat.100KB.WGS$Value[Rates.Cat.100KB.WGS$Rate == "FN"], Rates.Cat.100KB.WGS$Value[Rates.Cat.100KB.WGS$Rate == "FP"]))
row.names(stackedRates.Cat.100KB.WGS) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.Cat.100KB.WGS) = PercSeq.WGS

#Merge WGS, arrays and RAD
stackedRates.Cat.100KB = cbind(stackedRates.Cat.100KB.ARRAY, stackedRates.Cat.100KB.WGS)

#rm the rest
rm(dtaTFPN.Cat.100KB.ARRAY); rm(dtaTFPN.Cat.100KB.50k); rm(dtaTFPN.Cat.100KB.700k); rm(dtaTFPN.Cat.100KB.WGS)
rm(FN.Cat.100KB.ARRAY);rm(FP.Cat.100KB.ARRAY); rm(TN.Cat.100KB.ARRAY); rm(TP.Cat.100KB.ARRAY)
rm(FN.Cat.100KB.WGS); rm(FP.Cat.100KB.WGS); rm(TN.Cat.100KB.WGS); rm(TP.Cat.100KB.WGS)
rm(Rates.Cat.100KB.ARRAY);rm(Rates.Cat.100KB.WGS)
rm(stackedRates.Cat.100KB.ARRAY);rm(stackedRates.Cat.100KB.WGS)


## 1K RZooRoH

PercSeq.ARRAY = c("SMALL\nARRAY", "LARGE\nARRAY")

#Read the data
dtaTFPN.Cat.Rzoo.50k = read.table("./CattlePop/Analyses/SmallArray_RZooRoH_TFPN_Rates_GEN100.txt", header = T)
dtaTFPN.Cat.Rzoo.50k = dtaTFPN.Cat.Rzoo.50k[rowSums(dtaTFPN.Cat.Rzoo.50k[,5:8]) > 0.999999,]
dtaTFPN.Cat.Rzoo.700k = read.table("./CattlePop/Analyses/LargeArray_RZooRoH_TFPN_Rates_GEN100.txt", header = T)
dtaTFPN.Cat.Rzoo.ARRAY = rbind(dtaTFPN.Cat.Rzoo.50k, dtaTFPN.Cat.Rzoo.700k)
#Get a vector with mean Rates per RAD_Frag
TP.Cat.Rzoo.ARRAY = aggregate(dtaTFPN.Cat.Rzoo.ARRAY$TP, by = list(dtaTFPN.Cat.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
TN.Cat.Rzoo.ARRAY = aggregate(dtaTFPN.Cat.Rzoo.ARRAY$TN, by = list(dtaTFPN.Cat.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
FP.Cat.Rzoo.ARRAY = aggregate(dtaTFPN.Cat.Rzoo.ARRAY$FP, by = list(dtaTFPN.Cat.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
FN.Cat.Rzoo.ARRAY = aggregate(dtaTFPN.Cat.Rzoo.ARRAY$FN, by = list(dtaTFPN.Cat.Rzoo.ARRAY$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one df
Rates.Cat.Rzoo.ARRAY = rbind(cbind(TP.Cat.Rzoo.ARRAY, Rate=rep("TP", nrow(TP.Cat.Rzoo.ARRAY))), cbind(TN.Cat.Rzoo.ARRAY, Rate=rep("TN", nrow(TN.Cat.Rzoo.ARRAY))),
                             cbind(FP.Cat.Rzoo.ARRAY, Rate=rep("FP", nrow(FP.Cat.Rzoo.ARRAY))), cbind(FN.Cat.Rzoo.ARRAY, Rate=rep("FN", nrow(FN.Cat.Rzoo.ARRAY))))
#Change names
colnames(Rates.Cat.Rzoo.ARRAY) = c("PercSeq", "Value", "Rate")
#Create the datafarme
stackedRates.Cat.Rzoo.ARRAY = t(data.frame(Rates.Cat.Rzoo.ARRAY$Value[Rates.Cat.Rzoo.ARRAY$Rate == "TN"], Rates.Cat.Rzoo.ARRAY$Value[Rates.Cat.Rzoo.ARRAY$Rate == "TP"],
                                           Rates.Cat.Rzoo.ARRAY$Value[Rates.Cat.Rzoo.ARRAY$Rate == "FN"], Rates.Cat.Rzoo.ARRAY$Value[Rates.Cat.Rzoo.ARRAY$Rate == "FP"]))
row.names(stackedRates.Cat.Rzoo.ARRAY) = c("TN", "TP", "FN", "FP")
stackedRates.Cat.Rzoo.ARRAY = as.matrix(as.data.frame(cbind(stackedRates.Cat.Rzoo.ARRAY[,2],stackedRates.Cat.Rzoo.ARRAY[,1])))
colnames(stackedRates.Cat.Rzoo.ARRAY) = PercSeq.ARRAY

## Read the data TFPN Rates Ne 1'000 PLINK Rzoo WGS ##

PercSeq.WGS = "WGS"

#Read the data
dtaTFPN.Cat.Rzoo.WGS = read.table("./CattlePop/Analyses/WGS_RZooRoH_TFPN_Rates_GEN100.txt", header = T)
#Get a vector with mean Rates per RAD_Frag
TP.Cat.Rzoo.WGS = aggregate(dtaTFPN.Cat.Rzoo.WGS$TP, by = list(dtaTFPN.Cat.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
TN.Cat.Rzoo.WGS = aggregate(dtaTFPN.Cat.Rzoo.WGS$TN, by = list(dtaTFPN.Cat.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
FP.Cat.Rzoo.WGS = aggregate(dtaTFPN.Cat.Rzoo.WGS$FP, by = list(dtaTFPN.Cat.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
FN.Cat.Rzoo.WGS = aggregate(dtaTFPN.Cat.Rzoo.WGS$FN, by = list(dtaTFPN.Cat.Rzoo.WGS$NB_RAD_FRAG), FUN = mean)
#Fuse all rates in one dfWGS
Rates.Cat.Rzoo.WGS = rbind(cbind(TP.Cat.Rzoo.WGS, Rate=rep("TP", nrow(TP.Cat.Rzoo.WGS))), cbind(TN.Cat.Rzoo.WGS, Rate=rep("TN", nrow(TN.Cat.Rzoo.WGS))),
                           cbind(FP.Cat.Rzoo.WGS, Rate=rep("FP", nrow(FP.Cat.Rzoo.WGS))), cbind(FN.Cat.Rzoo.WGS, Rate=rep("FN", nrow(FN.Cat.Rzoo.WGS))))
#Change names
colnames(Rates.Cat.Rzoo.WGS) = c("PercSeq", "Value", "Rate")
#Create the datafarme
stackedRates.Cat.Rzoo.WGS = t(data.frame(Rates.Cat.Rzoo.WGS$Value[Rates.Cat.Rzoo.WGS$Rate == "TN"], Rates.Cat.Rzoo.WGS$Value[Rates.Cat.Rzoo.WGS$Rate == "TP"],
                                         Rates.Cat.Rzoo.WGS$Value[Rates.Cat.Rzoo.WGS$Rate == "FN"], Rates.Cat.Rzoo.WGS$Value[Rates.Cat.Rzoo.WGS$Rate == "FP"]))
row.names(stackedRates.Cat.Rzoo.WGS) = c("TN", "TP", "FN", "FP")
colnames(stackedRates.Cat.Rzoo.WGS) = PercSeq.WGS

#Merge WGS, arrays and RAD
stackedRates.Cat.Rzoo = cbind(stackedRates.Cat.Rzoo.ARRAY, stackedRates.Cat.Rzoo.WGS)

#rm the rest
rm(dtaTFPN.Cat.Rzoo.ARRAY); rm(dtaTFPN.Cat.Rzoo.50k); rm(dtaTFPN.Cat.Rzoo.700k); rm(dtaTFPN.Cat.Rzoo.WGS)
rm(FN.Cat.Rzoo.ARRAY);rm(FP.Cat.Rzoo.ARRAY); rm(TN.Cat.Rzoo.ARRAY); rm(TP.Cat.Rzoo.ARRAY)
rm(FN.Cat.Rzoo.WGS); rm(FP.Cat.Rzoo.WGS); rm(TN.Cat.Rzoo.WGS); rm(TP.Cat.Rzoo.WGS)
rm(Rates.Cat.Rzoo.ARRAY);rm(Rates.Cat.Rzoo.WGS)
rm(stackedRates.Cat.Rzoo.ARRAY);rm(stackedRates.Cat.Rzoo.WGS)

#Define ROHs classes
classesROH = c("< 2", "2 - 4","4 - 6","6 - 10","10 - 16", "> 16")

#PLOT
pdf("~/PhD/ROH/Which_data/Manuscript/SUBMITTED_VERSION/Review_round_I/FIGURES/SUPPLEMENTARY_FIGURES/FIGURES10.pdf", width = 80, height = 60)

#Layout for having all plots
layout(matrix(c(rep(0,9),0,0,0,2,0,0,0,5,0,0,1,0,2,0,4,0,5,0,rep(0,9),0,3,3,3,0,6,6,6,0,rep(0,9)), 6, 9, byrow = T), heights = c(1,.5,1,1,2,.5), widths = c(.7,1.5,1.5,1,2.5,1.5,1.5,1,.5))

#Par for margin and outer margings
par(mar = c(3.5,3.5,3.5,3.5), oma = c(.2,10,.2,.2))

#### PANNEL A: FROH CATTLEs PLINK 100KB ARRAYs ####

#The plot Cattle 50K PLINK100KB
plot(dtaFROH.Cat.100KB.ARRAY$FROH_sub ~ dtaFROH.Cat.100KB.ARRAY$Froh_TRUE_IBD_100GEN, type = 'n', ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n',
     col = "black", pch = 19, xlim = c(0,0.90), ylim = c(0,0.90), cex.axis = 4, frame = F)

#Add points SMALL ARRAY
points(dtaFROH.Cat.100KB.ARRAY$FROH_sub[dtaFROH.Cat.100KB.ARRAY$ARRAY == "SMALL"] ~ dtaFROH.Cat.100KB.ARRAY$Froh_TRUE_IBD_100GEN[dtaFROH.Cat.100KB.ARRAY$ARRAY == "SMALL"],
       pch = 19, col = "#3DED97")
#Add points LARGE ARRAY
points(dtaFROH.Cat.100KB.ARRAY$FROH_sub[dtaFROH.Cat.100KB.ARRAY$ARRAY == "LARGE"] ~ dtaFROH.Cat.100KB.ARRAY$Froh_TRUE_IBD_100GEN[dtaFROH.Cat.100KB.ARRAY$ARRAY == "LARGE"],
       pch = 17, col = "#028A0F")
#Add points WGS
points(dtaFROH.Cat.100KB.WGS$Froh_WGS ~ dtaFROH.Cat.100KB.WGS$Froh_TRUE_IBD_100GEN,
       pch = 17, col = "grey18")

axis(1,xlab = NULL, at=c(0:9)*0.1, padj = 1.5, cex.axis= 8, lwd.ticks = 3, las = 1, tck = -0.03)
axis(2,xlab = NULL, at=c(0:9)*0.1, hadj = 1.5, cex.axis= 8, lwd.ticks = 3, las = 1, tck = -0.03)

#Add the one-to-one line
abline(0,1, lwd = 4)

#Add axes titles
mtext(expression(F[IBD]), side = 1, outer = FALSE, line = 25, cex = 8)
mtext(expression(F[ROH]), side = 2, outer = FALSE, line = 25, cex = 8)
#Add panel lab
mtext(text = "A", side = 3, cex = 8, line = 50, at = -.66)
mtext(text = "PLINK", side = 3, cex = 10, line = 100, at = 1.2)

#### PANNEL B: TFPN CATTLEs PLINK 100KB ARRAYs ####

#Barplot
barplot(stackedRates.Cat.100KB, col = c("#17A398", "#78E33B", "firebrick4", "#F55014"), space = c(rep(0.3,2),1.5), yaxt = 'n',
        xaxt = 'n', axes = F, ylab = NA, xlab = NA, ylim = c(0,1))

#Add the X axis
axis(1, at = c(seq(0.8, to = (2*(1+0.3)-0.5), length.out = 2), (2*(1+0.3)-0.5 + 2.5)), labels = c(PercSeq.ARRAY, PercSeq.WGS),
     cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add thre y axis
axis(2, at = seq(0,10,1)/10, labels = seq(0,10,1)/10, cex.axis = 8, hadj = 1.5, las = 1, lwd.ticks = 5, tck = -0.03)

mtext("Fraction of genome", side = 2, outer = FALSE, line = 20, cex = 4)
#Add panel
mtext(text = "B", side = 3, cex = 8, line = 10, at = -3)

#### PANNEL C: ROHsDist CATTLEs PLINK ####

#Define colors for subsampled dataset
subCOLORS = c("#3DED97", "#028A0F", "grey18")

#Create empty plot
plot(dtaROHSDIST.Cat.PLINK$CLASS, dtaROHSDIST.Cat.PLINK$x/1000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0, 1000))

#Loop through classes to draw boxplots
for(class in as.numeric(sort(unique(dtaROHSDIST.Cat.PLINK$CLASS)))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.30,class-0.10, class+0.10,class+0.30)
  
  #Loop through reduced representation
  for(reduced in 1:3){
    
    #subset the data
    dta_sub = dtaROHSDIST.Cat.PLINK[as.numeric(dtaROHSDIST.Cat.PLINK$CLASS) == class & dtaROHSDIST.Cat.PLINK$PERC == levels(dtaROHSDIST.Cat.PLINK$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.Cat.PLINK$x/1000)[as.numeric(dtaROHSDIST.Cat.PLINK$CLASS) == class & dtaROHSDIST.Cat.PLINK$PERC == "TRUE_IBD_100"],
         ytop = dta_sub$x/1000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.09), x1 = (polyXcoord[reduced] + 0.09), y0 = ((dta_sub$x - dta_sub$sd)/1000), y1 = ((dta_sub$x + dta_sub$sd)/1000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.Cat.PLINK$x/1000)[dtaROHSDIST.Cat.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.Cat.PLINK$CLASS == class],
             y1 = (dtaROHSDIST.Cat.PLINK$x/1000)[dtaROHSDIST.Cat.PLINK$PERC == "TRUE_IBD_100" & dtaROHSDIST.Cat.PLINK$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 8, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add the y axis
axis(2, cex.axis = 4, hadj = 1.5, las = 1, cex.axis = 8, lwd.ticks = 5, tck = -0.03)
#Add y axis title
mtext(text = "Mean Mb ± SD (among individuals)", side = 2, cex = 5, line = 30)
#Add x axis title
mtext(text = "ROHs Length Classes [Mb]", side = 1, cex = 5, line = 25)
#Add panel
mtext(text = "C", side = 3, cex = 8, line = 15, at = -1.2)

#### PANNEL D: FROH CATTLE RZooRoH ####

#The plot Cattle 50K PLINK100KB
plot(dtaFROH.Cat.Rzoo.ARRAY$Froh_array ~ dtaFROH.Cat.Rzoo.ARRAY$Froh_trueIBD100Gen, type = 'n', ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n',
     col = "black", pch = 19, xlim = c(0,0.90), ylim = c(0,0.90), cex.axis = 4, frame = F)

#Add points SMALL ARRAY
points(dtaFROH.Cat.Rzoo.ARRAY$Froh_array[dtaFROH.Cat.Rzoo.ARRAY$ARRAY == "SMALL"] ~ dtaFROH.Cat.Rzoo.ARRAY$Froh_trueIBD100Gen[dtaFROH.Cat.Rzoo.ARRAY$ARRAY == "SMALL"],
       pch = 19, col = "#3DED97")
#Add points LARGE ARRAY
points(dtaFROH.Cat.Rzoo.ARRAY$Froh_array[dtaFROH.Cat.Rzoo.ARRAY$ARRAY == "LARGE"] ~ dtaFROH.Cat.Rzoo.ARRAY$Froh_trueIBD100Gen[dtaFROH.Cat.Rzoo.ARRAY$ARRAY == "LARGE"],
       pch = 17, col = "#028A0F")
#Add points WGS
points(dtaFROH.Cat.RZoo.WGS$FROH_wgs ~ dtaFROH.Cat.RZoo.WGS$Froh_TRUE_IBD_100GEN,
       pch = 17, col = "grey18")

axis(1,xlab = NULL, at=c(0:9)*0.1, padj = 1.5, cex.axis= 8, lwd.ticks = 3, las = 1, tck = -0.03)
axis(2,xlab = NULL, at=c(0:9)*0.1, hadj = 1.5, cex.axis= 8, lwd.ticks = 3, las = 1, tck = -0.03)

#Add the one-to-one line
abline(0,1, lwd = 4)

#Add axes titles
mtext(expression(F[IBD]), side = 1, outer = FALSE, line = 25, cex = 8)
mtext(expression(F[HBD]), side = 2, outer = FALSE, line = 25, cex = 8)
#Add panel lab
mtext(text = "D", side = 3, cex = 8, line = 50, at = -.66)
mtext(text = "RZooRoH", side = 3, cex = 10, line = 100, at = 1.2)

#### PANNEL E: TFPN CATTLEs PLINK 100KB ARRAYs ####

#Barplot
barplot(stackedRates.Cat.Rzoo, col = c("#17A398", "#78E33B", "firebrick4", "#F55014"), space = c(rep(0.3,2),1.5), yaxt = 'n',
        xaxt = 'n', axes = F, ylab = NA, xlab = NA)

#Add the X axis
axis(1, at = c(seq(0.8, to = (2*(1+0.3)-0.5), length.out = 2), (2*(1+0.3)-0.5 + 2.5)), labels = c(PercSeq.ARRAY, PercSeq.WGS),
     cex.axis = 5, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add thre y axis
axis(2, at = seq(0,10,1)/10, labels = seq(0,10,1)/10, cex.axis = 8, hadj = 1.5, las = 1, lwd.ticks = 5, tck = -0.03)

mtext("Fraction of genome", side = 2, outer = FALSE, line = 20, cex = 4)
#Add panel
mtext(text = "E", side = 3, cex = 8, line = 10, at = -3)
#### PANNEL C: ROHsDist CATTLEs PLINK 100KB ARRAYs ####

#Define colors for subsampled dataset
subCOLORS = c("#3DED97", "#028A0F", "grey18")

#Create empty plot
plot(dtaROHSDIST.Cat.Rzoo$CLASS, dtaROHSDIST.Cat.Rzoo$x/1000000, type = 'n', yaxt = 'n', xaxt = 'n', axes = F, ylab = NA, xlab = NA,
     xlim = c(.5,6.5), ylim = c(0, 1000))

#Loop through classes to draw boxplots
for(class in as.numeric(sort(unique(dtaROHSDIST.Cat.Rzoo$CLASS)))){
  
  #Create the polygon (BARPLOTS) X coordinates for different reduced levels
  polyXcoord = c(class-0.30,class-0.10, class+0.10,class+0.30)
  
  #Loop through reduced representation
  for(reduced in 1:3){
    
    #subset the data
    dta_sub = dtaROHSDIST.Cat.Rzoo[as.numeric(dtaROHSDIST.Cat.Rzoo$CLASS) == class & dtaROHSDIST.Cat.Rzoo$PERC == levels(dtaROHSDIST.Cat.Rzoo$PERC)[reduced],]
    
    #draw polygon --> rectangle for this sub in this class
    rect(xleft = polyXcoord[reduced], xright = polyXcoord[reduced + 1],
         ybottom = (dtaROHSDIST.Cat.Rzoo$x/1000)[as.numeric(dtaROHSDIST.Cat.Rzoo$CLASS) == class & dtaROHSDIST.Cat.Rzoo$PERC == "TRUE_IBD_100"],
         ytop = dta_sub$x/1000000, col = subCOLORS[reduced])
    
    #Add error bars
    segments(x0 = (polyXcoord[reduced] + 0.09), x1 = (polyXcoord[reduced] + 0.09), y0 = ((dta_sub$x - dta_sub$sd)/1000), y1 = ((dta_sub$x + dta_sub$sd)/1000000))
    
    #Draw the WGS line (at the end so we see them PERFECTLY)
    segments(y0 = (dtaROHSDIST.Cat.Rzoo$x/1000)[dtaROHSDIST.Cat.Rzoo$PERC == "TRUE_IBD_100" & dtaROHSDIST.Cat.Rzoo$CLASS == class],
             y1 = (dtaROHSDIST.Cat.Rzoo$x/1000)[dtaROHSDIST.Cat.Rzoo$PERC == "TRUE_IBD_100" & dtaROHSDIST.Cat.Rzoo$CLASS == class],
             x0 = (class - 0.4), x1 = (class + 0.4), lwd = 3)
    
  }
  
}

#Add the x axis
axis(1, at = seq(1,6), labels = classesROH, cex.axis = 8, padj = 1.5, lwd.ticks = 5, tck = -0.03)
#Add the y axis
axis(2, hadj = 1.5, las = 1, cex.axis = 8, lwd.ticks = 5, tck = -0.03)
#Add y axis title
mtext(text = "Mean Mb ± SD (among individuals)", side = 2, cex = 5, line = 30)
#Add x axis title
mtext(text = "HBD Length Classes [Mb]", side = 1, cex = 5, line = 25)
#Add panel
mtext(text = "F", side = 3, cex = 8, line = 15, at = -1.2)

dev.off()

############################################
############ FIGURE S11 Cattle FHOM ########
############################################

dtaCatFHOM = read.table("./CattlePop/Analyses/FHOM_WGS_bothArrays.txt", header = T)

#PLOT
pdf("~/PhD/ROH/Which_data/Manuscript/SUBMITTED_VERSION/Review_round_I/FIGURES/SUPPLEMENTARY_FIGURES/FIGURES11.pdf", width = 20, height = 20)

#Par for margin and outer margings
par(mar = c(3.5,3.5,3.5,3.5), oma = c(30,30,2,2))

#### PANNEL A: FROH CATTLEs PLINK 100KB ARRAYs ####

#The plot Cattle 50K PLINK100KB
plot(dtaCatFHOM$FHOM_WGS ~ dtaCatFHOM$Froh_TRUE_IBD_1000GEN, type = 'n', ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n',
     col = "black", pch = 19, xlim = c(0,0.90), ylim = c(-0.3,0.90), cex.axis = 4, frame = F)

#Add points SMALL ARRAY
points(dtaCatFHOM$FHOM_smallArray ~ dtaCatFHOM$Froh_TRUE_IBD_1000GEN, pch = 19, col = "#3DED97")
#Add points LARGE ARRAY
points(dtaCatFHOM$FHOM_largeArray ~ dtaCatFHOM$Froh_TRUE_IBD_1000GEN, pch = 17,col = "#028A0F")
#Add points WGS
points(dtaCatFHOM$FHOM_WGS ~ dtaCatFHOM$Froh_TRUE_IBD_1000GEN, pch = 19, col = "grey18")

axis(1,xlab = NULL, at=c(0:9)*0.1, padj = 1.5, cex.axis= 8, lwd.ticks = 3, las = 1, tck = -0.03)
axis(2,xlab = NULL, at=c(-3:9)*0.1, hadj = 1.5, cex.axis= 8, lwd.ticks = 3, las = 1, tck = -0.03)

#Add the one-to-one line
abline(0,1, lwd = 4)

#Add axes titles
mtext(expression(F[IBD]), side = 1, outer = FALSE, line = 25, cex = 8)
mtext(expression(F[HOM]), side = 2, outer = FALSE, line = 25, cex = 8)

dev.off()

