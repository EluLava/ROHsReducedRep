#THIS SCRIPT RUNS ON MY COMPUTER

library(foreach)
library(doParallel)

#Listing the files produced by fregene
files = list.files(path = "./Fregene/Output", pattern = "RecMap_")

#Choosing the nb of cores
registerDoParallel(cores=5)

foreach(i=1:length(files)) %dopar% {

  #Read TOTAL Fregene df
  dta = read.table(files[i], header = T)
  #Subsample the subregions background rates
  subregion_lev = dta[dta$Level == 1,]
  #Nb of subregions
  nb_sub = nrow(subregion_lev)
  #Subsample the hotspots rates
  hotspot_rates = dta[dta$Level == 2,]
  #Subsample the regions rates
  region = dta[dta$Level == 0,]

  #Create the empty df
  rec_rate = as.data.frame(matrix(ncol = 2, nrow = 2))
  colnames(rec_rate) = c("positions", "rates")
  #Fill the first row with position 0 and rate from background rate from subregion 1
  rec_rate[1,] = c(0,subregion_lev$Intensity[1])
  #Fill the second row with position = first row in the hotspot df (because we loop through it from line 2)
  rec_rate[2,] = c(hotspot_rates$Posn_init[1],hotspot_rates$Intensity[1])
  rec_rate$positions = as.integer(rec_rate$positions)

  #Looping through the hotspots df to extract the positions and rates
  for (lin in 2:nrow(hotspot_rates)) {

    #If there is no missing position
    if(hotspot_rates[lin,2] == hotspot_rates[lin - 1,3]){
      rec_rate = rbind(rec_rate, c(hotspot_rates[lin,2],hotspot_rates[lin,4]))
    } else { #If position with bacjground rate is missing
      vec1 = c(hotspot_rates[lin - 1,3], subregion_lev[subregion_lev$Posn_init < hotspot_rates$Posn_final[lin -1] &
                                                       subregion_lev$Posn_final > hotspot_rates$Posn_final[lin - 1],4])
      vec2 = c(hotspot_rates[lin,2], hotspot_rates[lin,4])
      rec_rate = rbind(rec_rate,vec1,vec2)
    }

  }

  #Add the final rate til the last position
  rec_rate = rbind(rec_rate, c((1e8)-1,subregion_lev$Intensity[nrow(subregion_lev)]))

  rec_rate = rec_rate[order(rec_rate$positions),]
  write.table(rec_rate, paste0("./Fregene/Output/Slim_Friendly_", files[i]), row.names = F, quote = F, col.names = F)
}
