rec = read.csv("./12711_2020_593_MOESM2_ESM.csv")

#Create output directory
system("mkdir -p ./RecMaps")

#rm rows where distance = 0
rec = rec[!is.na(rec$Mbp_inter_marker_distance),]
rec = rec[rec$Mbp_inter_marker_distance != 0,]

#RM EXISTING FILES
FILES = list.files(path = "./RecMaps/", full.names = T)
system2("rm", args = FILES)

#Save CHR_length
chr_lengths = c(159000000,137500000,121900000,121300000,121600000,120000000,113200000,114100000,106200000,104700000,107600000,91600000,84600000,85000000,85700000,82000000,75400000,66300000,64300000,72300000,71900000,61700000,52700000,62900000,43100000,51900000,45700000,46500000,51700000)
#Loop throug CHR
for(chr in unique(rec$Chr)){
  rec_sub = rec[rec$Chr == chr,]
  #ORDER by marker position
  rec_sub = rec_sub[order(rec_sub$Mbp_position),]
  #NO header
  #Loop through lines to create new recMAP SLIM FRIENDLY per CHR
  for(row in 1:(nrow(rec_sub) - 1)){
    #If first row, add pos = 0 !
    if(row == 1){
      #creating first line, POS = 0
      writing_row = cbind(0,
            (rec_sub$recrate_adjacent_deterministic[row] * rec_sub$Mbp_inter_marker_distance[row])/1000000)
      #Add the first line
      write.table(writing_row, paste0("./RecMaps/SlimFriendly_RecMap_CHR_", chr),
                  quote = F, col.names = F, row.names = F, append = T, sep = ",")
    #If last row
    } else if(row == (nrow(rec_sub) - 1)){
      #Last line
      writing_row = cbind(as.numeric(chr_lengths[chr])-1,
             (rec_sub$recrate_adjacent_deterministic[row] * rec_sub$Mbp_inter_marker_distance[row])/1000000)
      #Add the line to the RECMAPFILE
      write.table(writing_row, paste0("./RecMaps/SlimFriendly_RecMap_CHR_", chr),
                  quote = F, col.names = F, row.names = F, append = T, sep = ",")
    } else {
      writing_row = cbind(rec_sub$Mbp_position[row+1]*1000000,
            (rec_sub$recrate_adjacent_deterministic[row] * rec_sub$Mbp_inter_marker_distance[row])/1000000)
       #Add the line to the RECMAPFILE
       write.table(writing_row, paste0("./RecMaps/SlimFriendly_RecMap_CHR_", chr),
                   quote = F, col.names = F, row.names = F, append = T, sep = ",")
    }
  }
}




