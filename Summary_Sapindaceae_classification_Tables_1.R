###~~~
# Produce Table 1 summarizing Sapindaceae infra-generic classification
###~~~

#Author: Sven Buerki, Boise State University 25th of October 2020

#Load data on all genera of Sapindaceae
csv <- read.csv("All_genera_Sapindaceae–24oct2020_SB.csv")


#TABLE 1: The distribution col had to be manually adjusted

#Vector with phylogenetic order (clade 22 = unplaced genera)
PhylOrd <- sort(unique(csv$Phylo_order))

#Create empty matrix for tab1
tab1 <- matrix(ncol=8, nrow=length(PhylOrd))
colnames(tab1) <- c("Clade","This_study", "Buerki_et_al_2009", "Radlkofer", "N_genera", "Genera", "N_species","Distribution")
tab1[,1] <- PhylOrd

#Populate tab1
for(i in 1:length(PhylOrd)){
  #Subset
  sub <- subset(csv, csv$Phylo_order == PhylOrd[i])
  
  #This study
  tab1[i,2] <- unique(as.vector(sub$Buerki.et.al...2020))
  
  #Buerki groupings
  ordgrp <- sort(rowSums(table(as.vector(sub$Group..Buerki.et.al..), as.vector(sub$Genus))), decreasing = T)
  tab1[i,3] <- paste(paste(paste(names(ordgrp)), " (", paste(ordgrp), ")", sep=""), collapse = ", ")
  
  #radlk
  ordgrpRad <- sort(rowSums(table(as.vector(sub$Tribe.Radlk.), as.vector(sub$Genus))), decreasing = T)
  tab1[i,4] <- paste(paste(paste(names(ordgrpRad)), " (", paste(ordgrpRad), ")", sep=""), collapse = ", ")
  
  #N genera
  tab1[i,5] <- nrow(sub)
  
  #List genera (with sp richness)
  tab1[i,6] <- paste(sort(paste0(sub$Genus, " (", sub$N.species, ")")), collapse=", ")
  
  #N species
  tab1[i,7] <- sum(as.numeric(sub$N.species))
  
  #Distr
  tab1[i,8] <- paste(paste(unique(as.vector(sub$Distribution)), sep=', '), collapse = ', ')
}
#Write tab1
write.csv(tab1, file="Tables/Table_1_Sapindaceae_tribes_MS.csv", row.names = F)
