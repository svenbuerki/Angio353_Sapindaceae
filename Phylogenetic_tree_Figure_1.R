###
# Phylogenetic trees: Produce genus level phylogeny (= Figure 1)
###

#Author: Sven Buerki, Boise State University 18th of January 2021

#Load packages
library(ape)
library(phytools)

#Load RAxML tr
tr <- read.tree("Phylogenetic_trees/Sapindaceae_RAxML_100_bs_dna_species_tree_USE_THIS.nwk")
tr$tip.label <- paste(sapply(strsplit(tr$tip.label, split="_"), "[[", 3), sapply(strsplit(tr$tip.label, split="_"), "[[", 4), sep = "_")
tr <- root(tr, outgroup = "Xanthoceras_sorbifolium", resolve.root=T)
tr <- drop.tip(tr, tip = c("Smelophyllum_capense", "Averrhoidium_gardnerianum", "Distichostemon_dodecarpus","Dimocarpus_longanmaliensis","Tina_coursii"))

# Produce a matrix to select only 1 tip per genus
tips <- tr$tip.label 
tipdrop <- matrix(ncol=4, nrow=length(tips))
colnames(tipdrop) <- c("tipID", "Genus", "Species", "ToKeep")
#Populate matrix
#Tips
tipdrop[,1] <- tips
#Genera
tipdrop[,2] <- sapply(strsplit(tips, split="_"), "[[",3)
#Species
tipdrop[,3] <- sapply(strsplit(tips, split="_"), "[[",4)

#Drop tip to include only 1 accession per genus
#List of genera
genera <- unique(tipdrop[,2])

#Find genus and edit toKeep with 1s
# --> To have a genus level phylogeny
for(i in 1:length(genera)){
  tipdrop[which(tipdrop[,2] == genera[i])[1],4] <- 1
}

#Droptrip
trGenus <- drop.tip(tr, tipdrop[which(is.na(as.vector(tipdrop[,4])) == T),1])

#Rename tips
trGenus$tip.label <- tipdrop[match(trGenus$tip.label, tipdrop[,1]),2]

#Plot tree
pdf("RAxML_genus_level_radial_clades_Figure_1_Jan21.pdf")
#Ladderize tree
x <- ladderize(trGenus, FALSE)

#Create label ID (phylogenetic clades) == TRIBES
labels <- seq(from=1, to=21, by=1)
#Find MRCA of clades
nodesClades <- c(grep("Xanthoceras", trGenus$tip.label), 
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Acer"), which(trGenus$tip.label == "Dipteronia"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Billia"), which(trGenus$tip.label == "Aesculus"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Hypelate"), which(trGenus$tip.label == "Filicium"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Magonia"), which(trGenus$tip.label == "Dodonaea"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Delavaya"), which(trGenus$tip.label == "Ungnadia"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Koelreuteria"), which(trGenus$tip.label == "Stocksia"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Phyllotrichum"), which(trGenus$tip.label == "Amesiodendron"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Cubilia"), which(trGenus$tip.label == "Dimocarpus"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Tristira"), which(trGenus$tip.label == "Eriocoelum"))),
                 which(trGenus$tip.label == "Tristiropsis"),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Haplocoelum"), which(trGenus$tip.label == "Blighiopsis"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Melicoccus"), which(trGenus$tip.label == "Dilodendron"))),
                 which(trGenus$tip.label == "Blomia"),
                 which(trGenus$tip.label == "Guindilia"),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Athyana"), which(trGenus$tip.label == "Diatenopteryx"))),
                 which(trGenus$tip.label == "Bridgesia"),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Thouinia"), which(trGenus$tip.label == "Allophylus"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Thinouia"), which(trGenus$tip.label == "Serjania"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Tsingya"), which(trGenus$tip.label == "Beguea"))),
                 getMRCA(trGenus, tip=c(which(trGenus$tip.label == "Alectryon"), which(trGenus$tip.label == "Tina")))
)

#Create label ID (phylogenetic clades) == Subfamilies
labelsSubfam <- c("I", "II", "III", "IV")

#Plot tree (radial)
plot(x, type="radial", use.edge.length = F, show.node.label = F, cex=.4, label.offset = 0.05)
#Plot nodes with clades IDs
nodes <- labelnodes(text = labels, node = nodesClades,
                    shape = "ellipse", cex = 0.5, interactive = FALSE)
#Plot subfamilies on edges
edgelabels(text = labelsSubfam, edge = c(1,3,13,50), bg='white', frame = "circle")

dev.off()

