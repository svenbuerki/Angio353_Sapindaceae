###
# Phylogenetic trees: Produce RAxML and ASTRAL species level phylogenies (= Supp Mat Figures 1 and 2)
###

#Author: Sven Buerki, Boise State University 25th of October 2020

#Load packages
library(ape)
library(phytools)

###~~~
##SUPP MAT FIG 1
###~~~
#Load RAxML tr
ML <- read.tree("Phylogenetic_trees/Sapindaceae_RAxML_100_bs_dna_species_tree_USE_THIS.nwk")
ML$tip.label <- paste(sapply(strsplit(ML$tip.label, split="_"), "[[", 3), sapply(strsplit(ML$tip.label, split="_"), "[[", 4), sep = "_")
ML <- root(ML, outgroup = "Xanthoceras_sorbifolium", resolve.root=T)
ML <- drop.tip(ML, tip = c("Smelophyllum_capense", "Averrhoidium_gardnerianum"))

#Plot tree
pdf("Phylogenetic_trees/RAxML_radial_clades_SuppMatFigure_1.pdf")
#Ladderize tree
x <- ladderize(ML, FALSE)

#Create label ID (phylogenetic clades)
labels <- seq(from=1, to=21, by=1)
#Find MRCA of clades
tips <- ML$tip.label

nodesClades <- c(grep("Xanthoceras", ML$tip.label), 
                 getMRCA(ML, tip=c(which(ML$tip.label == "Acer_campestre"), which(ML$tip.label == "Dipteronia_sinensis"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Billia_hippocastanum"), which(ML$tip.label == "Aesculus_pavia"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Hypelate_trifoliata"), which(ML$tip.label == "Filicium_decipiens"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Magonia_pubescens"), which(ML$tip.label == "Dodonaea_viscosa"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Delavaya_yunnanensis"), which(ML$tip.label == "Ungnadia_speciosa"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Koelreuteria_paniculata"), which(ML$tip.label == "Stocksia_brahuica"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Phyllotrichum_mekongense"), which(ML$tip.label == "Amesiodendron_chinense"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Cubilia_cubili"), which(ML$tip.label == "Euphoria_malaiensis"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Tristira_triptera"), which(ML$tip.label == "Eriocoelum_rubiginosum"))),
                 which(ML$tip.label == "Tristiropsis_NULL"),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Haplocoelum_inopleum"), which(ML$tip.label == "Blighiopsis_pseudostipularis"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Talisia_nervosa"), which(ML$tip.label == "Dilodendron_costaricense"))),
                 which(ML$tip.label == "Blomia_prisca"),
                 which(ML$tip.label == "Guindilia_trinervis"),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Athyana_weinmannifolia"), which(ML$tip.label == "Diatenopteryx_sorbifolia"))),
                 which(ML$tip.label == "Bridgesia_incisifolia"),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Thouinia_acuminata"), which(ML$tip.label == "Allophylus_NULL"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Thinouia_myriantha"), which(ML$tip.label == "Serjania_communis"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Tsingya_bemarana"), which(ML$tip.label == "Beguea_apetala"))),
                 getMRCA(ML, tip=c(which(ML$tip.label == "Alectryon_carinatum"), which(ML$tip.label == "Tina_tamatavensis")))
)
#Plot tree (radial)
plot(x, type="radial", use.edge.length = F, show.node.label = T, cex=.4, label.offset = 0.01)
#Plot nodes with clades IDs
nodes <- labelnodes(text=labels,node=nodesClades,
                    shape="ellipse",cex=0.5,interactive=FALSE)

dev.off()



###~~~
##SUP MAT FIG 2
###~~~
#Load ASTRAL tr
ASTRAL <- read.tree("Phylogenetic_trees/Sapindaceae_astral_dna_species_tree_bs_less10pc_USE_THIS.nwk")
ASTRAL$tip.label <- paste(sapply(strsplit(ASTRAL$tip.label, split="_"), "[[", 3), sapply(strsplit(ASTRAL$tip.label, split="_"), "[[", 4), sep = "_")
ASTRAL <- root(ASTRAL, outgroup = "Xanthoceras_sorbifolium", resolve.root=T)
ASTRAL <- drop.tip(ASTRAL, tip = c("Smelophyllum_capense", "Averrhoidium_gardnerianum"))

#Plot tree
pdf("Phylogenetic_trees/ASTRAL_radial_clades_SuppMatFigure_2.pdf")
#Ladderize tree
x <- ladderize(ASTRAL, FALSE)

#Create label ID (phylogenetic clades)
labels <- seq(from=1, to=21, by=1)
#Find MRCA of clades
tips <- ASTRAL$tip.label

nodesClades <- c(grep("Xanthoceras", ASTRAL$tip.label), 
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Acer_campestre"), which(ASTRAL$tip.label == "Dipteronia_sinensis"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Billia_hippocastanum"), which(ASTRAL$tip.label == "Aesculus_pavia"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Hypelate_trifoliata"), which(ASTRAL$tip.label == "Filicium_decipiens"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Magonia_pubescens"), which(ASTRAL$tip.label == "Dodonaea_viscosa"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Delavaya_yunnanensis"), which(ASTRAL$tip.label == "Ungnadia_speciosa"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Koelreuteria_paniculata"), which(ASTRAL$tip.label == "Stocksia_brahuica"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Phyllotrichum_mekongense"), which(ASTRAL$tip.label == "Amesiodendron_chinense"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Cubilia_cubili"), which(ASTRAL$tip.label == "Euphoria_malaiensis"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Tristira_triptera"), which(ASTRAL$tip.label == "Eriocoelum_rubiginosum"))),
                 which(ASTRAL$tip.label == "Tristiropsis_NULL"),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Haplocoelum_inopleum"), which(ASTRAL$tip.label == "Blighiopsis_pseudostipularis"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Talisia_nervosa"), which(ASTRAL$tip.label == "Dilodendron_costaricense"))),
                 which(ASTRAL$tip.label == "Blomia_prisca"),
                 which(ASTRAL$tip.label == "Guindilia_trinervis"),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Athyana_weinmannifolia"), which(ASTRAL$tip.label == "Diatenopteryx_sorbifolia"))),
                 which(ASTRAL$tip.label == "Bridgesia_incisifolia"),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Thouinia_acuminata"), which(ASTRAL$tip.label == "Allophylus_NULL"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Thinouia_myriantha"), which(ASTRAL$tip.label == "Serjania_communis"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Tsingya_bemarana"), which(ASTRAL$tip.label == "Beguea_apetala"))),
                 getMRCA(ASTRAL, tip=c(which(ASTRAL$tip.label == "Alectryon_carinatum"), which(ASTRAL$tip.label == "Tina_tamatavensis")))
)
#Plot tree (radial)
plot(x, type="radial", use.edge.length = F, show.node.label = T, cex=.4, label.offset = 0.01)
#Plot nodes with clades IDs
nodes <- labelnodes(text=labels,node=nodesClades,
                    shape="ellipse",cex=0.5,interactive=FALSE)

dev.off()

