# Introduction

This repository provides the data and R code associated to Buerki et al. (in prep.). This study aims at providing an updated infra-familial classification of Sapindaceae based on targeted enrichment data (more precisely using the [Angiosperms 353 baiting kit](https://arborbiosci.com/genomics/targeted-sequencing/mybaits/mybaits-expert/mybaits-expert-angiosperms-353/)).

# Contact information

Sven Buerki, Ph.D

Assistant Professor

Department of Biological Sciences

Boise State University

Email: svenbuerki@boisestate.edu

# Project structure

The repository contains:

1. List of all genera of Sapindaceae (`All_genera_Sapindaceae–24oct2020_SB.xlsx`) with data on:
  - Infra-familial classifications: Radlkofer (1933), Buerki et al. (2009, 2011), Acevedo-Rodriguez et al. (2011).
  - Phylogenetic clades inferred by Angiosperms 353 data analyses (clade 22 corresponds to unplaced genera).
  - Updated infra-familial classification (at subfamily and tribal levels).
  - Distribution and species richness (at genus level).
  - Sampling: species sampled for this study and associated vouchers (and where those are deposited).
2. Phylogenetic trees in newick format (in folder `Phylogenetic_trees`):
  - [ASTRAL](https://github.com/smirarab/ASTRAL): `Sapindaceae_astral_dna_species_tree_bs_less10pc_USE_THIS.nwk`
  - [RAxML](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/0471250953.bi0614s51): `Sapindaceae_RAxML_100_bs_dna_species_tree_USE_THIS.nwk`
  - These trees can be opened in R (using e.g. `ape`; see R scripts below) or [FigTree](http://tree.bio.ed.ac.uk/software/figtree/).
3. R scripts to:
  - Produce Table 1 summarizing Sapindaceae infra-familial classification (based on `All_genera_Sapindaceae–24oct2020_SB.csv` and `RAxML_Sapindaeae_25Oct2020_newick.txt`): `Summary_Sapindaceae_classification_Tables_1.R`. The folder `Tables` contains the edited Table in `csv` and `xlsx`.
  - Produce a genus-level Sapindaceae phylogeny based on the RAxML tree (Figure 1): `Phylogenetic_tree_Figure_1.R`. To run this script the R packages [*ape*](https://cran.r-project.org/web/packages/ape/index.html) and [*phytools*](https://cran.r-project.org/web/packages/phytools/index.html) have to be installed on your computer.
  - Produce species-level RAxML and ASTRAL phylogenetic trees (Supplementary Material Figures 1 & 2): `Phylogenetic_trees_Supp_Material_Figures.R`.
  