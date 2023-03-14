################################################
# R compute Phylogenetic Species Diversity Metrics
# by Simonetta Selva
#
# Created: March 9th, 2023
# Project: Bumblebee 2022
################################################
rm(list=ls())

# Preparation ----

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(V.PhyloMaker2)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

setwd(input)

# import phylo tree from (file "Goal 0")
tree <- ape::read.tree("BB22_plant_tree.tre")

# the following code has to be run 4 times with each time different input variables i and j 
# to calculate FD for each bumblebee species and on a site or ID level
# i = "B.lapidarius" or "B.pascuorum"
# j = "site" or "ID.short"
# It can also be performed in a loop, but due to computational time it is easier to them individually. 

# import community matrices (file "compute FD package")
# i <- "B.lapidarius"
# j <- "ID.short"

species <- c("B.lapidarius", "B.pascuorum")
level <- c("site", "ID.short")

for (i in species) {
  for (j in level) {
    sp.ab <- read.delim(paste("./FD/community_matrix_", i, "_", j, ".txt", sep=""), sep = ",", check.names = F)
    colnames(sp.ab) <- sub(" ", "_", colnames(sp.ab))
    
    library(picante)
    PD_var <- psv(sp.ab,tree,compute.var=TRUE,scale.vcv=TRUE)
    PD_ric <- psr(sp.ab,tree,compute.var=TRUE,scale.vcv=TRUE)
    PD_eve <- pse(sp.ab,tree,scale.vcv=TRUE)
    PD_clu <- psc(sp.ab,tree,scale.vcv=TRUE)
    
    # assign(paste("PD", i, j, sep = "_"), 
    #        data.frame(ID = rownames(PD_var),
    #                   nbsp = PD_var$SR,
    #                   PVar = PD_var$PSVs,
    #                   PRic = PD_ric$PSR,
    #                   PEve = PD_eve$PSEs,
    #                   PClu = PD_clu$PSCs))
    
    write.csv(assign(paste("PD", i, j, sep = "_"), 
                     data.frame(ID = rownames(PD_var),
                                nbsp = PD_var$SR,
                                PVar = PD_var$PSVs,
                                PRic = PD_ric$PSR,
                                PEve = PD_eve$PSEs,
                                PClu = PD_clu$PSCs)), 
              file = paste("./PD/PD_", i, "_", j, ".csv", sep = ""))
  }
}


