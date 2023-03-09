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

# import communty matrices (file "compute FD package")
pasc_id <- read.delim("./FD/community_matrix_ B.pascuorum _ ID.short .txt", sep = ",", check.names = F)
colnames(pasc_id) <- sub(" ", "_", colnames(pasc_id))

library(psd)
library(picante)
PD_var <- psv(pasc_id,tree,compute.var=TRUE,scale.vcv=TRUE)
PD_ric <- psr(pasc_id,tree,compute.var=TRUE,scale.vcv=TRUE)
PD_eve <- pse(pasc_id,tree,scale.vcv=TRUE)
PD_clu <- psc(pasc_id,tree,scale.vcv=TRUE)

PD <- data.frame(ID = rownames(PD_var),
                  nbsp = PD_var$SR,
                  PVar = PD_var$PSVs,
                  PRic = PD_ric$PSR,
                  PEve = PD_eve$PSEs,
                  PClu = PD_clu$PSCs)
