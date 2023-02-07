################################################
# R compute Functional Diversity Bumblebees
# by Simonetta Selva
#
# Created: December 14, 2022
# Project: Bumblebee 2022
################################################
rm(list=ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

setwd(input)

# load data
BB22.bb.traits <- read_csv("BB22_traits.csv")

# data preparation -----
# replace the names for urban and rural
for (i in 1:nrow(BB22.bb.traits)) {
  if (BB22.bb.traits$landscape[i] == "urban") {
    BB22.bb.traits$landscape[i] = "U"
  } else {
    BB22.bb.traits$landscape[i] = "R"
  }
}

# rename column site to match other dataframes
BB22.bb.traits <- BB22.bb.traits %>%
  mutate(site = paste(location, landscape, replicate, sep = ""),
         region = paste(location, landscape, sep = "")) %>%
  dplyr::select(-NrSpecies, -Shannon)
BB22.bb.traits$site <- as.factor(BB22.bb.traits$site)
BB22.bb.traits$ID <- as.factor(BB22.bb.traits$ID)

# B.pascuorum ----
# select only B.pascuorum
BB22.bb.traits.species <- BB22.bb.traits%>%
  filter(bbspecies == "B.pascuorum")

# prepare functional matrix
traits <- scale(BB22.bb.traits.species[, c(11,13,15)],center=TRUE,scale=TRUE)
rownames(traits)  <- BB22.bb.traits.species$ID
colnames(BB22.bb.traits.species)

# prepare ecological matrix
library(reshape2)
sp.pa <- dcast(BB22.bb.traits.species, site ~ ID, length)[,-1]
rownames(sp.pa)  <- levels(BB22.bb.traits.species$site)
sp.pa[is.na(sp.pa)] <- 0
sp.pa <- as.matrix(sp.pa)

# compute FDs
library(FD)
fd.list <- FD::dbFD(x = traits , a = sp.pa) # not weighted
fd.bb <- data_frame(nbsp = fd.list$nbsp,
                    FRic = fd.list$FRic,
                    FEve = fd.list$FEve,
                    FDiv = fd.list$FDiv)
# relationship bumblebee FD and plant FD
fd.plants <-  read_csv(paste("./FD/FD_package_B.pascuorum_site.csv", sep = ""))









