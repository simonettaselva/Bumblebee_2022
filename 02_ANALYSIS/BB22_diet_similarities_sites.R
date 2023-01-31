################################################
# Consistency of Diet IN a site
# by Simonetta Selva
#
# Created: January, 31th, 2023
# Project: Bumblebee 2022
################################################

rm(list=ls())

# preparation and load data ----

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
BB22.full <- read_csv("BB22_full.csv") %>%
  mutate(site = as_factor(site),
         species = as_factor(species),
         region = substr(site, 1, 3)) 
#remove information on leg and body from ID
BB22_full$ID.short = as.factor(substring(BB22_full$ID,1, nchar(BB22_full$ID)-1)) 


# B.PACUORUM ----
# only use B.pascuorum
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.pascuorum",]

# select a 

# prepare list element to store correlation matrices for species, genus and family
pasc.ID <- list()

# ON FAMILY
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("ID.short", "family")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$family <- rownames(BB22.full.table)














